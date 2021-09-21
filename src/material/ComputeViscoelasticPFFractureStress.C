// Computes stress due to fracture and viscoelasticity from deviatoric strain

#include "ComputeViscoelasticPFFractureStress.h"

registerMooseObject("TensorMechanicsApp", ComputeViscoelasticPFFractureStress);

template <> InputParameters validParams<ComputeViscoelasticPFFractureStress>()
{
  InputParameters params = validParams<ComputeStressBase>();
  params.addRequiredCoupledVar("c", "Name of damage variable");
  params.addClassDescription("Computes the stress and free energy derivatives for the phase field fracture model, with small strain");
  params.addParam<bool>("use_current_history_variable", false, "Use the current value of the history variable.");
  params.addParam<MaterialPropertyName>("barrier_energy", "Name of material property for fracture energy barrier.");
  params.addParam<MaterialPropertyName>("E_name", "elastic_energy", "Name of material property for elastic energy");
  params.addParam<MaterialPropertyName>("D_name", "degradation", "Name of material property for energetic degradation function.");
  params.addParam<MaterialPropertyName>("F_name", "local_fracture_energy", "Name of material property for local fracture energy function.");
  params.addRequiredParam<Real>("bulk_modulus", "bulk modulus of the incompressible material");
  params.addRequiredParam<Real>("infinity_shear_modulus", "Long term shear modulus of the material");
  params.addRequiredParam<Real>("relaxation_modulus", "shear modulus of the spring in the maxwell unit");
  params.addRequiredParam<Real>("relaxation_time", "characteristic time of the dashpot in the maxwell unit");
  params.addParam<Real>("beta", 0.0, "portion of viscous energy that will contribute to fracture");
  return params;
}

ComputeViscoelasticPFFractureStress::ComputeViscoelasticPFFractureStress(const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _c(coupledValue("c")),
    _grad_c(coupledGradient("c")),
    _l(getMaterialProperty<Real>("l")),
    _gc(getMaterialProperty<Real>("gc_prop")),
    _use_current_hist(getParam<bool>("use_current_history_variable")),
    _H(declareProperty<Real>("hist")),
    _H_old(getMaterialPropertyOld<Real>("hist")),
    _barrier(getDefaultMaterialProperty<Real>("barrier_energy")),
    _E(declareProperty<Real>(getParam<MaterialPropertyName>("E_name"))),
    _dEdc(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("E_name"), getVar("c", 0)->name())),
    _d2Ed2c(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("E_name"), getVar("c", 0)->name(), getVar("c", 0)->name())),
    _dstress_dc(declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name())),
    _d2Fdcdstrain(declareProperty<RankTwoTensor>("d2Fdcdstrain")),
    _D(getMaterialProperty<Real>("D_name")),
    _dDdc(getMaterialPropertyDerivative<Real>("D_name", getVar("c", 0)->name())),
    _d2Dd2c(getMaterialPropertyDerivative<Real>("D_name", getVar("c", 0)->name(), getVar("c", 0)->name())),
    _Kv(getParam<Real>("bulk_modulus")),
    _uinf(getParam<Real>("infinity_shear_modulus")),
    _ui(getParam<Real>("relaxation_modulus")),
    _Ti(getParam<Real>("relaxation_time")),
    _beta(getParam<Real>("beta")),
    _creep_strain(declareProperty<RankTwoTensor>(isParamValid("base_name") ? _base_name + "_creep_strain" : "creep_strain")),
    _creep_strain_old(getMaterialPropertyOld<RankTwoTensor>(isParamValid("base_name") ? _base_name + "_creep_strain" : "creep_strain")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _viscous_dissipation(declareProperty<Real>("viscous_dissipation")),
    _viscous_dissipation_old(getMaterialPropertyOldByName<Real>("viscous_dissipation")),
    _psi_pos(declareProperty<Real>("psi_pos")),
    _psi_neg(declareProperty<Real>("psi_neg")),
    _psi_f(declareProperty<Real>("psi_f"))
{
}

void
ComputeViscoelasticPFFractureStress::initQpStatefulProperties()
{
  _H[_qp] = 0.0;
  const double B=1E10;
  const double PI=3.141592653589793238463;
  if ( std::abs(_q_point[_qp](1)-250.0) < 50.5 ) {
    if ( std::abs(_q_point[_qp](0)-250.0) < (0.5*_l[_qp]) ) {
      _H[_qp] = B*(_gc[_qp]/4/(0.5*_l[_qp]))* (1 - std::abs(_q_point[_qp](0)-250.0)/(0.5*_l[_qp]) ); } }

  _creep_strain[_qp].zero();
  _viscous_dissipation[_qp] = 0.0;
}

void
ComputeViscoelasticPFFractureStress::computeQpStress()
{
  // Distribute mechanical strain into two parts
  RankTwoTensor straindev, straindev_old;
  straindev = _mechanical_strain[_qp].deviatoric();
  straindev_old = _mechanical_strain_old[_qp].deviatoric();

  // Calculate creep strain in the dashpot
  Real hvis = 0.5 * ( 1.0 - exp(-_dt/_Ti) );
  _creep_strain[_qp] = exp(-_dt/_Ti) * _creep_strain_old[_qp] + hvis * (straindev + straindev_old);

  // Assign value for elastic strain
  _elastic_strain[_qp] = _mechanical_strain[_qp] - _creep_strain[_qp];

  // Calculating elastic energy components
  Real straintr, straintr_neg, straintr_pos, we_pos, we_neg;
  RankTwoTensor straindev_elastic;
  straintr = _mechanical_strain[_qp].trace();
  straintr_neg = std::min(straintr, 0.0);
  straintr_pos = straintr - straintr_neg;
  straindev_elastic = straindev - _creep_strain[_qp];
  we_pos = 0.5 * _Kv * 3.0 * std::pow(straintr_pos,2.0) + _uinf * (straindev).doubleContraction(straindev) + _ui * (straindev_elastic).doubleContraction(straindev_elastic);
  we_neg = 0.5 * _Kv * 3.0 * std::pow(straintr_neg,2.0);
  _psi_pos[_qp] = we_pos;
  _psi_neg[_qp] = we_neg;

  // Calculating viscous energy components
  RankTwoTensor si, creep_increment;
  Real psi_v;
  si = 2.0 * _ui * straindev_elastic;
  creep_increment = _creep_strain[_qp] - _creep_strain_old[_qp];
  _viscous_dissipation[_qp] = _viscous_dissipation_old[_qp] + 0.5 * (si).doubleContraction(creep_increment);
  psi_v = _beta *  _viscous_dissipation[_qp];

  // Stress calculation
  RankTwoTensor stress_pos, stress_neg;
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  stress_neg = _Kv * straintr_neg * I2;
  stress_pos = _Kv * straintr_pos * I2 + 2.0 * _uinf * straindev + 2.0 * _ui * straindev_elastic;
  _stress[_qp] = stress_pos * _D[_qp] + stress_neg;

  // Energy split
  Real F_pos, F_neg;
  F_pos = we_pos + psi_v;
  F_neg = we_neg;

  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = stress_pos * _dDdc[_qp];

  // Used in StressDivergencePFFracTensors off-diagonal Jacobian
  _dstress_dc[_qp] = stress_pos * _dDdc[_qp];

  // Jacobian calculation - dstress_dstrain
  RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);
  RankFourTensor Ivol, Idev, Jacobian_neg, Jacobian_pos;
  Ivol = I2.outerProduct(I2) / 3.0;
  Idev = I4sym - Ivol;
  if (straintr >= 0)
    _Jacobian_mult[_qp] = _D[_qp] * (_Kv * Ivol + 2.0 * _uinf * I4sym + 2.0 * _ui * (1.0 - hvis) * Idev);
  else
    _Jacobian_mult[_qp] = _Kv * Ivol + _D[_qp] * (2.0 * _uinf * I4sym + 2.0 * _ui * (1.0 - hvis) * Idev);

  // // Assign history variable
  if (F_pos > _H_old[_qp])
    _H[_qp] = F_pos;
  else
    _H[_qp] = _H_old[_qp];

  Real hist_variable = _H_old[_qp];
  if (_use_current_hist)
    hist_variable = _H[_qp];

  if (hist_variable < _barrier[_qp])
    hist_variable = _barrier[_qp];

  // Elastic free energy density
  _E[_qp] = (hist_variable - psi_v) * _D[_qp] + F_neg;
  _dEdc[_qp] = (hist_variable - psi_v) * _dDdc[_qp];
  _d2Ed2c[_qp] = (hist_variable - psi_v) * _d2Dd2c[_qp];

  _psi_f[_qp] = _gc[_qp] * (std::pow(_c[_qp],2.0) / (4 * _l[_qp]) + _l[_qp] * (_grad_c[_qp] * _grad_c[_qp]));
}
