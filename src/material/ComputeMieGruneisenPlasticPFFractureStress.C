/// Calculates stress for anisortopic crack propagation
/// Includes artificial viscosity and Mie Gruneisen Equation of State
/// Calculated plastic strain based on J2 plasticity

#include "ComputeMieGruneisenPlasticPFFractureStress.h"

registerMooseObject("TensorMechanicsApp", ComputeMieGruneisenPlasticPFFractureStress);

InputParameters
ComputeMieGruneisenPlasticPFFractureStress::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addRequiredCoupledVar("c", "Name of damage variable");
  params.addClassDescription("Computes the stress and free energy derivatives for the phase field fracture model, with small strain"
                             "Considers Mie Gruneisen EOS and artificial viscosity damping");
  params.addParam<bool>("use_current_history_variable", false, "Use the current value of the history variable.");
  params.addParam<bool>("use_snes_vi_solver",false,"Use PETSc's SNES variational inequalities solver to enforce damage "
                        "irreversibility condition and restrict damage value <= 1.");
  params.addParam<MaterialPropertyName>("barrier_energy", "Name of material property for fracture energy barrier.");
  params.addParam<MaterialPropertyName>("E_name", "elastic_energy", "Name of material property for elastic energy");
  params.addParam<MaterialPropertyName>("D_name", "degradation", "Name of material property for energetic degradation function.");
  params.addParam<MaterialPropertyName>("I_name", "indicator", "Name of material property for damage indicator function.");
  params.addParam<MaterialPropertyName>("F_name", "local_fracture_energy", "Name of material property for local fracture energy function.");
  params.addRequiredParam<Real>("Gamma", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  params.addParam<MaterialPropertyName>("density", "density", "Name of Material Property that provides the density");
  params.addParam<MaterialPropertyName>("specific_heat", "specific_heat", "Name of Material Property that provides the density");
  params.addRequiredCoupledVar("temperature","Temperature");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addRequiredParam<Real>("slope_UsUp", "Us-Up slope in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  params.addRequiredParam<Real>("Le","Maximum element size");
  params.addRequiredParam<Real>("sound_speed","Speed of sound in the material");
  params.addRequiredParam<std::vector<Real>>("yield_stress","Input data as pairs of equivalent plastic strain and yield stress: Should start with equivalent plastic strain 0");
  params.addParam<Real>("rtol", 1e-8, "Plastic strain NR tolerance");
  params.addParam<Real>("ftol", 1e-4, "Consistency condition NR tolerance");
  params.addParam<Real>("eptol", 1e-7, "Equivalent plastic strain NR tolerance");
  return params;
}

ComputeMieGruneisenPlasticPFFractureStress::ComputeMieGruneisenPlasticPFFractureStress(const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    _c(coupledValue("c")),
    _l(getMaterialProperty<Real>("l")),
    _pressure(getDefaultMaterialProperty<Real>("fracture_pressure")),
    _gc(getMaterialProperty<Real>("gc_prop")),
    _use_current_hist(getParam<bool>("use_current_history_variable")),
    _use_snes_vi_solver(getParam<bool>("use_snes_vi_solver")),
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
    _I(getDefaultMaterialProperty<Real>("I_name")),
    _dIdc(getMaterialPropertyDerivative<Real>("I_name", getVar("c", 0)->name())),
    _d2Id2c(getMaterialPropertyDerivative<Real>("I_name", getVar("c", 0)->name(), getVar("c", 0)->name())),
    _Gamma(getParam<Real>("Gamma")),
    _density(getMaterialProperty<Real>("density")),
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _temperature(coupledValue("temperature")),
    _ref_temperature(getParam<Real>("reference_temperature")),
    _s(getParam<Real>("slope_UsUp")),
    _strain_increment(getMaterialPropertyByName<RankTwoTensor>(_base_name + "strain_increment")),
    _elastic_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "elastic_strain")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _Le(getParam<Real>("Le")),
    _sound_speed(getParam<Real>("sound_speed")),
    _bulk_modulus(declareProperty<Real>("bulk_modulus")),
    _pressure_eos(declareProperty<Real>("pressure_eos")),
    _yield_stress_vector(getParam<std::vector<Real>>("yield_stress")), // Read from input file
    _plastic_strain(declareProperty<RankTwoTensor>(_base_name + "plastic_strain")),
    _plastic_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "plastic_strain")),
    _eqv_plastic_strain(declareProperty<Real>(_base_name + "eqv_plastic_strain")),
    _eqv_plastic_strain_old(getMaterialPropertyOld<Real>(_base_name + "eqv_plastic_strain")),
    _rotation_increment(getMaterialProperty<RankTwoTensor>(_base_name + "rotation_increment")),
    _rtol(getParam<Real>("rtol")),
    _ftol(getParam<Real>("ftol")),
    _eptol(getParam<Real>("eptol")),
    _deltaOuter(RankTwoTensor::Identity().outerProduct(RankTwoTensor::Identity())),
    _deltaMixed(RankTwoTensor::Identity().mixedProductIkJl(RankTwoTensor::Identity())),
    _W0p(declareProperty<Real>("W0p")),
    _stress_eos_elastic(declareProperty<RankTwoTensor>(_base_name + "stress_eos_elastic")),
    _stress_cpl_elastic(declareProperty<RankTwoTensor>(_base_name + "stress_cpl_elastic")),
    _stress_eos_plastic(declareProperty<RankTwoTensor>(_base_name + "stress_eos_plastic")),
    _stress_cpl_plastic(declareProperty<RankTwoTensor>(_base_name + "stress_cpl_plastic"))
{
}

void
ComputeMieGruneisenPlasticPFFractureStress::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();
  _plastic_strain[_qp].zero();
  _eqv_plastic_strain[_qp] = 0.0;

  _H[_qp] = 0.0;
  const double B=1E10;
  const double PI=3.141592653589793238463;
  // verticle crack
  if ( std::abs(_q_point[_qp](1)-100E-3) < 40E-3 ) {
    if ( std::abs(_q_point[_qp](0)-125E-3) < (0.1*_l[_qp]) ) {
      _H[_qp] = B*(_gc[_qp]/4/(0.5*_l[_qp]))* (1 - std::abs(_q_point[_qp](0)-125E-3)/(0.3*_l[_qp]) ); } }

  // verticle crack 40 nm
  //if ( std::abs(_q_point[_qp](1)-100E-3) < 20E-3 ) {
    //if ( std::abs(_q_point[_qp](0)-125E-3) < (0.1*_l[_qp]) ) {
      //_H[_qp] = B*(_gc[_qp]/4/(0.5*_l[_qp]))* (1 - std::abs(_q_point[_qp](0)-125E-3)/(0.3*_l[_qp]) ); } }

  // verticle crack 60 nm
  //if ( std::abs(_q_point[_qp](1)-100E-3) < 30E-3 ) {
    //if ( std::abs(_q_point[_qp](0)-125E-3) < (0.1*_l[_qp]) ) {
      //_H[_qp] = B*(_gc[_qp]/4/(0.5*_l[_qp]))* (1 - std::abs(_q_point[_qp](0)-125E-3)/(0.3*_l[_qp]) ); } }

// verticle crack lowered
  //if ( std::abs(_q_point[_qp](1)-60E-3) < 40E-3 ) {
    //if ( std::abs(_q_point[_qp](0)-125E-3) < (0.1*_l[_qp]) ) {
      //_H[_qp] = B*(_gc[_qp]/4/(0.5*_l[_qp]))* (1 - std::abs(_q_point[_qp](0)-125E-3)/(0.3*_l[_qp]) ); } }

// horizontal crack
  //if ( std::abs(_q_point[_qp](0)-125E-3) < 40E-3 ) {
    //if ( std::abs(_q_point[_qp](1)-100E-3) < (0.1*_l[_qp]) ) {
      //_H[_qp] = B*(_gc[_qp]/4/(0.5*_l[_qp]))* (1 - std::abs(_q_point[_qp](1)-100E-3)/(0.3*_l[_qp]) ); } }

  // horizontal crack center shift 1 unit (working)
  //if ( std::abs(_q_point[_qp](0)-125E-3) < 40E-3 ) {
    //if ( std::abs(_q_point[_qp](1)-101E-3) < (0.1*_l[_qp]) ) {
      //_H[_qp] = B*(_gc[_qp]/4/(0.5*_l[_qp]))* (1 - std::abs(_q_point[_qp](1)-101E-3)/(0.3*_l[_qp]) ); } }

  // horizontal crack lower for debugging
  //if ( std::abs(_q_point[_qp](0)-125E-3) < 40E-3 ) {
    //if ( std::abs(_q_point[_qp](1)-20E-3) < (0.1*_l[_qp]) ) {
      //_H[_qp] = B*(_gc[_qp]/4/(0.5*_l[_qp]))* (1 - std::abs(_q_point[_qp](1)-20E-3)/(0.3*_l[_qp]) ); } }

// debug log: changed 0.5 to 0.3 in abs denominator, makes no differnece, verified, all switched to 0.3
// since this is a strain term, 1.3 might be due to anisotropy? check 45 deg, other deg to verify this claim
// check 45 deg case, start at 1.3, slowly into 1

//  if ( std::sqrt( std::pow(_q_point[_qp](0)-0.125,2.0) + std::pow(_q_point[_qp](1)-0.080,2.0)) < 0.040 ) {
//    if ( std::sqrt( std::pow(_q_point[_qp](0)-0.125,2.0) + std::pow(_q_point[_qp](1)-0.080,2.0)) > 0.035 ) {
//      if ( _q_point[_qp](1)-0.100 > 0.0 ) {
//        _H[_qp] = B*_gc[_qp]/4/(0.5*_l[_qp]); } } }

  // 45 degree crack narrow(normal width) might still have 1.3 issue, run to verify
  //if ( std::abs(_q_point[_qp](1)+_q_point[_qp](0)-225E-3) < 40E-3 ) {
    //if ( std::abs(_q_point[_qp](1)-_q_point[_qp](0)+25E-3) < (0.1*_l[_qp]) ) {
      //_H[_qp] = B*(_gc[_qp]/4/(0.5*_l[_qp]))* (1 - std::abs(_q_point[_qp](1)-_q_point[_qp](0)+25E-3)/(0.3*_l[_qp]) ); } }

  //if ( std::abs(_q_point[_qp](1)+_q_point[_qp](0)-225E-3) < 40E-3 ) {
    //if ( std::abs(_q_point[_qp](1)-_q_point[_qp](0)+26E-3) < (0.1*_l[_qp]) ) {
      //_H[_qp] = B*(_gc[_qp]/4/(0.5*_l[_qp]))* (1 - std::abs(_q_point[_qp](1)-_q_point[_qp](0)+26E-3)/(0.3*_l[_qp]) ); } }

// test 60 deg case (not wrking)
  //if ( std::abs(_q_point[_qp](1)+_q_point[_qp](0)/std::sqrt(3.0)-(125E-3/std::sqrt(3.0)+100)) < 40E-3 ) {
    //if ( std::abs(_q_point[_qp](1)-std::sqrt(3.0)*_q_point[_qp](0)-100E-3+std::sqrt(3.0)*125E-3) < (0.1*_l[_qp]) ) {
      //_H[_qp] = B*(_gc[_qp]/4/(0.5*_l[_qp]))* (1 - std::abs(_q_point[_qp](1)-std::sqrt(3.0)*_q_point[_qp](0)-100E-3+std::sqrt(3.0)*125E-3)/(0.3*_l[_qp]) ); } }
// working, but location off
  //if ( std::abs(std::sqrt(3.0)*_q_point[_qp](1)+_q_point[_qp](0)-225E-3) < 40E-3 ) {
    //if ( std::abs(_q_point[_qp](1)-std::sqrt(3.0)*_q_point[_qp](0)+25E-3) < (0.1*_l[_qp]) ) {
      //_H[_qp] = B*(_gc[_qp]/4/(0.5*_l[_qp]))* (1 - std::abs(_q_point[_qp](1)-std::sqrt(3.0)*_q_point[_qp](0)+25E-3)/(0.3*_l[_qp]) ); } }
// working 60 deg, draw a picture of the two lines to visulaize
  //if ( std::abs(_q_point[_qp](1)+_q_point[_qp](0)/std::sqrt(3.0)-((125E-3)/std::sqrt(3.0)+100E-3)) < 40E-3 ) {
    //if ( std::abs(_q_point[_qp](1)-std::sqrt(3.0)*_q_point[_qp](0)-100E-3+std::sqrt(3.0)*(125E-3)) < (0.1*_l[_qp]) ) {
      //_H[_qp] = B*(_gc[_qp]/4/(0.5*_l[_qp]))* (1 - std::abs(_q_point[_qp](1)-std::sqrt(3.0)*_q_point[_qp](0)-100E-3+std::sqrt(3.0)*(125E-3))/(0.3*_l[_qp]) ); } }
}

void
ComputeMieGruneisenPlasticPFFractureStress::computeQpStress()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);

  //  Calculate "elastic" stresses due to equation of state
  // Calculate pressure from Mie Gruneisen (Menon, 2014), (Zhang, 2011)
  // https://en.wikipedia.org/wiki/Mie%E2%80%93Gruneisen_equation_of_state
  Real K0, delta, eta, temperature, peos;
  RankTwoTensor stress_eos, stress, stress_cpl;
  K0 = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  _bulk_modulus[_qp] = K0;
  delta = _mechanical_strain[_qp].trace();
  eta = - delta;
  temperature = _temperature[_qp];
  peos = - K0 * eta * (1.0 - (_Gamma * eta / 2.0)) / std::pow((1.0 - _s * eta), 2.0) - _Gamma * _density[_qp] * _specific_heat[_qp] * (temperature - _ref_temperature);
  _pressure_eos[_qp] = peos;
  stress_eos = peos * I2;
  _stress_eos_elastic[_qp] = stress_eos;
  stress_cpl = _elasticity_tensor[_qp] * (_elastic_strain_old[_qp]+_strain_increment[_qp]) - K0 * delta * I2;
  _stress_cpl_elastic[_qp] = stress_cpl;
  stress = stress_eos + stress_cpl;

  //  Calculate plastic strain update
  Real f, flow_incr, rep, dflow_incr, deqvpstrain, fq, yield_stress, err1, err2, err3, toly, eqvpstrain;
  RankTwoTensor flow_dirn, resid, ddsig, delta_dp, df_dsig, plastic_strain;
  RankFourTensor dr_dsig, dr_dsig_inv;

  flow_incr = 0.0;
  dflow_incr = 0.0;
  deqvpstrain = 0.0;
  toly = 1.0e-8;

  unsigned int iter = 0;
  unsigned int maxiter = 100;

  eqvpstrain = _eqv_plastic_strain_old[_qp];
  plastic_strain = _plastic_strain_old[_qp];

  yield_stress = getYieldStress(eqvpstrain); // yield stress at this equivalent plastic strain
  if (yieldFunction(stress, yield_stress) > toly)
  {
    // the sig just calculated is inadmissable.  We must return to the yield surface.
    // This is done iteratively, using a Newton-Raphson process.
    delta_dp.zero();
    flow_dirn = flowPotential(stress);

    resid = flow_dirn * flow_incr - delta_dp; // Residual 1 - refer Hughes Simo
    f = yieldFunction(stress, yield_stress);
    rep = -eqvpstrain + _eqv_plastic_strain_old[_qp] - flow_incr * internalPotential(); // Residual 3 rep=0

    err1 = resid.L2norm();
    err2 = std::abs(f);
    err3 = std::abs(rep);

    while ((err1 > _rtol || err2 > _ftol || err3 > _eptol) &&
           iter < maxiter) // Stress update iteration (hardness fixed)
    {
      iter++;

      df_dsig = dyieldFunction_dstress(stress);
      getJac(stress, _elasticity_tensor[_qp], flow_incr, dr_dsig);   // gets dr_dsig = d(resid_ij)/d(sig_kl)
      fq = dyieldFunction_dinternal(eqvpstrain); // d(f)/d(eqvpstrain)
      dr_dsig_inv = dr_dsig.invSymm();

      dflow_incr = (f - df_dsig.doubleContraction(dr_dsig_inv * resid) + fq * rep) / (df_dsig.doubleContraction(dr_dsig_inv * flow_dirn) - fq);
      ddsig = dr_dsig_inv * (-resid - flow_dirn * dflow_incr); // from solving the top row of linear system, given dflow_incr
      deqvpstrain = rep + dflow_incr; // from solving the bottom row of linear system, given dflow_incr

      // update the variables
      flow_incr += dflow_incr;
      delta_dp -= _elasticity_tensor[_qp].invSymm() * ddsig;
      stress += ddsig;
      eqvpstrain += deqvpstrain;

      // evaluate the RHS equations ready for next Newton-Raphson iteration
      flow_dirn = flowPotential(stress);
      resid = flow_dirn * flow_incr - delta_dp;
      f = yieldFunction(stress, yield_stress);
      rep = -eqvpstrain + _eqv_plastic_strain_old[_qp] - flow_incr * internalPotential();

      err1 = resid.L2norm();
      err2 = std::abs(f);
      err3 = std::abs(rep);
    }

    if (iter >= maxiter)
      mooseError("Constitutive failure");

    plastic_strain += delta_dp;
  }

  // Rotate plastic strain tensor to the current configuration
  _plastic_strain[_qp] = plastic_strain;

  // Calculate the elastic strain_increment
  _elastic_strain[_qp] = _mechanical_strain[_qp] - _plastic_strain[_qp];

  //  Store plasticity updated stress
  peos = stress.trace()/3;
  stress_eos = peos * I2;
  _stress_eos_plastic[_qp] = stress_eos;
  stress_cpl = stress - stress_eos;
  _stress_cpl_plastic[_qp] = stress_cpl;

  // Create the positive and negative projection tensors
  std::vector<Real> eigval;
  RankTwoTensor eigvec;
  RankFourTensor Ppos = stress_cpl.positiveProjectionEigenDecomposition(eigval, eigvec);

  // Project the positive and negative stresses
  Real peos_pos;
  RankTwoTensor stress_eos_pos, stress_eos_neg, stress_cpl_pos, stress_cpl_neg, stress0pos, stress0neg;
  peos_pos = (std::abs(peos) + peos) / 2.0;
  stress_eos_pos = peos_pos * I2;
  stress_eos_neg = stress_eos - stress_eos_pos;
  stress_cpl_pos = Ppos * stress_cpl;
  stress_cpl_neg = stress_cpl - stress_cpl_pos;
  stress0pos = stress_eos_pos + stress_cpl_pos;
  stress0neg = stress - stress0pos;

  // Compute the positive and negative elastic energies
  Real F_pos, F_neg;
  Real A, B, C;
  A = K0 * (_Gamma * (1.0/(2.0 * std::pow(_s,2.0)) + 0.5 - 1.0 / _s) - 1.0);
  B = K0 * (_Gamma * (1.0 / _s - 0.5) + 1.0);
  C = - K0 * _Gamma / (2.0 * std::pow(_s,2.0));
  delta = _mechanical_strain[_qp].trace();
  eta = - _elastic_strain[_qp].trace();
  temperature = _temperature[_qp];
  if (delta>=0.0) {
    F_pos = A / (_s - std::pow(_s,2.0) * eta) - B * std::log(1.0 - _s * eta) / _s + C * eta - A / _s -
            _Gamma * _density[_qp] * _specific_heat[_qp] * (temperature - _ref_temperature) * eta;
    F_neg = 0.0;
  }
  else {
    F_pos = 0.0;
    F_neg = A / (_s - std::pow(_s,2.0) * eta) - B * std::log(1.0 - _s * eta) / _s + C * eta - A / _s -
            _Gamma * _density[_qp] * _specific_heat[_qp] * (temperature - _ref_temperature) * eta;
  }

  F_pos += (stress_cpl_pos).doubleContraction(_elastic_strain[_qp]) / 2.0;
  F_neg += (stress_cpl_neg).doubleContraction(_elastic_strain[_qp]) / 2.0;

  _stress[_qp] = stress0pos * _D[_qp] - _pressure[_qp] * I2 * _I[_qp] + stress0neg;

  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = stress0pos * _dDdc[_qp];

  // Used in StressDivergencePFFracTensors off-diagonal Jacobian
  _dstress_dc[_qp] = stress0pos * _dDdc[_qp] - _pressure[_qp] * I2 * _dIdc[_qp];

  _Jacobian_mult[_qp] = (I4sym - (1 - _D[_qp]) * Ppos) * _elasticity_tensor[_qp];

  // Calculate bulk-viscosity stress term
  Real trD, jacob, q_bv;
  trD = (_mechanical_strain[_qp].trace() - _mechanical_strain_old[_qp].trace()) / _dt;
  jacob = 1.0 + _mechanical_strain[_qp].trace();
  q_bv = 0.0;
  if (jacob < 1.0) {
    q_bv = ( _C0 * _density[_qp] * trD * std::abs(trD) * std::pow(_Le,2.0) / std::pow(jacob,2.0) ) + ( _C1 * _density[_qp] * _sound_speed * trD * _Le / jacob );
  }
  _stress[_qp] += q_bv * I2;


  // // Assign history variable
  Real hist_variable = _H_old[_qp];
  if (_use_snes_vi_solver)
  {
    _H[_qp] = F_pos;

    if (_use_current_hist)
      hist_variable = _H[_qp];
  }
  else
  {
    if (F_pos > _H_old[_qp])
      _H[_qp] = F_pos;
    else
      _H[_qp] = _H_old[_qp];

    if (_use_current_hist)
      hist_variable = _H[_qp];

    if (hist_variable < _barrier[_qp])
      hist_variable = _barrier[_qp];
  }

  // Elastic free energy density
  _E[_qp] = hist_variable * _D[_qp] + F_neg - _pressure[_qp] * _elastic_strain[_qp].trace() * _I[_qp];
  _dEdc[_qp] = hist_variable * _dDdc[_qp] - _pressure[_qp] * _elastic_strain[_qp].trace() * _dIdc[_qp];
  _d2Ed2c[_qp] = hist_variable * _d2Dd2c[_qp] - _pressure[_qp] * _elastic_strain[_qp].trace() * _d2Id2c[_qp];
}

Real
ComputeMieGruneisenPlasticPFFractureStress::yieldFunction(const RankTwoTensor & stress, const Real yield_stress)
{
  return getSigEqv(stress) - yield_stress;
}

RankTwoTensor
ComputeMieGruneisenPlasticPFFractureStress::dyieldFunction_dstress(const RankTwoTensor & sig)
{
  RankTwoTensor deriv = sig.dsecondInvariant();
  deriv *= std::sqrt(3.0 / sig.secondInvariant()) / 2.0;
  return deriv;
}

Real
ComputeMieGruneisenPlasticPFFractureStress::dyieldFunction_dinternal(const Real equivalent_plastic_strain)
{
  return -getdYieldStressdPlasticStrain(equivalent_plastic_strain);
}

RankTwoTensor
ComputeMieGruneisenPlasticPFFractureStress::flowPotential(const RankTwoTensor & sig)
{
  return dyieldFunction_dstress(sig); // this plasticity model assumes associative flow
}

Real
ComputeMieGruneisenPlasticPFFractureStress::internalPotential()
{
  return -1;
}

Real
ComputeMieGruneisenPlasticPFFractureStress::getSigEqv(const RankTwoTensor & stress)
{
  return std::sqrt(3 * stress.secondInvariant());
}

// Jacobian for stress update algorithm
void
ComputeMieGruneisenPlasticPFFractureStress::getJac(const RankTwoTensor & sig,
                                    const RankFourTensor & E_ijkl,
                                    Real flow_incr,
                                    RankFourTensor & dresid_dsig)
{
  RankTwoTensor sig_dev, df_dsig, flow_dirn;
  RankTwoTensor dfi_dft, dfi_dsig;
  RankFourTensor dft_dsig, dfd_dft, dfd_dsig;
  Real sig_eqv;
  Real f1, f2, f3;
  RankFourTensor temp;

  sig_dev = sig.deviatoric();
  sig_eqv = getSigEqv(sig);
  df_dsig = dyieldFunction_dstress(sig);
  flow_dirn = flowPotential(sig);

  f1 = 3.0 / (2.0 * sig_eqv);
  f2 = f1 / 3.0;
  f3 = 9.0 / (4.0 * Utility::pow<3>(sig_eqv));

  dft_dsig = f1 * _deltaMixed - f2 * _deltaOuter - f3 * sig_dev.outerProduct(sig_dev);

  dfd_dsig = dft_dsig;
  dresid_dsig = E_ijkl.invSymm() + dfd_dsig * flow_incr;
}

// Obtain yield stress for a given equivalent plastic strain (input)
Real
ComputeMieGruneisenPlasticPFFractureStress::getYieldStress(const Real eqpe)
{
  unsigned nsize;
  nsize = _yield_stress_vector.size();

  if (_yield_stress_vector[0] > 0.0 || nsize % 2 > 0) // Error check for input inconsitency
    mooseError("Error in yield stress input: Should be a vector with eqv plastic strain and yield "
               "stress pair values.\n");

  unsigned int ind = 0;
  Real tol = 1e-8;
  while (ind < nsize) {
    if (std::abs(eqpe - _yield_stress_vector[ind]) < tol)
      return _yield_stress_vector[ind + 1];

    if (ind + 2 < nsize) {
      if (eqpe > _yield_stress_vector[ind] && eqpe < _yield_stress_vector[ind + 2])
        return _yield_stress_vector[ind + 1] +
               (eqpe - _yield_stress_vector[ind]) /
                   (_yield_stress_vector[ind + 2] - _yield_stress_vector[ind]) *
                   (_yield_stress_vector[ind + 3] - _yield_stress_vector[ind + 1]);
    }
    else
      return _yield_stress_vector[nsize - 1];
    ind += 2;
  }
  return 0.0;
}

Real
ComputeMieGruneisenPlasticPFFractureStress::getdYieldStressdPlasticStrain(const Real eqpe)
{
  unsigned nsize;
  nsize = _yield_stress_vector.size();

  if (_yield_stress_vector[0] > 0.0 || nsize % 2 > 0) // Error check for input inconsitency
    mooseError("Error in yield stress input: Should be a vector with eqv plastic strain and yield "
               "stress pair values.\n");

  unsigned int ind = 0;
  while (ind < nsize) {
    if (ind + 2 < nsize) {
      if (eqpe >= _yield_stress_vector[ind] && eqpe < _yield_stress_vector[ind + 2])
        return (_yield_stress_vector[ind + 3] - _yield_stress_vector[ind + 1]) /
               (_yield_stress_vector[ind + 2] - _yield_stress_vector[ind]);
    }
    else
      return 0.0;
    ind += 2;
  }
  return 0.0;
}
