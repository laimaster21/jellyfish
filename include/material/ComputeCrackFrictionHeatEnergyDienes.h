/// Computes the energy from crack friction

#pragma once

#include "DerivativeMaterialInterface.h"
#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

/// If currentlyComputingJacobian, then the derivative of this quantity wrt total strain

class ComputeCrackFrictionHeatEnergyDienes : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  ComputeCrackFrictionHeatEnergyDienes(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// optional parameter that allows multiple mechanics materials to be defined
  std::string _base_name;

  const VariableValue & _c;
  const VariableValue & _dcdx;
  const VariableValue & _dcdy;

  MaterialProperty<std::vector<Real>> & _crack_normal;
  MaterialProperty<Real> & _crack_normal_norm;

  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankTwoTensor> & _strain_rate;
  const MaterialProperty<RankFourTensor> & _Jacobian_mult;

  const MaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old;

  const Real _friction_coefficient;
  const Real _l;

  MaterialProperty<std::vector<Real>> & _friction_force;
  MaterialProperty<Real> & _friction_normal_force;
  MaterialProperty<std::vector<Real>> & _slide_velocity;
  MaterialProperty<Real> & _slide_velocity_parallel;

  MaterialProperty<Real> & _crack_surface_density;
  MaterialProperty<Real> & _heat_source_rate;

  MaterialProperty<RankTwoTensor> & _velocity_gradient;

};
