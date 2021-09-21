// Computes stress due to fracture and viscoelasticity from deviatoric strain

#pragma once

#include "ComputeStressBase.h"
#include "MathUtils.h"
#include "ExternalPetscSolverApp.h"

class ComputeViscoelasticPFFractureStress;

template <>
InputParameters validParams<ComputeViscoelasticPFFractureStress>();

/**
 * Phase-field fracture
 * This class computes the stress and energy contribution for the
 * small strain Linear Elastic formulation of phase field fracture coupled
 * with Generalized Maxwell viscoelasticity applied only to deviatoric strain
 */
class ComputeViscoelasticPFFractureStress : public ComputeStressBase
{
public:
  ComputeViscoelasticPFFractureStress(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;

  /// Coupled order parameter defining the crack
  const VariableValue & _c;
  const VariableGradient & _grad_c;
  
  /// Material property defining crack width, declared elsewhere
  const MaterialProperty<Real> & _l;

  /// Material property defining gc parameter, declared elsewhere
  const MaterialProperty<Real> & _gc;

  /// Use current value of history variable
  bool _use_current_hist;

  /// History variable that prevents crack healing, declared in this material
  MaterialProperty<Real> & _H;

  /// Old value of history variable
  const MaterialProperty<Real> & _H_old;

  /// material property for fracture energy barrier
  const MaterialProperty<Real> & _barrier;

  /// Material property for elastic energy
  MaterialProperty<Real> & _E;

  /// Derivative of elastic energy w.r.t damage variable
  MaterialProperty<Real> & _dEdc;

  /// Second-order derivative of elastic energy w.r.t damage variable
  MaterialProperty<Real> & _d2Ed2c;

  /// Derivative of stress w.r.t damage variable
  MaterialProperty<RankTwoTensor> & _dstress_dc;

  /// Second-order derivative of elastic energy w.r.t damage variable and strain
  MaterialProperty<RankTwoTensor> & _d2Fdcdstrain;

  /// Material property for energetic degradation function
  const MaterialProperty<Real> & _D;

  /// Derivative of degradation function w.r.t damage variable
  const MaterialProperty<Real> & _dDdc;

  /// Second-order derivative of degradation w.r.t damage variable
  const MaterialProperty<Real> & _d2Dd2c;

  /// Material properties
  const Real _Kv;
  const Real _uinf;
  const Real _ui;
  const Real _Ti;
  const Real _beta;

  /// Strains
  MaterialProperty<RankTwoTensor> & _creep_strain;
  const MaterialProperty<RankTwoTensor> & _creep_strain_old;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;

  /// Viscous dissiption work
  MaterialProperty<Real> & _viscous_dissipation;
  const MaterialProperty<Real> & _viscous_dissipation_old;

  MaterialProperty<Real> & _psi_pos;
  MaterialProperty<Real> & _psi_neg;
  MaterialProperty<Real> & _psi_f;

  virtual void computeQpStress() override;
};
