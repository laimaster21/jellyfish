################################################################################
# Material - HMX
# Units: Length-micrometer, Mass-picogram, Time-nanosecond
# Domain size - 2D - 250 x 300 nm
# Element - Quad of 2 nm
# L0 = 8 nm
# Plasticity - J2
# Yield stress - 0.26 GPa
# Loading velocity - 500 m/s
# Boundary conditions:
#     Left - periodic
#     Bottom - displacement loading
#     Right - periodic
#     Top - fixed
################################################################################

[GlobalParams]
  displacements = 'ux uy'
[]
################################################################################

[Mesh]
  type = GeneratedMesh
  dim = 2
  xmax = 250E-3
  nx = 125
  ymax = 300E-3
  ny = 150
[]
################################################################################

[Variables]
  [./ux]
  [../]
  [./uy]
  [../]
  [./c]
  [../]
  [./temperature]
  [../]
	[./dirac_switch]
	order = CONSTANT
	family = MONOMIAL
	[../]
[]
################################################################################

[AuxVariables]
  [./stressxx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stressyy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stresszz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stressxy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stressyz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stressxz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stresseosplasticxx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stresscplelasticxx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stresscplelasticyy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stresscplelasticzz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stresscplelasticxy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stresscplelasticyz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stresscplelasticxz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stresscplplasticxx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stresscplplasticyy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stresscplplasticzz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stresscplplasticxy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stresscplplasticyz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stresscplplasticxz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./volstress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmisesstress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./totalstrainxx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./totalstrainyy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./totalstrainzz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./totalstrainxy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./totalstrainyz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./totalstrainxz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elasticstrainxx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elasticstrainyy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elasticstrainzz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elasticstrainxy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elasticstrainyz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elasticstrainxz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plasticstrainxx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plasticstrainyy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plasticstrainzz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plasticstrainxy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plasticstrainyz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plasticstrainxz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vx]
  [../]
  [./ax]
  [../]
  [./vy]
  [../]
  [./ay]
  [../]
  [./dcdx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./dcdy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pressure_eos]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pressure]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plastic_heat]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./friction_heat]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./shock_heat]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./expansion_heat]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vis_heat]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./jetting_heat]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]
################################################################################

[Functions]
  [./loading]
    type = ParsedFunction
    value = '0.5*t'
  [../]
  [./dirac_switch_induction_time]
    type = ParsedFunction
    value = '1. / 0.001'
  [../]
[]
################################################################################

[Kernels]
  [./DynamicTensorMechanics]
  [../]
  [./ACbulk]
    type = AllenCahn
    variable = c
    f_name = F
  [../]
  [./ACInterfaceCleavageFracture]
    type = ACInterfaceCleavageFracture
    variable = c
    beta_penalty = 10
    cleavage_plane_normal = '0 1 1'
  [../]
  [./solid_x]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = ux
    component = 0
    c = c
  [../]
  [./solid_y]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = uy
    component = 1
    c = c
  [../]
  [./ACOffDiag]
    type = AllenCahnElasticEnergyOffDiag
    variable = c
  [../]
  [./dcdt]
    type = ADTimeDerivative
    variable = c
  [../]
  [./inertiax]
    type = InertialForce
    variable = ux
    velocity = vx
    acceleration = ax
    beta = 0.3025
    gamma = 0.6
  [../]
  [./inertiay]
    type = InertialForce
    variable = uy
    velocity = vy
    acceleration = ay
    beta = 0.3025
    gamma = 0.6
  [../]
  [./hc]
    type = HeatConduction
    variable = temperature
  [../]
  [./hct]
    type = HeatConductionTimeDerivative
    variable = temperature
  [../]
  [./friction]
    type = CrackFrictionHeatSource
    variable = temperature
    friction_coefficient = 0.0
    dcdx = dcdx
    dcdy = dcdy
  [../]

  [./molecular_jetting]
    type = MolecularJettingHeatSource
    variable = temperature
		c=c
  [../]

  [./thermoelastic_heat_source]
    type = ThermalExpansionHeatSourceSmallStrainMieGruneisen
    variable = temperature
    Gamma = 0.7
    c = c
    specific_heat = specific_heat
    density = density
  [../]
  [./plasticheat]
    type = PlasticHeatEnergy
    coeff = 0.5
    variable = temperature
  [../]
  [./d_dirac_switch_dt]
    type = ADTimeDerivative
    variable = dirac_switch
    use_displaced_mesh = true #false
  [../]
  [./dirac_switch_constant_omega_tau]
    type = MaskedBodyForce
    variable = dirac_switch
    value = 1.0 
    function = dirac_switch_induction_time
    mask = dirac_switch_pressure
    use_displaced_mesh = true #false
  [../]
[]
################################################################################

[AuxKernels]
  [./stressxx]
    type = RankTwoAux
    variable = stressxx
    rank_two_tensor = stress
    index_j = 0
    index_i = 0
  [../]
  [./stressyy]
    type = RankTwoAux
    variable = stressyy
    rank_two_tensor = stress
    index_j = 1
    index_i = 1
  [../]
  [./stresszz]
    type = RankTwoAux
    variable = stresszz
    rank_two_tensor = stress
    index_j = 2
    index_i = 2
  [../]
  [./stressxy]
    type = RankTwoAux
    variable = stressxy
    rank_two_tensor = stress
    index_j = 1
    index_i = 0
  [../]
  [./stressyz]
    type = RankTwoAux
    variable = stressyz
    rank_two_tensor = stress
    index_j = 2
    index_i = 1
  [../]
  [./stressxz]
    type = RankTwoAux
    variable = stressxz
    rank_two_tensor = stress
    index_j = 2
    index_i = 0
  [../]
  [./stresseosplasticxx]
    type = RankTwoAux
    variable = stresseosplasticxx
    rank_two_tensor = stress_eos_plastic
    index_j = 0
    index_i = 0
  [../]
  [./stresscplelasticxx]
    type = RankTwoAux
    variable = stresscplelasticxx
    rank_two_tensor = stress_cpl_elastic
    index_j = 0
    index_i = 0
  [../]
  [./stresscplelasticyy]
    type = RankTwoAux
    variable = stresscplelasticyy
    rank_two_tensor = stress_cpl_elastic
    index_j = 1
    index_i = 1
  [../]
  [./stresscplelasticzz]
    type = RankTwoAux
    variable = stresscplelasticzz
    rank_two_tensor = stress_cpl_elastic
    index_j = 2
    index_i = 2
  [../]
  [./stresscplelasticxy]
    type = RankTwoAux
    variable = stresscplelasticxy
    rank_two_tensor = stress_cpl_elastic
    index_j = 1
    index_i = 0
  [../]
  [./stresscplelasticyz]
    type = RankTwoAux
    variable = stresscplelasticyz
    rank_two_tensor = stress_cpl_elastic
    index_j = 2
    index_i = 1
  [../]
  [./stresscplelasticxz]
    type = RankTwoAux
    variable = stresscplelasticxz
    rank_two_tensor = stress_cpl_elastic
    index_j = 2
    index_i = 0
  [../]
  [./stresscplplasticxx]
    type = RankTwoAux
    variable = stresscplplasticxx
    rank_two_tensor = stress_cpl_plastic
    index_j = 0
    index_i = 0
  [../]
  [./stresscplplasticyy]
    type = RankTwoAux
    variable = stresscplplasticyy
    rank_two_tensor = stress_cpl_plastic
    index_j = 1
    index_i = 1
  [../]
  [./stresscplplasticzz]
    type = RankTwoAux
    variable = stresscplplasticzz
    rank_two_tensor = stress_cpl_plastic
    index_j = 2
    index_i = 2
  [../]
  [./stresscplplasticxy]
    type = RankTwoAux
    variable = stresscplplasticxy
    rank_two_tensor = stress_cpl_plastic
    index_j = 1
    index_i = 0
  [../]
  [./stresscplplasticyz]
    type = RankTwoAux
    variable = stresscplplasticyz
    rank_two_tensor = stress_cpl_plastic
    index_j = 2
    index_i = 1
  [../]
  [./stresscplplasticxz]
    type = RankTwoAux
    variable = stresscplplasticxz
    rank_two_tensor = stress_cpl_plastic
    index_j = 2
    index_i = 0
  [../]
  [./volstress]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = volstress
    scalar_type = Hydrostatic
  [../]
  [./vonmisesstress]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmisesstress
    scalar_type = VonMisesStress
  [../]
  [./totalstrainxx]
    type = RankTwoAux
    variable = totalstrainxx
    rank_two_tensor = total_strain
    index_j = 0
    index_i = 0
  [../]
  [./totalstrainyy]
    type = RankTwoAux
    variable = totalstrainyy
    rank_two_tensor = total_strain
    index_j = 1
    index_i = 1
  [../]
  [./totalstrainzz]
    type = RankTwoAux
    variable = totalstrainzz
    rank_two_tensor = total_strain
    index_j = 2
    index_i = 2
  [../]
  [./totalstrainxy]
    type = RankTwoAux
    variable = totalstrainxy
    rank_two_tensor = total_strain
    index_j = 1
    index_i = 0
  [../]
  [./totalstrainyz]
    type = RankTwoAux
    variable = totalstrainyz
    rank_two_tensor = total_strain
    index_j = 2
    index_i = 1
  [../]
  [./totalstrainxz]
    type = RankTwoAux
    variable = totalstrainxz
    rank_two_tensor = total_strain
    index_j = 2
    index_i = 0
  [../]
  [./elasticstrainxx]
    type = RankTwoAux
    variable = elasticstrainxx
    rank_two_tensor = elastic_strain
    index_j = 0
    index_i = 0
  [../]
  [./elasticstrainyy]
    type = RankTwoAux
    variable = elasticstrainyy
    rank_two_tensor = elastic_strain
    index_j = 1
    index_i = 1
  [../]
  [./elasticstrainzz]
    type = RankTwoAux
    variable = elasticstrainzz
    rank_two_tensor = elastic_strain
    index_j = 2
    index_i = 2
  [../]
  [./elasticstrainxy]
    type = RankTwoAux
    variable = elasticstrainxy
    rank_two_tensor = elastic_strain
    index_j = 1
    index_i = 0
  [../]
  [./elasticstrainyz]
    type = RankTwoAux
    variable = elasticstrainyz
    rank_two_tensor = elastic_strain
    index_j = 2
    index_i = 1
  [../]
  [./elasticstrainxz]
    type = RankTwoAux
    variable = elasticstrainxz
    rank_two_tensor = elastic_strain
    index_j = 2
    index_i = 0
  [../]
  [./plasticstrainxx]
    type = RankTwoAux
    variable = plasticstrainxx
    rank_two_tensor = plastic_strain
    index_j = 0
    index_i = 0
  [../]
  [./plasticstrainyy]
    type = RankTwoAux
    variable = plasticstrainyy
    rank_two_tensor = plastic_strain
    index_j = 1
    index_i = 1
  [../]
  [./plasticstrainzz]
    type = RankTwoAux
    variable = plasticstrainzz
    rank_two_tensor = plastic_strain
    index_j = 2
    index_i = 2
  [../]
  [./plasticstrainxy]
    type = RankTwoAux
    variable = plasticstrainxy
    rank_two_tensor = plastic_strain
    index_j = 1
    index_i = 0
  [../]
  [./plasticstrainyz]
    type = RankTwoAux
    variable = plasticstrainyz
    rank_two_tensor = plastic_strain
    index_j = 2
    index_i = 1
  [../]
  [./plasticstrainxz]
    type = RankTwoAux
    variable = plasticstrainxz
    rank_two_tensor = plastic_strain
    index_j = 2
    index_i = 0
  [../]
  [./ax]
    type = NewmarkAccelAux
    variable = ax
    displacement = ux
    velocity = vx
    beta = 0.3025
  [../]
  [./vx]
    type = NewmarkVelAux
    variable = vx
    acceleration = ax
    gamma = 0.6
  [../]
  [./ay]
    type = NewmarkAccelAux
    variable = ay
    displacement = uy
    velocity = vy
    beta = 0.3025
  [../]
  [./vy]
    type = NewmarkVelAux
    variable = vy
    acceleration = ay
    gamma = 0.6
  [../]
  [./dcdx]
    type = VariableGradientComponent
    variable = dcdx
    gradient_variable = c
    component = 'x'
  [../]
  [./dcdy]
    type = VariableGradientComponent
    variable = dcdy
    gradient_variable = c
    component = 'y'
  [../]
  [./pressure_eos]
    type = MaterialRealAux
    variable = pressure_eos
    property = pressure_eos
  [../]
  [./pressure]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = pressure
    scalar_type = Hydrostatic
  [../]
  [./plastic_heat]
    type = MaterialRealAux
    variable = plastic_heat
    property = plastic_heat
  [../]
  [./friction_heat]
    type = MaterialRealAux
    variable = friction_heat
    property = heat_source_rate
  [../]
  [./shock_heat]
    type = MaterialRealAux
    variable = shock_heat
    property = shock_heat
  [../]
  [./expansion_heat]
    type = MaterialRealAux
    variable = expansion_heat
    property = expansion_heat
  [../]
  [./vis_heat]
    type = MaterialRealAux
    variable = vis_heat
    property = vis_heat
  [../]
  [./jetting_heat]
    type = MaterialRealAux
    variable = jetting_heat
    property = molecular_jetting_heat
  [../]
[]
################################################################################

[ICs]
  [./temperature_ic]
    type = ConstantIC
    variable = temperature
    value = 300.0
  [../]
  [./dirac_switch_ic]
    type = ConstantIC
    variable = dirac_switch
    block = 0
    value = 0.0
  [../]
[]
################################################################################

[BCs]
  [./comp_bottom]
    type = FunctionDirichletBC
    preset = true
    variable = uy
    boundary = bottom
    function = loading
  [../]
  [./topfixed]
    type = DirichletBC
    preset = true
    boundary = top
    variable = uy
    value = 0.0
  [../]
  [./Periodic]
    [./y]
      primary = left
      secondary = right
      translation = '250E-3 0 0'
    [../]
  [../]
  [./TBC]
    type = NeumannBC
    variable = temperature
    boundary = 'right top left bottom'
    value = 0.0
  [../]
[]
################################################################################

[Materials]
  [./mat_prop]
    type = GenericConstantMaterial
    prop_names = 'density gc_prop l visco thermal_conductivity specific_heat'
    prop_values = '1.85 0.002 8E-3 10 0.31e-6 2357.3e-6'
  [../]
  [./mobility]
    type = ParsedMaterial
    material_property_names = 'visco gc_prop'
    f_name = L
    function = '1.0/(visco*gc_prop)'
  [../]
  [./kappa]
    type = ParsedMaterial
    material_property_names = 'gc_prop l'
    f_name = kappa_op
    function = 'gc_prop*l'
  [../]
  [./damage_stress_hmx]
    type = ComputeMieGruneisenPlasticPFFractureStress
    c = c
    E_name = 'elastic_energy'
    D_name = 'degradation'
    I_name = 'indicator'
    F_name = 'local_elastic_energy'
    Gamma = 0.7
    density = density
    specific_heat = specific_heat
    temperature = temperature
    reference_temperature = 300.0
    slope_UsUp = 2.29
    C0 = 0.1
    C1 = 1.0
    Le = 2E-3
    sound_speed = 2.77
    yield_stress='0. 0.26 100. 0.26'
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '25.1 9.7 12.8 0.0 -1.3 0.0 22.3 11.8 0.0 4.6 0.0 21.8 0.0 1.4 0.0 9.7 0.0 3.18 11.036 0.0 8.66'
    fill_method = symmetric21
    euler_angle_1 = 0   #deg
    euler_angle_2 = 0   #deg
    euler_angle_3 = 0   #deg
  [../]
  [./strain]
    type = ComputeIncrementalSmallStrain
  [../]
  [./degradation]
    type = DerivativeParsedMaterial
    f_name = degradation
    args = 'c'
    function = '(1.0-c)^2 + kdamage'
    constant_names = 'kdamage'
    constant_expressions = '1E-6'
    derivative_order = 2
  [../]
  [./indicator]
    type = DerivativeParsedMaterial
    f_name = indicator
    args = 'c'
    function = 'c'
    derivative_order = 2
  [../]
  [./local_fracture_energy]
    type = DerivativeParsedMaterial
    f_name = local_fracture_energy
    args = 'c'
    material_property_names = 'gc_prop l'
    function = 'c^2*gc_prop/2/l'
    derivative_order = 2
  [../]
  [./fracture_driving_energy]
    type = DerivativeSumMaterial
    f_name = F
    args = c
    sum_materials = 'elastic_energy local_fracture_energy'
    derivative_order = 2
  [../]
  [./crackfrictionheatenergy]
    type = ComputeCrackFrictionHeatEnergyDienes
    friction_coefficient = 0.0
    dcdx = dcdx
    dcdy = dcdy
    c = c
    l = 8E-3
  [../]
  [./plasticheatenergy]
    type = ComputePlasticHeatEnergy
  [../]
  [./shock_heat_source]
    type = ComputeThermalExpansionHeatSourceSmallStrainMieGruneisen
    Gamma = 0.7
    reference_temperature =300.0
    c = c
    temperature = temperature
    specific_heat = specific_heat
    density = density
    C0 = 0.1
    C1 = 1.0
    Le = 2E-3
    sound_speed = 2.77
    beta_v = 0.4
  [../]
  [./jetting_heat_source]
    type = ComputeMolecularJettingHeatSource
		c = c
    alpha=0
		up=0.5
		slope_UsUp=2.29
		sound_speed=3.33
		gap=0.08
		zeta=0.08
		dirac_switch=dirac_switch
		tau_induction=0.001
  [../]
  [./w_time_dirac_switch]
    type = ParsedMaterial
    args = 'pressure'
    f_name = dirac_switch_pressure
    function =  'if (pressure < -0.2 , 1.0 , 0.0)'
  [../]
[]
################################################################################

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]
################################################################################

[Executioner]
  type = Transient
  start_time = 0.0
  dt = 0.001E-3
  dtmin = 0.0001E-3
  end_time = 55E-3
  nl_rel_tol = 1E-6
  nl_abs_tol = 1E-6
  solve_type = PJFNK
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  petsc_options_value = '201 hypre boomeramg 20'
  line_search = 'none'
[]
################################################################################

[Outputs]
  interval = 100
  exodus = true
  checkpoint = true
[]
