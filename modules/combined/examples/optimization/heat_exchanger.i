# Heat Exchanger Topology Optimization with Brinkman-Stokes Flow
# 
# This example demonstrates topology optimization of a heat exchanger using
# the SIMP (Solid Isotropic Material with Penalization) method.
#
# The design minimizes thermal compliance (maximizes heat transfer) while
# constraining the volume fraction of solid material.
#
# Physics:
# - Brinkman-Stokes flow (PINSFV) through porous medium
# - Heat conduction + advection
# - Design density controls porosity, permeability, and thermal conductivity

vol_frac = 0.4
power = 3.0

# Material properties
k_solid = 10.0
k_fluid = 0.5
perm_solid = 1e-12
perm_fluid = 1e-7
mu = 0.001
rho = 1.0

[GlobalParams]
  advected_interp_method = 'average'
  velocity_interp_method = 'rc'
  rhie_chow_user_object = 'rc'
  porosity = porosity
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 40
    ny = 20
    xmin = 0
    xmax = 2
    ymin = 0
    ymax = 1
  []
[]

[Problem]
  nl_sys_names = 'u_system v_system pressure_system'
  previous_nl_solution_required = true
[]

[UserObjects]
  [rc]
    type = PINSFVRhieChowInterpolatorSegregated
    u = superficial_vel_x
    v = superficial_vel_y
    pressure = pressure
    porosity = porosity
  []
  [radial_filter]
    type = RadialAverage
    radius = 0.12
    weights = linear
    prop_name = sensitivity_filtered
    execute_on = 'TIMESTEP_END'
    variable = Dc
  []
  [density_update]
    type = DensityUpdate
    density_sensitivity = sensitivity_filtered
    design_density = mat_den
    volume_fraction = ${vol_frac}
    execute_on = 'TIMESTEP_BEGIN'
  []
[]

[Variables]
  [superficial_vel_x]
    type = PINSFVSuperficialVelocityVariable
    initial_condition = 1.0
    solver_sys = u_system
  []
  [superficial_vel_y]
    type = PINSFVSuperficialVelocityVariable
    initial_condition = 0.0
    solver_sys = v_system
  []
  [pressure]
    type = INSFVPressureVariable
    solver_sys = pressure_system
  []
  [temperature]
    family = LAGRANGE
    order = FIRST
  []
[]

[AuxVariables]
  [mat_den]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = ${vol_frac}
  []
  [porosity]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 0.5
  []
  [Dc]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = -1.0
  []
  [k_phys]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 1.0
  []
  [permeability]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 1e-10
  []
  [speed]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 1.0
  []
[]

[AuxKernels]
  [porosity]
    type = ParsedAux
    variable = porosity
    expression = 'mat_den'
    coupled_variables = 'mat_den'
    execute_on = 'INITIAL LINEAR'
  []
  [k_phys]
    type = ParsedAux
    variable = k_phys
    expression = '${k_fluid} + (mat_den ^ ${power}) * (${k_solid} - ${k_fluid})'
    coupled_variables = 'mat_den'
    execute_on = 'INITIAL LINEAR'
  []
  [permeability]
    type = ParsedAux
    variable = permeability
    expression = '${perm_solid} + (mat_den ^ ${power}) * (${perm_fluid} - ${perm_solid})'
    coupled_variables = 'mat_den'
    execute_on = 'INITIAL LINEAR'
  []
  [speed]
    type = ParsedAux
    variable = speed
    expression = 'sqrt(superficial_vel_x^2 + superficial_vel_y^2)'
    coupled_variables = 'superficial_vel_x superficial_vel_y'
    execute_on = 'LINEAR'
  []
  [Dc]
    type = MaterialRealAux
    variable = Dc
    material_property = thermal_flow_sensitivity
    execute_on = 'LINEAR TIMESTEP_END'
  []
[]

[FVKernels]
  # X-momentum
  [u_advection]
    type = PINSFVMomentumAdvection
    variable = superficial_vel_x
    rho = ${rho}
    porosity = porosity
    momentum_component = x
  []
  [u_diffusion]
    type = PINSFVMomentumDiffusion
    variable = superficial_vel_x
    mu = ${mu}
    porosity = porosity
    momentum_component = x
  []
  [u_pressure]
    type = PINSFVMomentumPressure
    variable = superficial_vel_x
    momentum_component = x
    pressure = pressure
    porosity = porosity
  []

  # Y-momentum
  [v_advection]
    type = PINSFVMomentumAdvection
    variable = superficial_vel_y
    rho = ${rho}
    porosity = porosity
    momentum_component = y
  []
  [v_diffusion]
    type = PINSFVMomentumDiffusion
    variable = superficial_vel_y
    mu = ${mu}
    porosity = porosity
    momentum_component = y
  []
  [v_pressure]
    type = PINSFVMomentumPressure
    variable = superficial_vel_y
    momentum_component = y
    pressure = pressure
    porosity = porosity
  []

  # Pressure equation
  [p_diffusion]
    type = FVAnisotropicDiffusion
    variable = pressure
    coeff = "Ainv"
    coeff_interp_method = 'average'
  []
  [p_source]
    type = FVDivergence
    variable = pressure
    vector_field = "HbyA"
    force_boundary_execution = true
  []
[]

[Kernels]
  [heat_conduction]
    type = HeatConduction
    variable = temperature
    thermal_conductivity = k_phys
  []
  [heat_advection]
    type = ConservativeAdvection
    variable = temperature
    velocity = 'superficial_vel_x superficial_vel_y'
  []
[]

[Materials]
  [k_phys_mat]
    type = DerivativeParsedMaterial
    expression = '${k_fluid} + (mat_den ^ ${power}) * (${k_solid} - ${k_fluid})'
    coupled_variables = 'mat_den'
    property_name = k_phys
    derivative_order = 1
  []
  [permeability_mat]
    type = DerivativeParsedMaterial
    expression = '${perm_solid} + (mat_den ^ ${power}) * (${perm_fluid} - ${perm_solid})'
    coupled_variables = 'mat_den'
    property_name = permeability
    derivative_order = 1
  []
  [thermal_sensitivity]
    type = ThermalFlowSensitivity
    temperature = temperature
    design_density = mat_den
    thermal_conductivity = k_phys
  []
  [thermal_compliance]
    type = ThermalCompliance
    temperature = temperature
    thermal_conductivity = k_phys
  []
[]

[FVBCs]
  [inlet_x]
    type = PINSFVDirichletBoundaryFlow
    variable = superficial_vel_x
    velocity = 1.0
    boundary = left
  []
  [inlet_y]
    type = PINSFVDirichletBoundaryFlow
    variable = superficial_vel_y
    velocity = 0.0
    boundary = left
  []
  [wall_no_slip]
    type = PINSFVFullySlippedWallBC
    variable = superficial_vel_x
    boundary = 'top bottom'
  []
  [wall_no_slip_y]
    type = PINSFVFullySlippedWallBC
    variable = superficial_vel_y
    boundary = 'top bottom'
  []
[]

[BCs]
  [cold_bc]
    type = DirichletBC
    variable = temperature
    boundary = left
    value = 0.0
  []
  [hot_bc]
    type = DirichletBC
    variable = temperature
    boundary = right
    value = 1.0
  []
  [insulation]
    type = NeumannBC
    variable = temperature
    boundary = 'top bottom'
    value = 0.0
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-6
[]

[Outputs]
  [exodus]
    type = Exodus
    execute_on = 'FINAL'
  []
  [csv]
    type = CSV
    execute_on = 'TIMESTEP_END'
  []
[]

[Postprocessors]
  [mesh_volume]
    type = VolumePostprocessor
  []
  [total_vol]
    type = ElementIntegralVariablePostprocessor
    variable = mat_den
  []
  [vol_frac]
    type = ParsedPostprocessor
    expression = 'total_vol / mesh_volume'
    pp_names = 'total_vol mesh_volume'
  []
  [thermal_compliance]
    type = ElementIntegralMaterialPropertyPostprocessor
    material_property = thermal_compliance
  []
  [avg_temp]
    type = AverageVariableValue
    variable = temperature
  []
  [max_velocity]
    type = ElementalL2Error
    variable = superficial_vel_x
    function = 1.0
  []
[]
