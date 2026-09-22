[Problem]
  kernel_coverage_check = false
[]

# Same convex element as in convex_escape.i, but the open geometry now exchanges
# radiation with a black body environment at 500 K instead of a vacuum.
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
[]

[Functions]
  [Tenv]
    type = ConstantFunction
    value = 500
  []
[]

[Variables]
  [temperature]
    initial_condition = 300
  []
[]

[GrayDiffuseRadiation]
  [open_world]
    boundary = 'bottom top left right front back'
    emissivity = '1 1 1 1 1 1'
    n_patches = '1 1 1 1 1 1'
    temperature = temperature
    adiabatic_boundary = 'top left right front back'
    fixed_temperature_boundary = 'bottom'
    fixed_boundary_temperatures = '1000'
    polar_quad_order = 16
    azimuthal_quad_order = 8
    view_factor_calculator = vacuum_ray_tracing
    environment = black_body
    environment_temperature = Tenv
  []
[]

[Postprocessors]
  [heat_flux_density_bottom]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_open_world
    return_type = HEAT_FLUX_DENSITY
    boundary = bottom
  []

  [temperature_top]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_open_world
    return_type = TEMPERATURE
    boundary = top
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Outputs]
  csv = true
[]