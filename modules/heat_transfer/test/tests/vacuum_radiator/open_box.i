[Problem]
  kernel_coverage_check = false
[]

# An open (hollow) box: a solid unit cube made of 10x10x10 cells whose interior
# x,z in [0.1, 0.9] and y in [0.1, 1.0] is deleted, leaving a cup with a bottom
# wall, four side walls and an opening at the top.  The walls facing the hollow
# (inner_bottom + the four inner walls) exchange radiation with each other
# through the hollow, while the outer faces and the rim around the opening only
# radiate into the void (escape fraction one).
[Mesh]
  [box]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 10
    nz = 10
  []

  [inner_region]
    type = SubdomainBoundingBoxGenerator
    input = box
    block_id = 1
    bottom_left = '0.1 0.1 0.1'
    top_right = '0.9 1.0 0.9'
  []

  # no new_boundary is requested: the cut surface between the kept cells and the
  # deleted interior is created directly below with SideSetsFromPointsGenerator
  [delete_inner]
    type = BlockDeletionGenerator
    input = inner_region
    block = 1
  []

  # split the cut surface (the walls of the hollow) into one sideset per wall
  [split_inner]
    type = SideSetsFromPointsGenerator
    input = delete_inner
    new_boundary = 'inner_bottom inner_front inner_back inner_left inner_right'
    points = '0.55 0.1 0.55
              0.55 0.55 0.1
              0.55 0.55 0.9
              0.1 0.55 0.55
              0.9 0.55 0.55'
  []
[]

[Variables]
  [temperature]
    initial_condition = 300
  []
[]

[GrayDiffuseRadiation]
  [open_box]
    boundary = 'bottom top left right front back inner_bottom inner_front inner_back inner_left inner_right'
    emissivity = '1 1 1 1 1 1 1 1 1 1 1'
    n_patches = '1 1 1 1 1 1 1 1 1 1 1'
    temperature = temperature
    adiabatic_boundary = 'bottom top left right front back inner_front inner_back inner_left inner_right'
    fixed_temperature_boundary = 'inner_bottom'
    fixed_boundary_temperatures = '1000'
    polar_quad_order = 16
    azimuthal_quad_order = 8
    view_factor_calculator = vacuum_ray_tracing
    environment = vacuum
  []
[]

[Postprocessors]
  [heat_flux_density_inner_bottom]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_open_box
    return_type = HEAT_FLUX_DENSITY
    boundary = inner_bottom
  []

  [temperature_inner_front]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_open_box
    return_type = TEMPERATURE
    boundary = inner_front
  []

  [temperature_inner_back]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_open_box
    return_type = TEMPERATURE
    boundary = inner_back
  []

  [temperature_inner_left]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_open_box
    return_type = TEMPERATURE
    boundary = inner_left
  []

  [temperature_inner_right]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_open_box
    return_type = TEMPERATURE
    boundary = inner_right
  []

  [temperature_top]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_open_box
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