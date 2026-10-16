[Problem]
  kernel_coverage_check = false
[]

# Two 10 cm cubes separated by a 10 cm gap along x: box1 occupies
# [0, 0.1]^3 and box2 occupies [0.2, 0.3] x [0, 0.1]^2.  The bottom faces
# are held at 4000 K (box1) and 400 K (box2), the two facing (inner) faces
# exchange radiation across the gap, and the eight remaining outer faces
# are adiabatic.  In steady state heat flows from the hot cube to the cold
# cube by conduction through each cube and radiation across the gap.
[Mesh]
  [box1]
    type = GeneratedMeshGenerator
    dim = 3
    xmax = '0.1'
    nx = 10
    ymax = '0.1'
    ny = 10
    zmax = '0.1'
    nz = 10
    subdomain_name = 'box1'
    subdomain_ids = 1
    boundary_id_offset = 0
    boundary_name_prefix = 'box1:'
  []
  [box2]
    type = GeneratedMeshGenerator
    dim = 3
    xmax = '0.1'
    nx = 10
    ymax = '0.1'
    ny = 10
    zmax = '0.1'
    nz = 10
    subdomain_name = 'box2'
    subdomain_ids = 2
    boundary_id_offset = 6
    boundary_name_prefix = 'box2:'
  []
  [transform]
    type = TransformGenerator
    input = box2
    transform = TRANSLATE
    vector_value = '0.2 0.0 0.0'
  []
  [cmbn]
    type = CombinerGenerator
    inputs = 'box1 transform' 
  []
[]

[Variables]
  [temperature]
    family = LAGRANGE
    order = FIRST
    initial_condition = 0
  []
[]

[BCs]
  [fixed_temp2]
    type = DirichletBC
    boundary = 'box1:_bottom'
    value = 4000
    variable = temperature
  []
  [fixed_temp]
    type = DirichletBC
    boundary = 'box2:_bottom'
    value = 400
    variable = temperature
  []  
[]
  
[Kernels]
  [heat_conduction]
    type = HeatConduction
    variable = temperature
  []
[]

[Materials]
  [thermal_conductivity]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = 390.0
  []
[]
  
[GrayDiffuseRadiation]
  [cube_gap]
    boundary = 'box1:_front box1:_back box1:_left box1:_right box1:_top box1:_bottom
                box2:_front box2:_back box2:_left box2:_right box2:_top box2:_bottom'
    emissivity = '1 1 1 1 1 1 1 1 1 1 1 1'
    n_patches =  '1 1 1 1 1 1 1 1 1 1 1 1'
    temperature = temperature
    adiabatic_boundary = 'box1:_front box1:_back box1:_left box1:_top 
                          box2:_front box2:_back box2:_right box2:_top'
    #variable_temperature = 
    polar_quad_order = 16
    azimuthal_quad_order = 8
    view_factor_calculator = vacuum_ray_tracing
    normalize_view_factor = false # radiation escaping to the vacuum cannot be redistributed
    environment = vacuum
    # print_view_factor_info = true
  []
[]

[Postprocessors]
  [heat_flux_density_a_in]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box1:_right'
  []

  [heat_flux_density_b_in]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box2:_left'
  []

  # the outer faces border only empty void, so they must exchange no
  # radiation and equilibrate to 0 K
  [view_factor_b1b2]
    type = ViewFactorPP
    from_boundary = 'box1:_right'
    to_boundary = 'box2:_left'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [view_factor_b2b1]
    type = ViewFactorPP
    from_boundary = 'box2:_left'
    to_boundary = 'box1:_right'
    view_factor_object_name = view_factor_uo_cube_gap
  []
[]

[Executioner]
  type = Steady
  petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type '
  petsc_options_value = '   hypre      boomeramg                                  0.7                             HMIS                           ext+i '
  nl_rel_tol = 1e-10
  nl_max_its = 10
  l_tol = 1e-8
  l_max_its = 100
[]

[Outputs]
  exodus = true
  csv = true
[]
