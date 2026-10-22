[Problem]
  kernel_coverage_check = false
[]

# Two 10 cm cubes separated by a 10 cm gap along x: cube A occupies
# [0, 0.1]^3 and cube B occupies [0.2, 0.3] x [0, 0.1]^2.  The two facing
# (inner) faces are held at 400 K and 300 K and are each subdivided into
# four patches (a 2x2 grid) by the PatchSidesetGenerator, exercising the
# per-patch radiosity model of the vacuum ray tracing view factor
# calculation.  Every other face is adiabatic and borders empty void.  With
# no conduction the temperature field is not updated by the solver; the
# radiative exchange is solved algebraically by the surface radiation
# userobject at each step, exactly like the other vacuum_radiator tests.
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
    initial_condition = 0
  []
[]

[GrayDiffuseRadiation]
  [cube_gap]
    boundary = 'box1:_front box1:_back box1:_left box1:_right box1:_top box1:_bottom
                box2:_front box2:_back box2:_left box2:_right box2:_top box2:_bottom'
    emissivity = '1 1 1 1 1 1 1 1 1 1 1 1'
    n_patches = '1 1 1 4 1 1 1 1 4 1 1 1'
    partitioners = 'metis metis metis grid metis metis metis metis grid metis metis metis'
    temperature = temperature
    adiabatic_boundary = 'box1:_front box1:_back box1:_left box1:_top box1:_bottom
                          box2:_front box2:_back box2:_right box2:_top box2:_bottom'
    fixed_temperature_boundary = 'box1:_right box2:_left'
    fixed_boundary_temperatures = '400 300'
    polar_quad_order = 16
    azimuthal_quad_order = 8
    view_factor_calculator = vacuum_ray_tracing
    create_patch_sidesets = true # materialize the patches as sidesets to query them by name
    normalize_view_factor = false # radiation escaping to the vacuum cannot be redistributed
    environment = vacuum
  []
[]

[Postprocessors]
  # per-patch heat flux densities and temperatures of the patched inner faces
  [flux_a_0]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box1:_right_0'
  []
  [flux_a_1]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box1:_right_1'
  []
  [flux_a_2]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box1:_right_2'
  []
  [flux_a_3]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box1:_right_3'
  []
  [flux_b_0]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box2:_left_0'
  []
  [flux_b_1]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box2:_left_1'
  []
  [flux_b_2]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box2:_left_2'
  []
  [flux_b_3]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box2:_left_3'
  []
  [temp_a_0]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = TEMPERATURE
    boundary = 'box1:_right_0'
  []
  [temp_a_1]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = TEMPERATURE
    boundary = 'box1:_right_1'
  []
  [temp_a_2]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = TEMPERATURE
    boundary = 'box1:_right_2'
  []
  [temp_a_3]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = TEMPERATURE
    boundary = 'box1:_right_3'
  []
  [temp_b_0]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = TEMPERATURE
    boundary = 'box2:_left_0'
  []
  [temp_b_1]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = TEMPERATURE
    boundary = 'box2:_left_1'
  []
  [temp_b_2]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = TEMPERATURE
    boundary = 'box2:_left_2'
  []
  [temp_b_3]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = TEMPERATURE
    boundary = 'box2:_left_3'
  []

  # view factors between every pair of patches of the two inner faces; the
  # patches are equal in area, so reciprocity requires F(a_i, b_j) = F(b_j, a_i)
  [vf_a0b0]
    type = ViewFactorPP
    from_boundary = 'box1:_right_0'
    to_boundary = 'box2:_left_0'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_a0b1]
    type = ViewFactorPP
    from_boundary = 'box1:_right_0'
    to_boundary = 'box2:_left_1'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_a0b2]
    type = ViewFactorPP
    from_boundary = 'box1:_right_0'
    to_boundary = 'box2:_left_2'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_a0b3]
    type = ViewFactorPP
    from_boundary = 'box1:_right_0'
    to_boundary = 'box2:_left_3'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_a1b0]
    type = ViewFactorPP
    from_boundary = 'box1:_right_1'
    to_boundary = 'box2:_left_0'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_a1b1]
    type = ViewFactorPP
    from_boundary = 'box1:_right_1'
    to_boundary = 'box2:_left_1'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_a1b2]
    type = ViewFactorPP
    from_boundary = 'box1:_right_1'
    to_boundary = 'box2:_left_2'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_a1b3]
    type = ViewFactorPP
    from_boundary = 'box1:_right_1'
    to_boundary = 'box2:_left_3'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_a2b0]
    type = ViewFactorPP
    from_boundary = 'box1:_right_2'
    to_boundary = 'box2:_left_0'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_a2b1]
    type = ViewFactorPP
    from_boundary = 'box1:_right_2'
    to_boundary = 'box2:_left_1'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_a2b2]
    type = ViewFactorPP
    from_boundary = 'box1:_right_2'
    to_boundary = 'box2:_left_2'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_a2b3]
    type = ViewFactorPP
    from_boundary = 'box1:_right_2'
    to_boundary = 'box2:_left_3'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_a3b0]
    type = ViewFactorPP
    from_boundary = 'box1:_right_3'
    to_boundary = 'box2:_left_0'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_a3b1]
    type = ViewFactorPP
    from_boundary = 'box1:_right_3'
    to_boundary = 'box2:_left_1'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_a3b2]
    type = ViewFactorPP
    from_boundary = 'box1:_right_3'
    to_boundary = 'box2:_left_2'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_a3b3]
    type = ViewFactorPP
    from_boundary = 'box1:_right_3'
    to_boundary = 'box2:_left_3'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_b0a0]
    type = ViewFactorPP
    from_boundary = 'box2:_left_0'
    to_boundary = 'box1:_right_0'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_b0a1]
    type = ViewFactorPP
    from_boundary = 'box2:_left_0'
    to_boundary = 'box1:_right_1'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_b0a2]
    type = ViewFactorPP
    from_boundary = 'box2:_left_0'
    to_boundary = 'box1:_right_2'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_b0a3]
    type = ViewFactorPP
    from_boundary = 'box2:_left_0'
    to_boundary = 'box1:_right_3'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_b1a0]
    type = ViewFactorPP
    from_boundary = 'box2:_left_1'
    to_boundary = 'box1:_right_0'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_b1a1]
    type = ViewFactorPP
    from_boundary = 'box2:_left_1'
    to_boundary = 'box1:_right_1'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_b1a2]
    type = ViewFactorPP
    from_boundary = 'box2:_left_1'
    to_boundary = 'box1:_right_2'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_b1a3]
    type = ViewFactorPP
    from_boundary = 'box2:_left_1'
    to_boundary = 'box1:_right_3'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_b2a0]
    type = ViewFactorPP
    from_boundary = 'box2:_left_2'
    to_boundary = 'box1:_right_0'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_b2a1]
    type = ViewFactorPP
    from_boundary = 'box2:_left_2'
    to_boundary = 'box1:_right_1'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_b2a2]
    type = ViewFactorPP
    from_boundary = 'box2:_left_2'
    to_boundary = 'box1:_right_2'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_b2a3]
    type = ViewFactorPP
    from_boundary = 'box2:_left_2'
    to_boundary = 'box1:_right_3'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_b3a0]
    type = ViewFactorPP
    from_boundary = 'box2:_left_3'
    to_boundary = 'box1:_right_0'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_b3a1]
    type = ViewFactorPP
    from_boundary = 'box2:_left_3'
    to_boundary = 'box1:_right_1'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_b3a2]
    type = ViewFactorPP
    from_boundary = 'box2:_left_3'
    to_boundary = 'box1:_right_2'
    view_factor_object_name = view_factor_uo_cube_gap
  []
  [vf_b3a3]
    type = ViewFactorPP
    from_boundary = 'box2:_left_3'
    to_boundary = 'box1:_right_3'
    view_factor_object_name = view_factor_uo_cube_gap
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Outputs]
  csv = true
[]
