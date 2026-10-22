[Problem]
  kernel_coverage_check = false
[]

# The same geometry and facet-mode radiation discretization as
# two_cubes_facet.i, except that the patches are indexed internally by the
# radiation objects (create_patch_sidesets defaults to false) and no patch
# sideset is added to the mesh.  The postprocessors therefore query the whole
# inner faces, and their area-weighted averages must reproduce the per-facet
# result of two_cubes_facet.i, where the four facets of each inner face see
# the same flux.
[Mesh]
  [box1]
    type = GeneratedMeshGenerator
    dim = 3
    xmax = '0.1'
    nx = 2
    ymax = '0.1'
    ny = 2
    zmax = '0.1'
    nz = 2
    subdomain_name = 'box1'
    subdomain_ids = 1
    boundary_id_offset = 0
    boundary_name_prefix = 'box1:'
  []
  [box2]
    type = GeneratedMeshGenerator
    dim = 3
    xmax = '0.1'
    nx = 2
    ymax = '0.1'
    ny = 2
    zmax = '0.1'
    nz = 2
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
    n_patches = '1 1 1 0 1 1 1 1 0 1 1 1'
    # n_patches = 0 is facet mode: one patch per face; the partitioner is not used
    temperature = temperature
    adiabatic_boundary = 'box1:_front box1:_back box1:_left box1:_top box1:_bottom
                          box2:_front box2:_back box2:_right box2:_top box2:_bottom'
    fixed_temperature_boundary = 'box1:_right box2:_left'
    fixed_boundary_temperatures = '400 300'
    polar_quad_order = 16
    azimuthal_quad_order = 8
    view_factor_calculator = vacuum_ray_tracing
    normalize_view_factor = false # radiation escaping to the vacuum cannot be redistributed
    environment = vacuum
  []
[]

[Postprocessors]
  # area-weighted averages over the facets of each inner face
  [flux_a]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box1:_right'
  []
  [flux_b]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box2:_left'
  []
  [temp_a]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = TEMPERATURE
    boundary = 'box1:_right'
  []
  [temp_b]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = TEMPERATURE
    boundary = 'box2:_left'
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Outputs]
  csv = true
[]
