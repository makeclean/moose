[Problem]
  kernel_coverage_check = false
[]

# A 10 cm box occupying [0, 0.1]^3 and a solid sphere of radius 5 cm
# centered at (0.2, 0.05, 0.05), leaving a 5 cm gap between the box right
# face and the sphere.  The box left face is held at 4000 K, the other
# four box faces are adiabatic, and the box right face and the sphere
# surface are free, each subdivided into one patch per facet (n_patches =
# 0) by the PatchSidesetGenerator.  In steady state heat flows from the box
# left face by conduction to the box right face and then by radiation
# across the gap to the sphere, which warms from its initial temperature
# and radiates to the surrounding vacuum.
[Mesh]
  [sphere]
    type = SphereMeshGenerator
    nr = 3
    radius = 0.05
    elem_type = TET4
  []
  [rename]
    type = RenameBoundaryGenerator
    input = sphere
    old_boundary = '0'
    new_boundary = 'sphere:_surface'
  []
  [box]
    type = GeneratedMeshGenerator
    dim = 3
    xmax = '0.1'
    nx = 2
    ymax = '0.1'
    ny = 2
    zmax = '0.1'
    nz = 2
    subdomain_name = 'box'
    subdomain_ids = 1
    boundary_id_offset = 1
    boundary_name_prefix = 'box:'
  []
  [transform]
    type = TransformGenerator
    input = rename
    transform = TRANSLATE
    vector_value = '0.2 0.05 0.05'
  []
  [cmbn]
    type = CombinerGenerator
    inputs = 'box transform'
  []
[]

[Variables]
  [temperature]
    initial_condition = 200
  []
[]

[BCs]
  [fixed_temp]
    type = DirichletBC
    boundary = 'box:_left'
    value = 4000
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
    boundary = 'box:_front box:_back box:_left box:_right box:_top box:_bottom
                sphere:_surface'
    emissivity = '1 1 1 1 1 1 1'
    n_patches = '1 1 1 0 1 1 0'
    # n_patches = 0 is facet mode: one patch per face; the partitioner is not used
    partitioners = 'metis metis metis metis metis metis metis'
    temperature = temperature
    adiabatic_boundary = 'box:_front box:_back box:_top box:_bottom'
    polar_quad_order = 16
    azimuthal_quad_order = 8
    view_factor_calculator = vacuum_ray_tracing
    create_patch_sidesets = true # materialize the patches as sidesets to query them by name
    normalize_view_factor = false # radiation escaping to the vacuum cannot be redistributed
    environment = vacuum
  []
[]

[Postprocessors]
  # temperature of the sphere: its center and the surface point closest to
  # the box
  [sphere_center_temp]
    type = PointValue
    variable = temperature
    point = '0.2 0.05 0.05'
  []
  [sphere_face_temp]
    type = PointValue
    variable = temperature
    point = '0.15 0.05 0.05'
  []

  # per-facet net radiative flux on the box right face, the only box face
  # that can exchange net heat
  [box_right_flux_0]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box:_right_0'
  []
  [box_right_flux_1]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box:_right_1'
  []
  [box_right_flux_2]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box:_right_2'
  []
  [box_right_flux_3]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box:_right_3'
  []

  # net radiative flux on the adiabatic box faces, which must be zero
  [adiabatic_flux_front]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box:_front'
  []
  [adiabatic_flux_back]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box:_back'
  []
  [adiabatic_flux_top]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box:_top'
  []
  [adiabatic_flux_bottom]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = view_factor_surface_radiation_cube_gap
    return_type = HEAT_FLUX_DENSITY
    boundary = 'box:_bottom'
  []
[]

[Executioner]
  type = Steady
  petsc_options_iname = '-pc_type'
  petsc_options_value = ' lu'
  nl_rel_tol = 1e-10
  nl_max_its = 20
  l_tol = 1e-8
  l_max_its = 100
[]

[Outputs]
  csv = true
[]
