# The 'right' sideset is copied into 'right_0', so the two sidesets share all of their
# faces. The radiation objects must reject such a boundary list instead of crediting the
# shared faces to one of the boundaries arbitrarily.
[Problem]
  kernel_coverage_check = false
[]

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 2
    nz = 2
  []
  [duplicate]
    type = PatchSidesetGenerator
    input = gmg
    boundary = 'right'
    n_patches = 1
  []
[]

[Variables]
  [temperature]
  []
[]

[GrayDiffuseRadiation]
  [cavity]
    boundary = 'left right right_0'
    emissivity = '1 1 1'
    n_patches = '1 1 1'
    temperature = temperature
    view_factor_calculator = ray_tracing
  []
[]

[Executioner]
  type = Steady
[]
