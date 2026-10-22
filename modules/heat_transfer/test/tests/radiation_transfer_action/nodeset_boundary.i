# 'inlet' is a node set, not a side set; it has no faces to exchange radiation.
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
  [inlet_nodes]
    type = ExtraNodesetGenerator
    input = gmg
    new_boundary = 'inlet'
    coord = '0 0 0'
  []
[]

[Variables]
  [temperature]
  []
[]

[GrayDiffuseRadiation]
  [cavity]
    boundary = 'left right inlet'
    emissivity = '1 1 1'
    n_patches = '1 1 1'
    temperature = temperature
    view_factor_calculator = ray_tracing
  []
[]

[Executioner]
  type = Steady
[]
