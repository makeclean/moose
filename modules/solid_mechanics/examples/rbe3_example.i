# Example input file for RBE3 constraint usage

[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 0
  xmax = 10
  ymin = 0
  ymax = 10
  nx = 10
  ny = 10
[]

[Variables]
  [displacement]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  [elasticity]
    type = Elasticity
    variable = displacement
  []
[]

[Constraints]
  # Example 1: Explicit node list with explicit weights
  [rbe3_explicit]
    type = RBE3Constraint
    primary_nodes = [1, 2, 3]
    secondary_nodes = [4]
    weights = [0.5, 0.3, 0.2]
  []

  # Example 2: Sideset-based with automatic equal weights  
  [rbe3_sideset_equal]
    type = RBE3Constraint
    primary_sideset = "primary_boundary"
    secondary_nodes = [5]
    weight_method = "equal"
  []

  # Example 3: Sideset-based with distance weights
  [rbe3_sideset_distance]
    type = RBE3Constraint
    primary_sideset = "secondary_boundary"
    secondary_nodes = [6]
    weight_method = "distance"
  []
[]

[BoundaryConditions]
  [fixed]
    type = DirichletBC
    variable = displacement
    boundary = "left"
    value = 0
  []
  
  [load]
    type = DirichletBC
    variable = displacement
    boundary = "right"
    value = 1
  []
[]

[Outputs]
  exodus = true
[]