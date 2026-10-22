# Diffusion with a source on the unit square mirrored by (x, y) -> (-x, y) with
# TransformGenerator: the mesh stays in the xy-plane and the node order of every element becomes
# clockwise. An element map that takes the signed determinant of the element Jacobian as the
# integration weight would give every element a negative weight and the mesh an area of -1; the
# host assembly gives every element a positive weight and the mesh an area of 1.

[Mesh]
  [square]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 8
    ny = 8
  []
  [mirror]
    type = TransformGenerator
    input = square
    transform = SCALE
    vector_value = '-1 1 1'
  []
[]

[Variables]
  [u]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  [one]
    initial_condition = 1
  []
[]

[ICs]
  [u]
    type = FunctionIC
    variable = u
    function = '1 + x + 2*y + 3*z'
  []
[]

[Kernels]
  [dt]
    type = KokkosTimeDerivative
    variable = u
  []
  [diff]
    type = KokkosDiffusion
    variable = u
  []
  [force]
    type = KokkosBodyForce
    variable = u
    value = 1
  []
[]

[BCs]
  [left]
    type = KokkosDirichletBC
    variable = u
    boundary = left
    value = 0
  []
  [right]
    type = KokkosNeumannBC
    variable = u
    boundary = right
    value = 1
  []
[]

[Postprocessors]
  [area]
    type = KokkosElementIntegralVariablePostprocessor
    variable = one
  []
  [mass]
    type = KokkosElementIntegralVariablePostprocessor
    variable = u
  []
  [l2]
    type = ElementL2Norm
    variable = u
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  dt = 0.1
  num_steps = 2
[]

[Outputs]
  exodus = true
[]
