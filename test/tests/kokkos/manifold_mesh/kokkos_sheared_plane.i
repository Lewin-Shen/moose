# As kokkos_tilted_plane.i, with the unit square sheared to a parallelogram by
# (x, y) -> (x + 0.3 y, y) before the rotation: the two tangents of each element are not
# orthogonal, so the surface metric of each element has an off-diagonal term. The area stays 1.

[Mesh]
  [square]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 8
    ny = 8
  []
  [shear]
    type = ParsedNodeTransformGenerator
    input = square
    x_function = 'x + 0.3*y'
  []
  [rotate]
    type = TransformGenerator
    input = shear
    transform = ROTATE
    vector_value = '20 50 10'
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
