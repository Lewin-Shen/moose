# Diffusion with a source on the unit square rotated out of the xy-plane. The Euler angles
# (20, 50, 10) degrees of TransformGenerator give a rotation with no zero entry, and the normal of
# the plane becomes (0.1330, -0.7544, 0.6428), at 50 degrees to the z-axis. An element map built
# from the x and y components of the element tangents alone would scale every quadrature weight
# by cos(50 degrees) and give wrong gradients.

[Mesh]
  [square]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 8
    ny = 8
  []
  [rotate]
    type = TransformGenerator
    input = square
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
