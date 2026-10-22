# Diffusion with a source on the unit segment rotated off the x-axis by the Euler angles
# (20, 50, 10) degrees of TransformGenerator: the direction of the segment becomes
# (0.8872, 0.3797, 0.2620), at 27.5 degrees to the x-axis. An element map built from the
# x component of the element tangent alone would scale every quadrature weight by 0.8872 and
# give wrong gradients.

[Mesh]
  [segment]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 8
  []
  [rotate]
    type = TransformGenerator
    input = segment
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
