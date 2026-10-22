# Diffusion with a source on an open cylinder of radius 1 about the z-axis: the strip
# [0, 3 pi / 2] x [0, 1] is wrapped by (x, y) -> (cos x, sin x, y). The tangent plane of every
# element contains the z direction, so an element map built from the x and y components of the
# element tangents alone would be singular. The area of the mesh is 12 * 2 * sin(pi / 16).

[Mesh]
  [strip]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 12
    ny = 4
    xmax = 4.71238898038469
  []
  [wrap]
    type = ParsedNodeTransformGenerator
    input = strip
    x_function = 'cos(x)'
    y_function = 'sin(x)'
    z_function = 'y'
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
