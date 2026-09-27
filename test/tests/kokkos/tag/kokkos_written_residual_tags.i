# Transient diffusion with Kokkos kernels. By default the kernels write only the time and non-time
# residual tags, so no Kokkos object of the residual pass writes the residual tag. The tests switch
# the diffusion kernel to the residual tag (Kernels/diff/vector_tags=residual) or switch on a Kokkos
# nodal boundary condition (BCs/active=left), which writes the residual tag in its own pass.

[Mesh]
  [square]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = 20
  []
[]

[Variables]
  [u]
    order = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  [ic]
    type = ParsedFunction
    expression = '0.5 + 0.4 * sin(2 * pi * x) * cos(2 * pi * y)'
  []
[]

[ICs]
  [u_ic]
    type = FunctionIC
    variable = u
    function = ic
  []
[]

[Kernels]
  [time]
    type = KokkosTimeDerivative
    variable = u
  []
  [diff]
    type = KokkosDiffusion
    variable = u
  []
[]

[BCs]
  active = ''
  [left]
    type = KokkosDirichletBC
    variable = u
    boundary = left
    value = 0
  []
[]

[Postprocessors]
  [u_l2]
    type = ElementL2Norm
    variable = u
  []
  [u_max]
    type = NodalExtremeValue
    variable = u
    value_type = max
  []
  [u_min]
    type = NodalExtremeValue
    variable = u
    value_type = min
  []
[]

[Executioner]
  type = Transient
  scheme = implicit-euler
  solve_type = PJFNK
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'jacobi'
  nl_rel_tol = 1e-10
  dt = 1e-3
  num_steps = 3
[]

[Outputs]
  csv = true
[]
