# No kernel: one column of four elements, so every node lies on the left or right boundary and
# carries a Kokkos Dirichlet boundary condition, and two processes can each own elements. The
# residual and Jacobian are computed together, and the combined Kokkos pass therefore runs only the
# nodal boundary conditions, which compute their Jacobian there and their residual in their own
# pass, so it activates no residual tag. The boundary conditions are not preset, so the solve takes
# one Newton step with the combined pass's Jacobian and the nodal pass's residual.

[Mesh]
  [square]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 1
    ny = 4
  []
[]

[Variables]
  [u]
    order = FIRST
    family = LAGRANGE
  []
[]

[BCs]
  [left]
    type = KokkosADDirichletBC
    variable = u
    boundary = left
    value = 0
    preset = false
  []
  [right]
    type = KokkosADDirichletBC
    variable = u
    boundary = right
    value = 1
    preset = false
  []
[]

[Problem]
  kernel_coverage_check = false
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  residual_and_jacobian_together = true
[]

[Outputs]
  exodus = true
[]
