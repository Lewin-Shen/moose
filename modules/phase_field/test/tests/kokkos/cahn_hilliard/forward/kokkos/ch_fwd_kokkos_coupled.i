# Cross-derivative Jacobian gate for the tier 2 kernels (implementation decisions I-3, I-4, I-10): the
# forward split on (c, w) plus a third variable eta with its own equation, coupled through the test-only
# two-variable energy f = W (c - ca)^2 (c - cb)^2 + a c eta^2 (KokkosCHFreeEnergyCoupledTest) and the
# degenerate mobility M = M0 c (1 - c) (KokkosCHMobility). Every Jacobian branch of the two kernels is
# exercised with a nonzero value:
#   c_MatDiffusion       : v = w, D = M, args = eta                -> dM/deta is not declared (unset handle)
#   w_CoupledMaterialDerivative : v = c, coupled_variables = eta    -> d^2f/dcdeta block (w, eta)
#   eta_MatDiffusion     : no v, D = M, args = c                     -> Laplacian on the diagonal, dM/dc block
#   eta_CrossDiffusion   : v = c, D = M                              -> dM/dc phi grad(c) + M grad(phi) block
#   eta_Bulk             : v = eta (the kernel variable), coupled_variables = c
#                                                                    -> diagonal d^2f/deta^2, block d^2f/dcdeta
# Unit-scaled coefficients as in ch_fwd_kokkos_degenerate.i. The eta equation is
# eta_t = div(M grad(eta)) + div(M grad(c)) - df/deta, stable (df/deta = 2 a c eta). CPU control:
# ../cpu_control/ch_fwd_cpu_coupled.i (DerivativeParsedMaterial twins), prefixes and material blocks only.
#
# Field scales for the exodiff floor: c in [0.3, 0.7], eta in [0.1, 0.5], w of order 0.1 to 1.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 16
  ny = 16
[]

[Variables]
  [c]
  []
  [w]
  []
  [eta]
  []
[]

[Functions]
  [c_ic]
    type = ParsedFunction
    expression = '0.5 + 0.2*cos(2*pi*x)*cos(2*pi*y)'
  []
  [eta_ic]
    type = ParsedFunction
    expression = '0.3 + 0.2*sin(pi*x)*sin(pi*y)'
  []
[]

[ICs]
  [c_ic]
    type = FunctionIC
    variable = c
    function = c_ic
  []
  [eta_ic]
    type = FunctionIC
    variable = eta
    function = eta_ic
  []
[]

[Kernels]
  [c_TimeDerivative]
    type = KokkosTimeDerivative
    variable = c
  []
  [c_MatDiffusion]
    type = KokkosMatDiffusion
    variable = c
    v = w
    diffusivity = M
    args = 'eta'
  []
  [w_MatDiffusion]
    type = KokkosMatDiffusion
    variable = w
    v = c
    diffusivity = kappa_c
  []
  [w_CoupledMaterialDerivative]
    type = KokkosCoupledMaterialDerivative
    variable = w
    v = c
    f_name = f
    coupled_variables = 'eta'
  []
  [w_CoefReaction]
    type = KokkosReaction
    variable = w
    rate = -1
  []
  [eta_TimeDerivative]
    type = KokkosTimeDerivative
    variable = eta
  []
  [eta_MatDiffusion]
    type = KokkosMatDiffusion
    variable = eta
    diffusivity = M
    args = 'c'
  []
  [eta_CrossDiffusion]
    type = KokkosMatDiffusion
    variable = eta
    v = c
    diffusivity = M
  []
  [eta_Bulk]
    type = KokkosCoupledMaterialDerivative
    variable = eta
    v = eta
    f_name = f
    coupled_variables = 'c'
  []
[]

[Materials]
  [kappa]
    type = KokkosGenericConstantMaterial
    prop_names = 'kappa_c'
    prop_values = '0.01'
  []
  [mobility]
    type = KokkosCHMobility
    property_name = M
    c = c
    form = DEGENERATE
    M0 = 1
  []
  [free_energy]
    type = KokkosCHFreeEnergyCoupledTest
    property_name = f
    c = c
    eta = eta
    W = 1
    a = 1
  []
[]

[Postprocessors]
  [c_total]
    type = KokkosElementIntegralVariablePostprocessor
    variable = c
    execute_on = 'initial timestep_end'
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  automatic_scaling = false
  # direct solve and full Newton steps: the backtracking line search creeps once the residual norm is
  # near its round-off floor (Campaign 02's E10 observation); the solver is not the point of the gates
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type -snes_linesearch_type'
  petsc_options_value = 'lu mumps basic'
  nl_rel_tol = 1e-11
  nl_abs_tol = 1e-12
  nl_max_its = 20
  l_max_its = 30
  dt = 0.05
  num_steps = 3
[]

[Outputs]
  exodus = true
  csv = true
[]
