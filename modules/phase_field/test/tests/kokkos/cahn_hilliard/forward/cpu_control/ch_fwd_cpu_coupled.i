# Explicit CPU control of ../kokkos/ch_fwd_kokkos_coupled.i: the same (c, w, eta) system with the
# two-variable energy f = W (c - ca)^2 (c - cb)^2 + a c eta^2 and the degenerate mobility both supplied by
# DerivativeParsedMaterial (which declares every first and second derivative, including d^2f/dcdeta).
# Differs from the Kokkos deck by the object prefixes and the two material blocks only.
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
    type = TimeDerivative
    variable = c
  []
  [c_MatDiffusion]
    type = MatDiffusion
    variable = c
    v = w
    diffusivity = M
    args = 'eta'
  []
  [w_MatDiffusion]
    type = MatDiffusion
    variable = w
    v = c
    diffusivity = kappa_c
  []
  [w_CoupledMaterialDerivative]
    type = CoupledMaterialDerivative
    variable = w
    v = c
    f_name = f
    coupled_variables = 'eta'
  []
  [w_CoefReaction]
    type = Reaction
    variable = w
    rate = -1
  []
  [eta_TimeDerivative]
    type = TimeDerivative
    variable = eta
  []
  [eta_MatDiffusion]
    type = MatDiffusion
    variable = eta
    diffusivity = M
    args = 'c'
  []
  [eta_CrossDiffusion]
    type = MatDiffusion
    variable = eta
    v = c
    diffusivity = M
  []
  [eta_Bulk]
    type = CoupledMaterialDerivative
    variable = eta
    v = eta
    f_name = f
    coupled_variables = 'c'
  []
[]

[Materials]
  [kappa]
    type = GenericConstantMaterial
    prop_names = 'kappa_c'
    prop_values = '0.01'
  []
  [mobility]
    type = DerivativeParsedMaterial
    property_name = M
    coupled_variables = c
    constant_names = 'M0'
    constant_expressions = '1'
    expression = 'M0*c*(1-c)'
  []
  [free_energy]
    type = DerivativeParsedMaterial
    property_name = f
    coupled_variables = 'c eta'
    constant_names = 'W ca cb a'
    constant_expressions = '1 0.27 0.83 1'
    expression = 'W*(c-ca)^2*(c-cb)^2 + a*c*eta^2'
  []
[]

[Postprocessors]
  [c_total]
    type = ElementIntegralVariablePostprocessor
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
