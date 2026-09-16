# Explicit CPU control of the coupled-argument variant of the reverse-split Cahn-Hilliard gate deck (port 02, tier 1; decisions I-3
# and I-4): the free energy gains a second nonlinear argument eta through the term a c eta^2, so
# SplitCHParsed carries the cross derivative d^2F/dc deta in its c<-eta block (coupled_variables = eta), and
# the mobility is the degenerate M0 c (1-c), so SplitCHWRes carries dM/dc in its chem_pot_c<-c block
# (coupled_variables = c). eta itself relaxes (a time derivative and a unit-rate reaction) from a smooth
# profile, so every cross term is nonzero and non-uniform at every Newton iteration of the Jacobian gate.
# Scaling: M0 = 1 (M(0.5) = 0.25) and dt = 0.1 s, so that the chemical potential row's blocks are O(1) like the
# concentration row's and a zeroed mobility-derivative term is visible in the Frobenius ratio of the gate.
# Same mesh, initial condition for c, boundaries (natural) and solver as ch_rev_cpu.i; two
# DerivativeParsedMaterial blocks (unclamped) supply the free energy and the mobility with their derivatives.
# The Kokkos twin ../kokkos/ch_rev_coupled_kokkos.i differs only by the Kokkos prefix and the two material blocks.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 16
  ny = 16
  xmax = 10
  ymax = 10
[]

[Variables]
  [c]
    order = FIRST
    family = LAGRANGE
  []
  [chem_pot_c]
    order = FIRST
    family = LAGRANGE
  []
  [eta]
    order = FIRST
    family = LAGRANGE
  []
[]

[ICs]
  [c_ic]
    type = FunctionIC
    variable = c
    function = '0.5 + 0.04*cos(2*pi*x/10)*cos(2*pi*y/10) + 0.02*cos(4*pi*x/10)'
  []
  [eta_ic]
    type = FunctionIC
    variable = eta
    function = '0.5 + 0.2*cos(pi*x/10)*cos(pi*y/10)'
  []
[]

# [BCs] block removed: no periodic pairing, natural (zero-flux) Neumann on all four edges

[Kernels]
  [c_CoupledTimeDerivative]
    type = CoupledTimeDerivative
    variable = chem_pot_c
    v = c
  []
  [c_SplitCHWRes]
    type = SplitCHWRes
    variable = chem_pot_c
    mob_name = M
    coupled_variables = c
  []
  [c_SplitCHParsed]
    type = SplitCHParsed
    variable = c
    w = chem_pot_c
    f_name = f_loc
    kappa_name = kappa_c
    coupled_variables = eta
  []
  [eta_TimeDerivative]
    type = TimeDerivative
    variable = eta
  []
  [eta_Reaction]
    type = Reaction
    variable = eta
    rate = 1
  []
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'kappa_c'
    prop_values = '${fparse 5e-16*6.24150934e+18*1e-9/7.1e-6}'
  []
  [mobility]
    # M = M0 c (1-c) with dM/dc from the parser; M0 = 1 so that the chemical potential row weighs in the Jacobian norm (a Jacobian gate, not the tutorial's mobility)
    type = DerivativeParsedMaterial
    property_name = M
    coupled_variables = c
    constant_names = 'M0'
    constant_expressions = '1'
    expression = 'M0*c*(1-c)'
  []
  [local_energy]
    # W (c-ca)^2 (c-cb)^2 + a c eta^2 with its derivatives from the parser, including d^2F/dc deta
    type = DerivativeParsedMaterial
    property_name = f_loc
    coupled_variables = 'c eta'
    constant_names = 'W ca cb a'
    constant_expressions = '16 0.27 0.83 1'
    expression = 'W*(c-ca)^2*(c-cb)^2 + a*c*eta^2'
  []
[]

[Postprocessors]
  # int f_loc + kappa_c/2 int |grad c|^2 (not a Lyapunov functional of this artificial coupled system: eta
  # relaxes on its own; recorded as the same scalar on both backends, not as a monotonicity check)
  [E_bulk]
    type = ElementIntegralMaterialProperty
    mat_prop = f_loc
    execute_on = 'initial timestep_end'
  []
  [grad_c_norm]
    type = ElementH1SemiError
    variable = c
    function = 0
    execute_on = 'initial timestep_end'
  []
  [total_energy]
    type = ParsedPostprocessor
    expression = 'E_bulk + 0.5*kappa_c*grad_c_norm^2'
    pp_names = 'E_bulk grad_c_norm'
    constant_names = 'kappa_c'
    constant_expressions = '${fparse 5e-16*6.24150934e+18*1e-9/7.1e-6}'
    execute_on = 'initial timestep_end'
  []
  [c_avg]
    type = ElementAverageValue
    variable = c
    execute_on = 'initial timestep_end'
  []
  [eta_avg]
    type = ElementAverageValue
    variable = eta
    execute_on = 'initial timestep_end'
  []
  [max_c]
    type = ElementExtremeValue
    variable = c
    value_type = max
    execute_on = 'initial timestep_end'
  []
  [min_c]
    type = ElementExtremeValue
    variable = c
    value_type = min
    execute_on = 'initial timestep_end'
  []
[]

[Preconditioning]
  [coupled]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = none
  # explicit: a deck with Kokkos kernels aborts under automatic scaling at the pin (the scaling Jacobian is a
  # libMesh DiagonalMatrix, NonlinearSystemBase.C:355, while the Kokkos assembly requires a PetscMatrix,
  # KokkosMatrix.K:21); the CPU control carries the same line so the decks stay identical
  automatic_scaling = false
  petsc_options_iname = '-pc_type -sub_pc_type -ksp_gmres_restart'
  petsc_options_value = 'asm lu 31'
  l_tol = 1e-8
  l_max_its = 50
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  nl_max_its = 25
  dt = 0.1
  num_steps = 5
[]

[Outputs]
  exodus = true
  csv = true
[]
