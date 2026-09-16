# Kokkos twin of the module's test/tests/phase_field_kernels/SplitCHWRes.i (two concentrations, off-diagonal Onsager
# mobilities): exercises the coupled-w branch of KokkosSplitCHWRes (w12: variable = w1, w = w2; w21: variable = w2,
# w = w1). The separable free energy 0.25 (1+c)^2 (1-c)^2 per component is two DOUBLE_WELL instances of
# KokkosCHFreeEnergy with W = 0.25, ca = -1, cb = 1. Implicit Euler in place of the module deck's bdf2; NEWTON,
# no line search, SMP full, lu, automatic_scaling = false, as the gate decks.
# Deck built by the tier 1 cross-review (port/reviews/FABLE-2026-09-16-t1-reverse-review.md, appendix A, L7).
# The CPU control ../cpu_control/ch_rev_onsager_cpu.i differs only by the Kokkos prefix and the two free-energy blocks.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  xmin = 0
  xmax = 60
  ymin = 0
  ymax = 60
  elem_type = QUAD4
[]

[Variables]
  [c1]
    [InitialCondition]
      type = FunctionIC
      function = 'cos(x/60*pi)'
    []
  []
  [c2]
    [InitialCondition]
      type = FunctionIC
      function = 'cos(y/60*pi)'
    []
  []
  [w1]
  []
  [w2]
  []
[]

[Kernels]
  [c1_res]
    type = KokkosSplitCHParsed
    variable = c1
    f_name = F1
    kappa_name = kappa_c
    w = w1
  []
  [w11_res]
    type = KokkosSplitCHWRes
    variable = w1
    mob_name = M11
  []
  [w12_res]
    type = KokkosSplitCHWRes
    variable = w1
    w = w2
    mob_name = M12
  []
  [c2_res]
    type = KokkosSplitCHParsed
    variable = c2
    f_name = F2
    kappa_name = kappa_c
    w = w2
  []
  [w22_res]
    type = KokkosSplitCHWRes
    variable = w2
    mob_name = M22
  []
  [w21_res]
    type = KokkosSplitCHWRes
    variable = w2
    w = w1
    mob_name = M21
  []
  [time1]
    type = KokkosCoupledTimeDerivative
    variable = w1
    v = c1
  []
  [time2]
    type = KokkosCoupledTimeDerivative
    variable = w2
    v = c2
  []
[]

[Materials]
  [pfmobility]
    type = KokkosGenericConstantMaterial
    prop_names = 'M11 M12 M21 M22 kappa_c'
    prop_values = '10  2.5 20  5   40'
  []
  [F1]
    type = KokkosCHFreeEnergy
    property_name = F1
    c = c1
    form = DOUBLE_WELL
    W = 0.25
    ca = -1
    cb = 1
  []
  [F2]
    type = KokkosCHFreeEnergy
    property_name = F2
    c = c2
    form = DOUBLE_WELL
    W = 0.25
    ca = -1
    cb = 1
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
  line_search = none
  # explicit: a deck with Kokkos kernels aborts under automatic scaling at the pin (the scaling Jacobian is a
  # libMesh DiagonalMatrix, NonlinearSystemBase.C:355, while the Kokkos assembly requires a PetscMatrix,
  # KokkosMatrix.K:21); the CPU control carries the same line so the decks stay identical
  automatic_scaling = false
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  l_max_its = 30
  l_tol = 1e-10
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  nl_max_its = 25
  start_time = 0.0
  num_steps = 2
  dt = 10
[]

[Outputs]
  exodus = true
[]
