# Forward-split Cahn-Hilliard with a degenerate mobility M = M0 c (1 - c) from KokkosCHMobility
# (implementation decisions I-3 and I-11): the Jacobian gate for the dD/du term of KokkosMatDiffusion
# (variable c, v = w, diffusivity = M with dM/dc declared) and for KokkosCHMobility's dM/dc.
# Unit-scaled coefficients (W = 1, kappa_c = 0.01, M0 = 1 on a unit box, dt = 0.05) and a 0.2 perturbation,
# so that every Jacobian term is within a few decades of the largest and the finite-difference ratio test
# can see each one (with the tutorial's 2.2e-5 mobility the dM/dc term would sit below the ratio tolerance).
# Kernel blocks as in ch_fwd_kokkos.i (ConservedAction.C:215-282 @ 6956830688). CPU control:
# ../cpu_control/ch_fwd_cpu_degenerate.i (DerivativeParsedMaterial mobility), differing by the prefixes and
# the two material blocks only.
#
# Field scales for the exodiff floor: c in [0.3, 0.7], w of order 0.1 to 1.

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
[]

[Functions]
  [c_ic]
    type = ParsedFunction
    expression = '0.5 + 0.2*cos(2*pi*x)*cos(2*pi*y)'
  []
[]

[ICs]
  [c_ic]
    type = FunctionIC
    variable = c
    function = c_ic
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
    f_name = f_loc
  []
  [w_CoefReaction]
    type = KokkosReaction
    variable = w
    rate = -1
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
  [local_energy]
    type = KokkosCHFreeEnergy
    property_name = f_loc
    c = c
    form = DOUBLE_WELL
    W = 1
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
