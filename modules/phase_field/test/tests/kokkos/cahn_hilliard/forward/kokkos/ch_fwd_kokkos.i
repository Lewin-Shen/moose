# Kokkos twin of the tutorial's forward-split Cahn-Hilliard deck (Default f_loc/spinodal.i), the port 02
# tier 2 gate deck with the REDLICH_KISTER (tutorial CALPHAD) free energy.
#
# Changes against the tutorial deck: [Modules/PhaseField/Conserved] solve_type = FORWARD_SPLIT is written
# out as the five kernel blocks ConservedAction adds (ConservedAction.C:215-282 @ 6956830688: TimeDerivative
# on c; MatDiffusion on c with v = w, diffusivity = M; MatDiffusion on w with v = c, diffusivity = kappa_c;
# CoupledMaterialDerivative on w with v = c; CoefReaction on w with coefficient = -1) with their Kokkos
# twins (KokkosReaction rate = -1 is CoefReaction coefficient = -1); the free energy through
# KokkosCHFreeEnergy (closed form, tutorial defaults) instead of DerivativeParsedMaterial; kappa_c and M
# through KokkosGenericConstantMaterial; no [BCs/Periodic] (D-022, natural boundaries); a deterministic
# function initial condition instead of RandomIC (I-2); a 16x16 mesh (I-2); NEWTON with SMP full = true and a
# direct solve (the solver is not the point of the gates); no adaptive stepping, no predictor; automatic_scaling = false
# explicitly (a deck with Kokkos kernels aborts under automatic scaling at the pin).
# The CPU control ../cpu_control/ch_fwd_cpu.i differs by the Kokkos prefix and the free-energy block only.
#
# Field scales for the exodiff floor: c in [0.52, 0.56] (order 1), w = f'(c) - kappa laplacian(c) of order 4.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 16
  ny = 16
  xmax = 25
  ymax = 25
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
    expression = '0.54 + 0.02*cos(4*pi*x/25)*cos(4*pi*y/25)'
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
  # ConservedAction FORWARD_SPLIT expansion, ConservedAction.C:215-282 @ 6956830688
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
    # spinodal.i: kappa_c = 5e-16 * eV_J * nm_m^2 * d, evaluated once by the parser on both backends
    type = KokkosGenericConstantMaterial
    prop_names = 'kappa_c'
    prop_values = '${fparse 5e-16*6.24150934e+18*1e-9/7.1e-6}'
  []
  [mobility]
    type = KokkosGenericConstantMaterial
    prop_names = 'M'
    prop_values = '2.2e-5'
  []
  [local_energy]
    # The tutorial's CALPHAD f_loc in closed form; every coefficient defaults to spinodal.i:90-95
    type = KokkosCHFreeEnergy
    property_name = f_loc
    c = c
    form = REDLICH_KISTER
  []
  [energy_density]
    # total free energy density f_loc + kappa_c / 2 |grad(c)|^2 as a material property, integrated at the
    # quadrature points by KokkosElementIntegralMaterialProperty on the device (a host object over a
    # Kokkos property reads a silent zero, so the energy of the Kokkos deck is formed on the device)
    type = KokkosCHEnergyDensity
    f_name = f_loc
    kappa_name = kappa_c
    c = c
  []
[]

[Postprocessors]
  [c_total]
    type = KokkosElementIntegralVariablePostprocessor
    variable = c
    execute_on = 'initial timestep_end'
  []
  [max_c]
    type = KokkosElementExtremeValue
    variable = c
    value_type = max
    execute_on = 'initial timestep_end'
  []
  [min_c]
    type = KokkosElementExtremeValue
    variable = c
    value_type = min
    execute_on = 'initial timestep_end'
  []
  [total_energy]
    type = KokkosElementIntegralMaterialProperty
    mat_prop = f_total
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
  # the relative tolerance is measured against each step's initial residual, which shrinks as the
  # dynamics slow; the absolute tolerance is the effective criterion on this O(100) residual scale
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-11
  nl_max_its = 20
  l_max_its = 30
  dt = 2000
  num_steps = 5
[]

[Outputs]
  exodus = true
  csv = true
[]
