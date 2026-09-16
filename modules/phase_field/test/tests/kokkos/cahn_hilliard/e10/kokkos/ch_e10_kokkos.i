# Port 02 tier 5: Kokkos E10 gate deck, the fully implicit split Cahn-Hilliard in (g, u, c) variables
# (implicit Euler for the forward split with S added and subtracted: g-row M g - F'(c) + S M c = 0,
# u-row (M + a1 K) u - M c_n + dt L K g = 0, c-row (M + a2 K) c - M u = 0; a1 + a2 = dt L S, a1 a2 = dt L kappa;
# S = 18.3 sets the preconditioning operator only; dt = 2000 s is baked into alpha1, alpha2, dtL).
# Twin of ../cpu_control/ch_e10_cpu.i (itself derived from Campaign 02's
# fimp2d_fspmboomergpunocf_dt02000_S18.3.i, de-periodized, function IC, 20 x 20 mesh, Chebyshev, nleqerr;
# the change list is in that file). Differences from the CPU control:
#   1. Kokkos object names for every kernel, aux kernel and element postprocessor (TimestepSize,
#      the iteration-count postprocessors were removed at integration: rank-dependent, and design trap T18 keeps them out of Exodus);
#   2. the u_cold term is KokkosOldCoupledForce (the old-state read is built in; no implicit parameter);
#   3. the g_s and c_u terms are KokkosMatCoupledForce with material_properties = one, a constant property
#      of value 1 from KokkosGenericConstantMaterial, and coef carrying the coefficient (Codex H04: the
#      kernel creates its property array only when material_properties is given, and the dispatcher's copy
#      of an unallocated array aborts);
#   4. the free energies f_loc and f_neg are two instances of tier 2's KokkosCHFreeEnergy in the
#      REDLICH_KISTER form with the campaign's constants; f_neg negates eVpJ, which negates the prefactor
#      eVpJ/Vm*length_scale^3 and therefore the whole expression, as the campaign's energy_neg does;
#   5. the kappa (kappa_c) and mobility (M) materials are absent from both decks: unused (the only consumer of
#      kappa_c was the removed TotalFreeEnergy; M is unused in the campaign deck as well, every coefficient is
#      baked into alpha1, alpha2, dtL) (tier 5 review, L4);
#   6. automatic_scaling is off in both decks: at the pin the scaling Jacobian is assembled into a libMesh
#      DiagonalMatrix (NonlinearSystemBase.C:355) and the Kokkos assembly requires a PetscMatrix
#      (KokkosMatrix.K:21), so any deck with Kokkos kernels aborts in preSolve (tier 5 review, 2026-09-16).
# UNTESTED UNTIL INTEGRATION: KokkosMatDiffusion, KokkosCoupledMaterialDerivative and KokkosCHFreeEnergy
# are tier 2 objects built in another worktree; the KokkosCHFreeEnergy parameter names below follow the
# tier 2 brief and must be checked against the staged header at integration.

[Mesh]
  type = GeneratedMesh
  dim = 2
  xmax = 25
  ymax = 25
  nx = 20
  ny = 20
[]

[Variables]
  [c]
  []
  [u]
  []
  [g]
  []
[]

[Functions]
  [c_ic]
    type = ParsedFunction
    expression = '0.5 + 0.04*cos(2*pi*x/25)*cos(2*pi*y/25) + 0.02*cos(pi*x/25)'
  []
[]

[ICs]
  [c_ic]
    type = FunctionIC
    variable = c
    function = c_ic
  []
[]

[AuxVariables]
  [precipitate_indicator]
  []
[]

[AuxKernels]
  [precipitate_indicator]
    type = KokkosParsedAux
    variable = precipitate_indicator
    coupled_variables = c
    expression = 'if(c>0.4,1.0,0)'
    execute_on = 'initial TIMESTEP_END'
  []
[]

[Materials]
  [local_energy]
    type = KokkosCHFreeEnergy
    property_name = f_loc
    c = c
    form = REDLICH_KISTER
    A = -2.45e+04
    B = -2.83e+04
    C = 4.17e+03
    D = 7.05e+03
    E = 1.21e+04
    F = 2.57e+03
    G = -2.35e+03
    length_scale = 1e-9
    eVpJ = 6.24150934e+18
    Vm = 7.1e-6
  []
  [energy_neg]
    type = KokkosCHFreeEnergy
    property_name = f_neg
    c = c
    form = REDLICH_KISTER
    A = -2.45e+04
    B = -2.83e+04
    C = 4.17e+03
    D = 7.05e+03
    E = 1.21e+04
    F = 2.57e+03
    G = -2.35e+03
    length_scale = 1e-9
    eVpJ = -6.24150934e+18
    Vm = 7.1e-6
  []
  [fimp_coefficients]
    type = KokkosGenericConstantMaterial
    prop_names = 'alpha1 alpha2 dtL'
    prop_values = '0.780418570089 0.0247814299108 0.044'
  []
  [one]
    type = KokkosGenericConstantMaterial
    prop_names = 'one'
    prop_values = '1'
  []
[]

[Kernels]
  [g_mass]
    type = KokkosReaction
    variable = g
    rate = 1
  []
  [g_f]
    type = KokkosCoupledMaterialDerivative
    variable = g
    v = c
    f_name = f_neg
  []
  [g_s]
    type = KokkosMatCoupledForce
    variable = g
    v = c
    coef = -18.3
    material_properties = one
  []
  [u_mass]
    type = KokkosReaction
    variable = u
    rate = 1
  []
  [u_stiff]
    type = KokkosMatDiffusion
    variable = u
    diffusivity = alpha1
  []
  [u_cold]
    type = KokkosOldCoupledForce
    variable = u
    v = c
    coef = 1
  []
  [u_g]
    type = KokkosMatDiffusion
    variable = u
    v = g
    diffusivity = dtL
  []
  [c_mass]
    type = KokkosReaction
    variable = c
    rate = 1
  []
  [c_stiff]
    type = KokkosMatDiffusion
    variable = c
    diffusivity = alpha2
  []
  [c_u]
    type = KokkosMatCoupledForce
    variable = c
    v = u
    coef = 1
    material_properties = one
  []
[]

[Postprocessors]
  [volume_fraction]
    type = KokkosElementAverageValue
    variable = precipitate_indicator
    execute_on = 'initial timestep_end'
  []
  [max_concentration]
    type = KokkosElementExtremeValue
    variable = c
    value_type = max
    execute_on = 'initial timestep_end'
  []
  [min_concentration]
    type = KokkosElementExtremeValue
    variable = c
    value_type = min
    execute_on = 'initial timestep_end'
  []
  [c_integral]
    type = KokkosElementIntegralVariablePostprocessor
    variable = c
    execute_on = 'initial timestep_end'
  []
  [dt]
    type = TimestepSize
  []
[]

[Preconditioning]
  [fsp]
    type = FSP
    full = true
    topsplit = 'guc'
    [guc]
      splitting = 'g u c'
      splitting_type = multiplicative
    []
    [g]
      vars = 'g'
      petsc_options_iname = '-ksp_type -pc_type -pc_jacobi_fixdiagonal'
      petsc_options_value = 'preonly jacobi false'
    []
    [u]
      vars = 'u'
      petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type -pc_hypre_boomeramg_relax_type_down -pc_hypre_boomeramg_relax_type_up -pc_hypre_boomeramg_no_CF'
      petsc_options_value = 'preonly hypre boomeramg PMIS ext+i Chebyshev Chebyshev true'
    []
    [c]
      vars = 'c'
      petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type -pc_hypre_boomeramg_relax_type_down -pc_hypre_boomeramg_relax_type_up -pc_hypre_boomeramg_no_CF'
      petsc_options_value = 'preonly hypre boomeramg PMIS ext+i Chebyshev Chebyshev true'
    []
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = default
  petsc_options_iname = '-snes_linesearch_type'
  petsc_options_value = 'nleqerr'
  l_tol = 1e-4
  nl_abs_tol = 1e-9
  nl_rel_tol = 1e-12
  nl_max_its = 20
  l_max_its = 1000
  num_steps = 5
  automatic_scaling = false
  [TimeStepper]
    type = ConstantDT
    dt = 2000
  []
[]

[Outputs]
  exodus = true
  csv = true
[]
