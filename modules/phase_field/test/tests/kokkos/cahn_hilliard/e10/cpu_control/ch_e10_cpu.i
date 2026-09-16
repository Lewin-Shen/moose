# Port 02 tier 5: CPU control of the E10 gate deck, the fully implicit split Cahn-Hilliard in (g, u, c)
# variables (implicit Euler for the forward split with S added and subtracted: g-row M g - F'(c) + S M c = 0,
# u-row (M + a1 K) u - M c_n + dt L K g = 0, c-row (M + a2 K) c - M u = 0; a1 + a2 = dt L S, a1 a2 = dt L kappa;
# S = 18.3 sets the preconditioning operator only; dt = 2000 s is baked into alpha1, alpha2, dtL).
# Derived from Campaign Runs/CPU/02-ch-kks-solvers/decks/fimp2d_fspmboomergpunocf_dt02000_S18.3.i (read-only,
# written by scripts/campaign.py on 2026-09-14) with these changes:
#   1. de-periodized: the [BCs]/[Periodic] block removed (natural zero-flux boundaries, D-022);
#   2. the SolutionUserObject/SolutionIC restart replaced by a deterministic function initial condition for
#      c inside (0, 1) (decision I-2); u and g start at zero as in the campaign (their rows are algebraic);
#   3. the mesh reduced to 20 x 20 elements on the 25 nm box and five steps of dt = 2000 s;
#   4. the hypre relax types on the u and c blocks set to Chebyshev (Campaign 03's E10 smoother, design 6.4)
#      and the nleqerr line search passed as a PETSc option, as Campaign 03 did on its command line (it is
#      not in MOOSE's line_search enum at this ref);
#   5. TotalFreeEnergy, f_density and total_energy removed (no Kokkos form of TotalFreeEnergy exists; the
#      pair keeps identical variable sets for the exodiff gate); c_integral is kept for the mass drift check;
#   6. exodus output added for the exodiff gate;
#   7. automatic_scaling = false (the campaign deck has true): the Kokkos twin cannot use automatic scaling
#      at the pin (KokkosMatrix.K:21 requires a PetscMatrix; the scaling matrix is a DiagonalMatrix), and the
#      pair must differ only by the Kokkos prefix. Newton counts and c statistics are unchanged (tier 5 review);
#   8. the kappa (kappa_c) and mobility (M) materials removed: unused in both decks (the only consumer of kappa_c
#      was the removed TotalFreeEnergy; M is unused in the campaign deck as well, every coefficient is baked into
#      alpha1, alpha2, dtL) (tier 5 review, L4).
# Kokkos twin: ../kokkos/ch_e10_kokkos.i (Kokkos object names; the free-energy and constant materials are
# KokkosCHFreeEnergy and KokkosGenericConstantMaterial; the u_cold term is KokkosOldCoupledForce; the g_s and
# c_u terms are KokkosMatCoupledForce with the constant property 'one').

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
    type = ParsedAux
    variable = precipitate_indicator
    coupled_variables = c
    expression = 'if(c>0.4,1.0,0)'
    execute_on = 'initial TIMESTEP_END'
  []
[]

[Materials]
  [local_energy]
    type = DerivativeParsedMaterial
    property_name = f_loc
    coupled_variables = c
    constant_names = 'A   B   C   D   E   F   G  length_scale  eVpJ  Vm'
    constant_expressions = '-2.45e+04 -2.83e+04 4.17e+03 7.05e+03
                            1.21e+04 2.57e+03 -2.35e+03
                            1e-9 6.24150934e+18 7.1e-6'
    expression = 'eVpJ/Vm*length_scale^3*(A*c+B*(1-c)+C*c*log(c)+D*(1-c)*log(1-c)+
                E*c*(1-c)+F*c*(1-c)*(2*c-1)+G*c*(1-c)*(2*c-1)^2)'
  []
  [energy_neg]
    type = DerivativeParsedMaterial
    property_name = f_neg
    coupled_variables = c
    constant_names = 'A   B   C   D   E   F   G  length_scale  eVpJ  Vm'
    constant_expressions = '-2.45e+04 -2.83e+04 4.17e+03 7.05e+03
                            1.21e+04 2.57e+03 -2.35e+03
                            1e-9 6.24150934e+18 7.1e-6'
    expression = '-eVpJ/Vm*length_scale^3*(A*c+B*(1-c)+C*c*log(c)+D*(1-c)*log(1-c)+
                E*c*(1-c)+F*c*(1-c)*(2*c-1)+G*c*(1-c)*(2*c-1)^2)'
  []
  [fimp_coefficients]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha2 dtL'
    prop_values = '0.780418570089 0.0247814299108 0.044'
  []
[]

[Kernels]
  [g_mass]
    type = Reaction
    variable = g
    rate = 1
  []
  [g_f]
    type = CoupledMaterialDerivative
    variable = g
    v = c
    f_name = f_neg
  []
  [g_s]
    type = CoupledForce
    variable = g
    v = c
    coef = -18.3
  []
  [u_mass]
    type = Reaction
    variable = u
    rate = 1
  []
  [u_stiff]
    type = MatDiffusion
    variable = u
    diffusivity = alpha1
  []
  [u_cold]
    type = CoupledForce
    variable = u
    v = c
    coef = 1
    implicit = false
  []
  [u_g]
    type = MatDiffusion
    variable = u
    v = g
    diffusivity = dtL
  []
  [c_mass]
    type = Reaction
    variable = c
    rate = 1
  []
  [c_stiff]
    type = MatDiffusion
    variable = c
    diffusivity = alpha2
  []
  [c_u]
    type = CoupledForce
    variable = c
    v = u
    coef = 1
  []
[]

[Postprocessors]
  [volume_fraction]
    type = ElementAverageValue
    variable = precipitate_indicator
    execute_on = 'initial timestep_end'
  []
  [max_concentration]
    type = ElementExtremeValue
    variable = c
    value_type = max
    execute_on = 'initial timestep_end'
  []
  [min_concentration]
    type = ElementExtremeValue
    variable = c
    value_type = min
    execute_on = 'initial timestep_end'
  []
  [c_integral]
    type = ElementIntegralVariablePostprocessor
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
