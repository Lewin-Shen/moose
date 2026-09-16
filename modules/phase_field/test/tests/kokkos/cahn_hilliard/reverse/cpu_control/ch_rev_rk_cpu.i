# Explicit CPU control of the reverse-split Cahn-Hilliard gate deck, CALPHAD (Redlich-Kister) energy (port 02, tier 1).
# Derived from the tutorial deck Default f_loc/spinodal.i (read-only) through the explicit
# REVERSE_SPLIT expansion of [Modules/PhaseField/Conserved] (ConservedAction.C:171-213 @ 6956830688:
# CoupledTimeDerivative on chem_pot_c with v = c; SplitCHWRes on chem_pot_c with mob_name; SplitCHParsed on c
# with w = chem_pot_c, f_name, kappa_name), with the CPU object names.
# Changes vs the tutorial deck: mesh 16x16 on 10x10 nm (decision I-2: small DOF count); a deterministic
# FunctionIC in place of RandomIC; [BCs/Periodic] REMOVED (constraints are silently ignored by the Kokkos
# assembly at the pin, D-022): natural zero-flux boundaries; the free energy through DerivativeParsedMaterial (unclamped, design section 5.5) with
# the tutorial deck's f_loc expression and its Fe-Cr coefficients, verbatim;
# NEWTON without line search, SMP full = true, asm + lu; dt = 500 s, five steps; automatic_scaling = false explicitly.
# The Kokkos twin ../kokkos/ch_rev_rk_kokkos.i differs only by the Kokkos prefix and the free-energy block.

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
    # f'(0.5) of the Redlich-Kister energy (3.69), so that the chemical potential does not start at zero: with w = 0
    # the Jacobian tester perturbs the w columns by 0.1 x step and the iteration-0 finite difference is roundoff
    # dominated on this energy (its f' is 15x the double well's); the converged solution does not depend on it
    initial_condition = 3.69
  []
[]

[ICs]
  [c_ic]
    type = FunctionIC
    variable = c
    function = '0.5 + 0.04*cos(2*pi*x/10)*cos(2*pi*y/10) + 0.02*cos(4*pi*x/10)'
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
  []
  [c_SplitCHParsed]
    type = SplitCHParsed
    variable = c
    w = chem_pot_c
    f_name = f_loc
    kappa_name = kappa_c
  []
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'M kappa_c'
    # M = 2.2e-5; kappa_c = 5e-16*6.24150934e+18*1e-9/7.1e-6 (the tutorial's expressions)
    prop_values = '2.2e-5 ${fparse 5e-16*6.24150934e+18*1e-9/7.1e-6}'
  []
  [local_energy]
    # the tutorial deck's [local_energy] block (Default f_loc/spinodal.i), verbatim
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
[]

[Postprocessors]
  # Total free energy = int f_loc + kappa_c/2 int |grad c|^2. A host TotalFreeEnergy aux kernel reads a
  # Kokkos material property as a silent zero (smoke run 2026-09-15), so the bulk part is integrated by the
  # property-integral postprocessor of the same backend and the gradient part comes from the variable alone.
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
  dt = 500
  num_steps = 5
[]

[Outputs]
  exodus = true
  csv = true
[]
