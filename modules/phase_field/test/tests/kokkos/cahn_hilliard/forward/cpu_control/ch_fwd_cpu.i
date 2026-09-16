# Explicit CPU control of ../kokkos/ch_fwd_kokkos.i: the tutorial's forward-split Cahn-Hilliard deck
# (Default f_loc/spinodal.i) with [Modules/PhaseField/Conserved] solve_type = FORWARD_SPLIT written out as the
# five kernel blocks ConservedAction adds (ConservedAction.C:215-282 @ 6956830688; Reaction rate = -1 is
# CoefReaction coefficient = -1, which the action adds), the tutorial's DerivativeParsedMaterial f_loc verbatim
# (unclamped), kappa_c and M through GenericConstantMaterial; no [BCs/Periodic] (D-022); a deterministic
# function initial condition (I-2); a 16x16 mesh; NEWTON with SMP full = true and a direct solve.
# Differs from the Kokkos deck by the object prefixes and the free-energy material block only.
#
# Field scales for the exodiff floor: c in [0.52, 0.56] (order 1), w of order 4.

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
    type = TimeDerivative
    variable = c
  []
  [c_MatDiffusion]
    type = MatDiffusion
    variable = c
    v = w
    diffusivity = M
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
    f_name = f_loc
  []
  [w_CoefReaction]
    type = Reaction
    variable = w
    rate = -1
  []
[]

[Materials]
  [kappa]
    type = GenericConstantMaterial
    prop_names = 'kappa_c'
    prop_values = '${fparse 5e-16*6.24150934e+18*1e-9/7.1e-6}'
  []
  [mobility]
    type = GenericConstantMaterial
    prop_names = 'M'
    prop_values = '2.2e-5'
  []
  [local_energy]
    # spinodal.i:86-96 verbatim (the tutorial's CALPHAD f_loc, unclamped)
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
  [energy_density]
    # total free energy density f_loc + kappa_c / 2 |grad(c)|^2 as a material property, integrated at the
    # quadrature points by ElementIntegralMaterialProperty (a host object over a
    # Kokkos property reads a silent zero, so the energy of the Kokkos deck is formed on the device)
    type = CHEnergyDensity
    f_name = f_loc
    kappa_name = kappa_c
    c = c
  []
[]

[Postprocessors]
  [c_total]
    type = ElementIntegralVariablePostprocessor
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
  [total_energy]
    type = ElementIntegralMaterialProperty
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
