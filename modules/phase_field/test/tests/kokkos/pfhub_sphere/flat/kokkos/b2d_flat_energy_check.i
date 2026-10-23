# Energy cross-check of port 05 (gate G5 (i)): the Kokkos twin b2d_flat_kokkos.i unchanged, plus the
# same total energy built a second way. Five single-variable KokkosCHEnergyDensity instances (c with
# f_name = f_loc and kappa_c; each eta with f_name = 0 and kappa_eta; each its own property_name, since the
# default f_total would be declared five times), five KokkosElementIntegralMaterialProperty, and their
# sum free_energy_sum, beside the list material's integral free_energy. Acceptance, fixed before the run:
# abs(free_energy - free_energy_sum) at most 1e-12 times the sum of the five absolute integrals (the
# energy is not sign-definite where h exceeds 1), at every step; the Terminator fails the run otherwise.
# Because every kappa of the deck is 3, this check cannot see a kappa-to-variable mispairing; the
# pairing fixture b2d_density_pairing_kokkos.i (distinct kappas 1 to 5) does.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
  xmin = 0
  xmax = 25
  ymin = 0
  ymax = 25
  elem_type = QUAD4
[]

[Variables]
  [c]
  []
  [chem_pot_c]
  []
  [eta1]
  []
  [eta2]
  []
  [eta3]
  []
  [eta4]
  []
[]

[Functions]
  [c_ic]
    type = ParsedFunction
    expression = '0.4 + 0.05*cos(2*pi*x/25)*cos(2*pi*y/25)'
  []
  [eta1_ic]
    type = ParsedFunction
    expression = '0.2 + 0.05*sin(2*pi*(x + 1*y)/25)'
  []
  [eta2_ic]
    type = ParsedFunction
    expression = '0.35 + 0.05*sin(2*pi*(x + 2*y)/25)'
  []
  [eta3_ic]
    type = ParsedFunction
    expression = '0.6 + 0.05*sin(2*pi*(x + 3*y)/25)'
  []
  [eta4_ic]
    type = ParsedFunction
    expression = '0.85 + 0.05*sin(2*pi*(x + 4*y)/25)'
  []
[]

[ICs]
  [c_ic]
    type = FunctionIC
    variable = c
    function = c_ic
  []
  [eta1_ic]
    type = FunctionIC
    variable = eta1
    function = eta1_ic
  []
  [eta2_ic]
    type = FunctionIC
    variable = eta2
    function = eta2_ic
  []
  [eta3_ic]
    type = FunctionIC
    variable = eta3
    function = eta3_ic
  []
  [eta4_ic]
    type = FunctionIC
    variable = eta4
    function = eta4_ic
  []
[]

[Kernels]
  # ConservedAction, FORWARD_SPLIT, on c and chem_pot_c
  [c_TimeDerivative]
    type = KokkosTimeDerivative
    variable = c
  []
  [c_MatDiffusion]
    type = KokkosMatDiffusion
    variable = c
    v = chem_pot_c
    diffusivity = M
  []
  [chem_pot_c_MatDiffusion]
    type = KokkosMatDiffusion
    variable = chem_pot_c
    v = c
    diffusivity = kappa_c
  []
  [chem_pot_c_CoupledMaterialDerivative]
    type = KokkosCoupledMaterialDerivative
    variable = chem_pot_c
    v = c
    f_name = f_loc
    coupled_variables = 'eta1 eta2 eta3 eta4'
  []
  [chem_pot_c_CoefReaction]
    type = KokkosReaction
    variable = chem_pot_c
    rate = -1
  []

  # NonconservedAction on eta1 to eta4
  [eta1_TimeDerivative]
    type = KokkosTimeDerivative
    variable = eta1
  []
  [eta1_AllenCahn]
    type = KokkosAllenCahn
    variable = eta1
    mob_name = L
    f_name = f_loc
    coupled_variables = 'c eta2 eta3 eta4'
  []
  [eta1_ACInterface]
    type = KokkosACInterface
    variable = eta1
    mob_name = L
    kappa_name = kappa_eta
    variable_L = true
  []

  [eta2_TimeDerivative]
    type = KokkosTimeDerivative
    variable = eta2
  []
  [eta2_AllenCahn]
    type = KokkosAllenCahn
    variable = eta2
    mob_name = L
    f_name = f_loc
    coupled_variables = 'c eta1 eta3 eta4'
  []
  [eta2_ACInterface]
    type = KokkosACInterface
    variable = eta2
    mob_name = L
    kappa_name = kappa_eta
    variable_L = true
  []

  [eta3_TimeDerivative]
    type = KokkosTimeDerivative
    variable = eta3
  []
  [eta3_AllenCahn]
    type = KokkosAllenCahn
    variable = eta3
    mob_name = L
    f_name = f_loc
    coupled_variables = 'c eta1 eta2 eta4'
  []
  [eta3_ACInterface]
    type = KokkosACInterface
    variable = eta3
    mob_name = L
    kappa_name = kappa_eta
    variable_L = true
  []

  [eta4_TimeDerivative]
    type = KokkosTimeDerivative
    variable = eta4
  []
  [eta4_AllenCahn]
    type = KokkosAllenCahn
    variable = eta4
    mob_name = L
    f_name = f_loc
    coupled_variables = 'c eta1 eta2 eta3'
  []
  [eta4_ACInterface]
    type = KokkosACInterface
    variable = eta4
    mob_name = L
    kappa_name = kappa_eta
    variable_L = true
  []
[]

[AuxVariables]
  [f_density]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [f_density]
    type = KokkosMaterialRealAux
    variable = f_density
    property = f_density
    execute_on = 'initial timestep_end'
  []
[]

[Materials]
  [constants_c]
    type = KokkosGenericConstantMaterial
    prop_names = 'M kappa_c'
    prop_values = '5 3'
  []
  [constants_eta]
    type = KokkosGenericConstantMaterial
    prop_names = 'L kappa_eta'
    prop_values = '5 3'
  []
  [local_energy]
    type = KokkosTwoPhaseParabolicFreeEnergy
    property_name = f_loc
    c = c
    etas = 'eta1 eta2 eta3 eta4'
  []
  [energy_density]
    type = KokkosTotalFreeEnergyDensity
    f_name = f_loc
    interfacial_vars = 'c eta1 eta2 eta3 eta4'
    kappa_names = 'kappa_c kappa_eta kappa_eta kappa_eta kappa_eta'
    property_name = f_density
  []
  [energy_c]
    type = KokkosCHEnergyDensity
    c = c
    f_name = f_loc
    kappa_name = kappa_c
    property_name = f_total_c
  []
  [energy_eta1]
    type = KokkosCHEnergyDensity
    c = eta1
    f_name = 0
    kappa_name = kappa_eta
    property_name = f_total_eta1
  []
  [energy_eta2]
    type = KokkosCHEnergyDensity
    c = eta2
    f_name = 0
    kappa_name = kappa_eta
    property_name = f_total_eta2
  []
  [energy_eta3]
    type = KokkosCHEnergyDensity
    c = eta3
    f_name = 0
    kappa_name = kappa_eta
    property_name = f_total_eta3
  []
  [energy_eta4]
    type = KokkosCHEnergyDensity
    c = eta4
    f_name = 0
    kappa_name = kappa_eta
    property_name = f_total_eta4
  []
[]

[Postprocessors]
  [free_energy]
    type = KokkosElementIntegralMaterialProperty
    mat_prop = f_density
    execute_on = 'initial timestep_end'
  []
  [e_c]
    type = KokkosElementIntegralMaterialProperty
    mat_prop = f_total_c
    execute_on = 'initial timestep_end'
  []
  [e_eta1]
    type = KokkosElementIntegralMaterialProperty
    mat_prop = f_total_eta1
    execute_on = 'initial timestep_end'
  []
  [e_eta2]
    type = KokkosElementIntegralMaterialProperty
    mat_prop = f_total_eta2
    execute_on = 'initial timestep_end'
  []
  [e_eta3]
    type = KokkosElementIntegralMaterialProperty
    mat_prop = f_total_eta3
    execute_on = 'initial timestep_end'
  []
  [e_eta4]
    type = KokkosElementIntegralMaterialProperty
    mat_prop = f_total_eta4
    execute_on = 'initial timestep_end'
  []
  [free_energy_sum]
    type = ParsedPostprocessor
    pp_names = 'e_c e_eta1 e_eta2 e_eta3 e_eta4'
    expression = 'e_c + e_eta1 + e_eta2 + e_eta3 + e_eta4'
    execute_on = 'initial timestep_end'
  []
  [energy_check_rel]
    type = ParsedPostprocessor
    pp_names = 'free_energy free_energy_sum e_c e_eta1 e_eta2 e_eta3 e_eta4'
    expression = 'abs(free_energy - free_energy_sum) / (abs(e_c) + abs(e_eta1) + abs(e_eta2) + abs(e_eta3) + abs(e_eta4))'
    execute_on = 'initial timestep_end'
  []
  # the acceptance as a 0/1 flag: ParsedPostprocessor compares with epsilon 0 (FunctionParserUtils),
  # a Terminator compares with fparser's default epsilon (about 1e-12), so the comparison lives here
  [energy_check_fail]
    type = ParsedPostprocessor
    pp_names = 'energy_check_rel'
    expression = 'if(energy_check_rel > 1e-12, 1, 0)'
    execute_on = 'initial timestep_end'
  []
[]

[UserObjects]
  [energy_check]
    type = Terminator
    expression = 'energy_check_fail > 0.5'
    fail_mode = HARD
    error_level = ERROR
    message = 'G5 (i): the list material integral differs from the five-instance sum by more than 1e-12 of the summed absolute integrals'
    execute_on = 'initial timestep_end'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  automatic_scaling = false
  petsc_options_iname = '-pc_type -sub_pc_type -snes_linesearch_type'
  petsc_options_value = 'asm lu basic'
  l_tol = 1e-8
  l_max_its = 100
  nl_rel_tol = 1e-11
  nl_abs_tol = 1e-12
  nl_max_its = 20
  num_steps = 5
  [TimeStepper]
    type = ConstantDT
    dt = 0.01
  []
[]

[Outputs]
  csv = true
  execute_on = 'initial timestep_end'
  [exodus]
    type = Exodus
    execute_input_on = none
  []
[]
