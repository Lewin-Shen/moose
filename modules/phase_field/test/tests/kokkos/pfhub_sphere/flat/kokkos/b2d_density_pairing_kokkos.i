# Pairing fixture of port 05 (gate G5 (i-b)): KokkosTotalFreeEnergyDensity pairs kappa_names[i] with
# interfacial_vars[i] by position. The deck's kappas are all 3, so no check on the deck can see a
# mispairing; here the five kappas are distinct (1 to 5) and the fields are linear in x with distinct
# slopes, c = 0.4 + 0.01 x and eta_k = b_k + 0.01 (k + 1) x, so |grad v_i|^2 = (0.01 i)^2 and, with
# f_name = 0, the density is 0.5 (1 0.01^2 + 2 0.02^2 + 3 0.03^2 + 4 0.04^2 + 5 0.05^2) = 0.01125 at
# every point (a reversed pairing gives 0.00525). No solve. Acceptance, fixed before the run: the list
# material's integral over the unit square, the sum of the five single-variable KokkosCHEnergyDensity
# integrals, and the largest and smallest element value of the f_density field all equal 0.01125
# within absolute 1e-12; the Terminator fails the run otherwise. The CPU twin
# ../cpu_control/b2d_density_pairing_cpu.i applies the same acceptance to TotalFreeEnergy.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 4
  ny = 4
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  elem_type = QUAD4
[]

[Problem]
  solve = false
[]

[Variables]
  [c]
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
    expression = '0.4 + 0.01*x'
  []
  [eta1_ic]
    type = ParsedFunction
    expression = '0.2 + 0.02*x'
  []
  [eta2_ic]
    type = ParsedFunction
    expression = '0.35 + 0.03*x'
  []
  [eta3_ic]
    type = ParsedFunction
    expression = '0.6 + 0.04*x'
  []
  [eta4_ic]
    type = ParsedFunction
    expression = '0.85 + 0.05*x'
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
  [kappas]
    type = KokkosGenericConstantMaterial
    prop_names = 'kappa_1 kappa_2 kappa_3 kappa_4 kappa_5'
    prop_values = '1 2 3 4 5'
  []
  [energy_density]
    type = KokkosTotalFreeEnergyDensity
    f_name = 0
    interfacial_vars = 'c eta1 eta2 eta3 eta4'
    kappa_names = 'kappa_1 kappa_2 kappa_3 kappa_4 kappa_5'
    property_name = f_density
  []
  # the same pairs, one variable per instance
  [energy_c]
    type = KokkosCHEnergyDensity
    c = c
    f_name = 0
    kappa_name = kappa_1
    property_name = f_pair_c
  []
  [energy_eta1]
    type = KokkosCHEnergyDensity
    c = eta1
    f_name = 0
    kappa_name = kappa_2
    property_name = f_pair_eta1
  []
  [energy_eta2]
    type = KokkosCHEnergyDensity
    c = eta2
    f_name = 0
    kappa_name = kappa_3
    property_name = f_pair_eta2
  []
  [energy_eta3]
    type = KokkosCHEnergyDensity
    c = eta3
    f_name = 0
    kappa_name = kappa_4
    property_name = f_pair_eta3
  []
  [energy_eta4]
    type = KokkosCHEnergyDensity
    c = eta4
    f_name = 0
    kappa_name = kappa_5
    property_name = f_pair_eta4
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
    mat_prop = f_pair_c
    execute_on = 'initial timestep_end'
  []
  [e_eta1]
    type = KokkosElementIntegralMaterialProperty
    mat_prop = f_pair_eta1
    execute_on = 'initial timestep_end'
  []
  [e_eta2]
    type = KokkosElementIntegralMaterialProperty
    mat_prop = f_pair_eta2
    execute_on = 'initial timestep_end'
  []
  [e_eta3]
    type = KokkosElementIntegralMaterialProperty
    mat_prop = f_pair_eta3
    execute_on = 'initial timestep_end'
  []
  [e_eta4]
    type = KokkosElementIntegralMaterialProperty
    mat_prop = f_pair_eta4
    execute_on = 'initial timestep_end'
  []
  [free_energy_sum]
    type = ParsedPostprocessor
    pp_names = 'e_c e_eta1 e_eta2 e_eta3 e_eta4'
    expression = 'e_c + e_eta1 + e_eta2 + e_eta3 + e_eta4'
    execute_on = 'initial timestep_end'
  []
  [f_density_max]
    type = KokkosElementExtremeValue
    variable = f_density
    value_type = max
    execute_on = 'initial timestep_end'
  []
  [f_density_min]
    type = KokkosElementExtremeValue
    variable = f_density
    value_type = min
    execute_on = 'initial timestep_end'
  []
[]

[UserObjects]
  [pairing_check]
    type = Terminator
    expression = 'abs(free_energy - 0.01125) > 1e-12 | abs(free_energy_sum - 0.01125) > 1e-12 | abs(f_density_max - 0.01125) > 1e-12 | abs(f_density_min - 0.01125) > 1e-12'
    fail_mode = HARD
    error_level = ERROR
    message = 'G5 (i-b): a density, its integral or the five-instance sum differs from 0.01125 by more than 1e-12'
    execute_on = 'initial timestep_end'
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
  csv = true
  execute_on = 'initial timestep_end'
[]
