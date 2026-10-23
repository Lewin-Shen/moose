# CPU twin of the pairing fixture of port 05 (gate G5 (i-b)): TotalFreeEnergy with distinct kappas
# 1 to 5 on linear fields, c = 0.4 + 0.01 x and eta_k = b_k + 0.01 (k + 1) x, f_name = 0, no solve.
# The correct pairing of kappa_names[i] with interfacial_vars[i] gives 0.01125 at every point (a
# reversed pairing gives 0.00525). Acceptance, fixed before the run: the integral over the unit square
# and the largest and smallest element value of the f_density field equal 0.01125 within absolute
# 1e-12; the Terminator fails the run otherwise. The Kokkos twin ../kokkos/b2d_density_pairing_kokkos.i
# differs by the Kokkos prefix, the energy block, and the five single-variable instances it adds.

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
    type = TotalFreeEnergy
    variable = f_density
    f_name = 0
    kappa_names = 'kappa_1 kappa_2 kappa_3 kappa_4 kappa_5'
    interfacial_vars = 'c eta1 eta2 eta3 eta4'
    execute_on = 'initial timestep_end'
  []
[]

[Materials]
  [kappas]
    type = GenericConstantMaterial
    prop_names = 'kappa_1 kappa_2 kappa_3 kappa_4 kappa_5'
    prop_values = '1 2 3 4 5'
  []
[]

[Postprocessors]
  [free_energy]
    type = ElementIntegralVariablePostprocessor
    variable = f_density
    execute_on = 'initial timestep_end'
  []
  [f_density_max]
    type = ElementExtremeValue
    variable = f_density
    value_type = max
    execute_on = 'initial timestep_end'
  []
  [f_density_min]
    type = ElementExtremeValue
    variable = f_density
    value_type = min
    execute_on = 'initial timestep_end'
  []
[]

[UserObjects]
  [pairing_check]
    type = Terminator
    expression = 'abs(free_energy - 0.01125) > 1e-12 | abs(f_density_max - 0.01125) > 1e-12 | abs(f_density_min - 0.01125) > 1e-12'
    fail_mode = HARD
    error_level = ERROR
    message = 'G5 (i-b): the TotalFreeEnergy density or its integral differs from 0.01125 by more than 1e-12'
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
