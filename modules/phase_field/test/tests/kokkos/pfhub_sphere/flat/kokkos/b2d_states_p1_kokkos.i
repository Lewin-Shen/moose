# Gate G2b of port 05: the 21 device outputs of KokkosTwoPhaseParabolicFreeEnergy at 17 fixed states
# (one state per element), written to CONSTANT MONOMIAL fields by 21 KokkosMaterialRealAux, no solve.
# Twin of ../cpu_control/b2d_states_p1_cpu.i (the DerivativeParsedMaterial of the benchmark deck at
# derivative_order = 2 and 21 MaterialRealAux); it differs by the Kokkos prefix and the free-energy
# block. Acceptance (design section 9): abs(device - host) < 1e-12 max(1, abs(device), abs(host)) on all
# 21 fields, exodiff -combined -t 1e-12 -F 0.
# The 17 states are copied verbatim from the first reader's fixed-state check (a grid of c times
# the cyclic shifts of four eta values, plus zeros, the wells, c = 0.5 and states outside [0, 1]);
# a 2D strip of 17 x 1 QUAD4, each state constant over its element (PiecewiseConstant with
# direction = left; the CONSTANT MONOMIAL initial condition is exact for a function constant on an element).

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 17
  ny = 1
  xmin = 0
  xmax = 17
  ymin = 0
  ymax = 1
  elem_type = QUAD4
[]

[Problem]
  solve = false
[]

[Variables]
  [c]
    family = MONOMIAL
    order = CONSTANT
  []
  [eta1]
    family = MONOMIAL
    order = CONSTANT
  []
  [eta2]
    family = MONOMIAL
    order = CONSTANT
  []
  [eta3]
    family = MONOMIAL
    order = CONSTANT
  []
  [eta4]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[Functions]
  [c_state]
    type = PiecewiseConstant
    axis = x
    direction = left
    x = '0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16'
    y = '0.2 0.2 0.2 0.2 0.45 0.45 0.45 0.45 0.65 0.65 0.65 0.65 0.5 0.3 0.7 -0.1 0.9'
  []
  [eta1_state]
    type = PiecewiseConstant
    axis = x
    direction = left
    x = '0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16'
    y = '0.1 0.3 0.6 0.9 0.1 0.3 0.6 0.9 0.1 0.3 0.6 0.9 0.0 1.0 0.25 1.2 0.9'
  []
  [eta2_state]
    type = PiecewiseConstant
    axis = x
    direction = left
    x = '0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16'
    y = '0.3 0.6 0.9 0.1 0.3 0.6 0.9 0.1 0.3 0.6 0.9 0.1 0.0 0.0 0.5 -0.2 0.1'
  []
  [eta3_state]
    type = PiecewiseConstant
    axis = x
    direction = left
    x = '0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16'
    y = '0.6 0.9 0.1 0.3 0.6 0.9 0.1 0.3 0.6 0.9 0.1 0.3 0.0 0.0 0.75 0.5 0.3'
  []
  [eta4_state]
    type = PiecewiseConstant
    axis = x
    direction = left
    x = '0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16'
    y = '0.9 0.1 0.3 0.6 0.9 0.1 0.3 0.6 0.9 0.1 0.3 0.6 0.0 0.5 0.05 0.0 0.6'
  []
[]

[ICs]
  [c_ic]
    type = FunctionIC
    variable = c
    function = c_state
  []
  [eta1_ic]
    type = FunctionIC
    variable = eta1
    function = eta1_state
  []
  [eta2_ic]
    type = FunctionIC
    variable = eta2
    function = eta2_state
  []
  [eta3_ic]
    type = FunctionIC
    variable = eta3
    function = eta3_state
  []
  [eta4_ic]
    type = FunctionIC
    variable = eta4
    function = eta4_state
  []
[]

[AuxVariables]
  [f]
    family = MONOMIAL
    order = CONSTANT
  []
  [d_c]
    family = MONOMIAL
    order = CONSTANT
  []
  [d2_c_c]
    family = MONOMIAL
    order = CONSTANT
  []
  [d_eta1]
    family = MONOMIAL
    order = CONSTANT
  []
  [d2_c_eta1]
    family = MONOMIAL
    order = CONSTANT
  []
  [d2_eta1_eta1]
    family = MONOMIAL
    order = CONSTANT
  []
  [d_eta2]
    family = MONOMIAL
    order = CONSTANT
  []
  [d2_c_eta2]
    family = MONOMIAL
    order = CONSTANT
  []
  [d2_eta2_eta2]
    family = MONOMIAL
    order = CONSTANT
  []
  [d_eta3]
    family = MONOMIAL
    order = CONSTANT
  []
  [d2_c_eta3]
    family = MONOMIAL
    order = CONSTANT
  []
  [d2_eta3_eta3]
    family = MONOMIAL
    order = CONSTANT
  []
  [d_eta4]
    family = MONOMIAL
    order = CONSTANT
  []
  [d2_c_eta4]
    family = MONOMIAL
    order = CONSTANT
  []
  [d2_eta4_eta4]
    family = MONOMIAL
    order = CONSTANT
  []
  [d2_eta1_eta2]
    family = MONOMIAL
    order = CONSTANT
  []
  [d2_eta1_eta3]
    family = MONOMIAL
    order = CONSTANT
  []
  [d2_eta1_eta4]
    family = MONOMIAL
    order = CONSTANT
  []
  [d2_eta2_eta3]
    family = MONOMIAL
    order = CONSTANT
  []
  [d2_eta2_eta4]
    family = MONOMIAL
    order = CONSTANT
  []
  [d2_eta3_eta4]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [f]
    type = KokkosMaterialRealAux
    variable = f
    property = 'f_loc'
    execute_on = 'initial timestep_end'
  []
  [d_c]
    type = KokkosMaterialRealAux
    variable = d_c
    property = 'df_loc/dc'
    execute_on = 'initial timestep_end'
  []
  [d2_c_c]
    type = KokkosMaterialRealAux
    variable = d2_c_c
    property = 'd^2f_loc/dc^2'
    execute_on = 'initial timestep_end'
  []
  [d_eta1]
    type = KokkosMaterialRealAux
    variable = d_eta1
    property = 'df_loc/deta1'
    execute_on = 'initial timestep_end'
  []
  [d2_c_eta1]
    type = KokkosMaterialRealAux
    variable = d2_c_eta1
    property = 'd^2f_loc/dcdeta1'
    execute_on = 'initial timestep_end'
  []
  [d2_eta1_eta1]
    type = KokkosMaterialRealAux
    variable = d2_eta1_eta1
    property = 'd^2f_loc/deta1^2'
    execute_on = 'initial timestep_end'
  []
  [d_eta2]
    type = KokkosMaterialRealAux
    variable = d_eta2
    property = 'df_loc/deta2'
    execute_on = 'initial timestep_end'
  []
  [d2_c_eta2]
    type = KokkosMaterialRealAux
    variable = d2_c_eta2
    property = 'd^2f_loc/dcdeta2'
    execute_on = 'initial timestep_end'
  []
  [d2_eta2_eta2]
    type = KokkosMaterialRealAux
    variable = d2_eta2_eta2
    property = 'd^2f_loc/deta2^2'
    execute_on = 'initial timestep_end'
  []
  [d_eta3]
    type = KokkosMaterialRealAux
    variable = d_eta3
    property = 'df_loc/deta3'
    execute_on = 'initial timestep_end'
  []
  [d2_c_eta3]
    type = KokkosMaterialRealAux
    variable = d2_c_eta3
    property = 'd^2f_loc/dcdeta3'
    execute_on = 'initial timestep_end'
  []
  [d2_eta3_eta3]
    type = KokkosMaterialRealAux
    variable = d2_eta3_eta3
    property = 'd^2f_loc/deta3^2'
    execute_on = 'initial timestep_end'
  []
  [d_eta4]
    type = KokkosMaterialRealAux
    variable = d_eta4
    property = 'df_loc/deta4'
    execute_on = 'initial timestep_end'
  []
  [d2_c_eta4]
    type = KokkosMaterialRealAux
    variable = d2_c_eta4
    property = 'd^2f_loc/dcdeta4'
    execute_on = 'initial timestep_end'
  []
  [d2_eta4_eta4]
    type = KokkosMaterialRealAux
    variable = d2_eta4_eta4
    property = 'd^2f_loc/deta4^2'
    execute_on = 'initial timestep_end'
  []
  [d2_eta1_eta2]
    type = KokkosMaterialRealAux
    variable = d2_eta1_eta2
    property = 'd^2f_loc/deta1deta2'
    execute_on = 'initial timestep_end'
  []
  [d2_eta1_eta3]
    type = KokkosMaterialRealAux
    variable = d2_eta1_eta3
    property = 'd^2f_loc/deta1deta3'
    execute_on = 'initial timestep_end'
  []
  [d2_eta1_eta4]
    type = KokkosMaterialRealAux
    variable = d2_eta1_eta4
    property = 'd^2f_loc/deta1deta4'
    execute_on = 'initial timestep_end'
  []
  [d2_eta2_eta3]
    type = KokkosMaterialRealAux
    variable = d2_eta2_eta3
    property = 'd^2f_loc/deta2deta3'
    execute_on = 'initial timestep_end'
  []
  [d2_eta2_eta4]
    type = KokkosMaterialRealAux
    variable = d2_eta2_eta4
    property = 'd^2f_loc/deta2deta4'
    execute_on = 'initial timestep_end'
  []
  [d2_eta3_eta4]
    type = KokkosMaterialRealAux
    variable = d2_eta3_eta4
    property = 'd^2f_loc/deta3deta4'
    execute_on = 'initial timestep_end'
  []
[]

[Materials]
  [local_energy]
    type = KokkosTwoPhaseParabolicFreeEnergy
    property_name = f_loc
    c = c
    etas = 'eta1 eta2 eta3 eta4'
  []
[]

[Debug]
  show_material_props = true
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
  csv = true
  execute_on = 'initial timestep_end'
[]
