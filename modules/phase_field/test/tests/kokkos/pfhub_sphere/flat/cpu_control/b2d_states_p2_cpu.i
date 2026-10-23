# CPU control of gate G2b of port 05: the 21 properties of the benchmark deck's DerivativeParsedMaterial
# (derivative_order = 2) at 24 fixed states, one state per element, written to CONSTANT MONOMIAL
# fields by 21 MaterialRealAux, no solve. The Kokkos twin is ../kokkos/b2d_states_p2_kokkos.i.
# The 24 states are copied verbatim from the second reader's fixed-state check (12 structured
# states and 12 random states in [-0.1, 1.1], seed 20261022);
# a 2D strip of 24 x 1 QUAD4, each state constant over its element (PiecewiseConstant with
# direction = left; the CONSTANT MONOMIAL initial condition is exact for a function constant on an element).

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 24
  ny = 1
  xmin = 0
  xmax = 24
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
    x = '0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23'
    y = '0.2 0.2 0.2 0.2 0.45 0.45 0.45 0.45 0.65 0.65 0.65 0.65 0.715699 0.176537 0.270241 0.080347 0.069118 0.448907 0.013952 0.448237 0.14959 -0.093082 -0.010691 0.611747'
  []
  [eta1_state]
    type = PiecewiseConstant
    axis = x
    direction = left
    x = '0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23'
    y = '0.1 0.3 0.6 0.9 0.1 0.3 0.6 0.9 0.1 0.3 0.6 0.9 0.595144 -0.096489 0.742462 0.038984 0.317245 0.998567 0.113224 0.818362 -0.030568 -0.009825 0.539589 1.089237'
  []
  [eta2_state]
    type = PiecewiseConstant
    axis = x
    direction = left
    x = '0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23'
    y = '0.3 0.6 0.9 0.1 0.3 0.6 0.9 0.1 0.3 0.6 0.9 0.1 0.92466 0.890599 0.689195 0.608288 -0.031539 0.575482 0.638336 -0.045026 -0.01336 0.163962 0.810739 -0.093599'
  []
  [eta3_state]
    type = PiecewiseConstant
    axis = x
    direction = left
    x = '0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23'
    y = '0.6 0.9 0.1 0.3 0.6 0.9 0.1 0.3 0.6 0.9 0.1 0.3 -0.041781 0.058004 0.027919 -0.093887 0.900465 0.515643 0.399688 0.924369 0.015082 1.092669 0.338763 0.476322'
  []
  [eta4_state]
    type = PiecewiseConstant
    axis = x
    direction = left
    x = '0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23'
    y = '0.9 0.1 0.3 0.6 0.9 0.1 0.3 0.6 0.9 0.1 0.3 0.6 0.888332 0.138622 0.132122 0.037899 0.078004 0.516377 0.006684 1.008706 0.505625 0.475424 0.825518 0.185497'
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
    type = MaterialRealAux
    variable = f
    property = 'f_loc'
    execute_on = 'initial timestep_end'
  []
  [d_c]
    type = MaterialRealAux
    variable = d_c
    property = 'df_loc/dc'
    execute_on = 'initial timestep_end'
  []
  [d2_c_c]
    type = MaterialRealAux
    variable = d2_c_c
    property = 'd^2f_loc/dc^2'
    execute_on = 'initial timestep_end'
  []
  [d_eta1]
    type = MaterialRealAux
    variable = d_eta1
    property = 'df_loc/deta1'
    execute_on = 'initial timestep_end'
  []
  [d2_c_eta1]
    type = MaterialRealAux
    variable = d2_c_eta1
    property = 'd^2f_loc/dcdeta1'
    execute_on = 'initial timestep_end'
  []
  [d2_eta1_eta1]
    type = MaterialRealAux
    variable = d2_eta1_eta1
    property = 'd^2f_loc/deta1^2'
    execute_on = 'initial timestep_end'
  []
  [d_eta2]
    type = MaterialRealAux
    variable = d_eta2
    property = 'df_loc/deta2'
    execute_on = 'initial timestep_end'
  []
  [d2_c_eta2]
    type = MaterialRealAux
    variable = d2_c_eta2
    property = 'd^2f_loc/dcdeta2'
    execute_on = 'initial timestep_end'
  []
  [d2_eta2_eta2]
    type = MaterialRealAux
    variable = d2_eta2_eta2
    property = 'd^2f_loc/deta2^2'
    execute_on = 'initial timestep_end'
  []
  [d_eta3]
    type = MaterialRealAux
    variable = d_eta3
    property = 'df_loc/deta3'
    execute_on = 'initial timestep_end'
  []
  [d2_c_eta3]
    type = MaterialRealAux
    variable = d2_c_eta3
    property = 'd^2f_loc/dcdeta3'
    execute_on = 'initial timestep_end'
  []
  [d2_eta3_eta3]
    type = MaterialRealAux
    variable = d2_eta3_eta3
    property = 'd^2f_loc/deta3^2'
    execute_on = 'initial timestep_end'
  []
  [d_eta4]
    type = MaterialRealAux
    variable = d_eta4
    property = 'df_loc/deta4'
    execute_on = 'initial timestep_end'
  []
  [d2_c_eta4]
    type = MaterialRealAux
    variable = d2_c_eta4
    property = 'd^2f_loc/dcdeta4'
    execute_on = 'initial timestep_end'
  []
  [d2_eta4_eta4]
    type = MaterialRealAux
    variable = d2_eta4_eta4
    property = 'd^2f_loc/deta4^2'
    execute_on = 'initial timestep_end'
  []
  [d2_eta1_eta2]
    type = MaterialRealAux
    variable = d2_eta1_eta2
    property = 'd^2f_loc/deta1deta2'
    execute_on = 'initial timestep_end'
  []
  [d2_eta1_eta3]
    type = MaterialRealAux
    variable = d2_eta1_eta3
    property = 'd^2f_loc/deta1deta3'
    execute_on = 'initial timestep_end'
  []
  [d2_eta1_eta4]
    type = MaterialRealAux
    variable = d2_eta1_eta4
    property = 'd^2f_loc/deta1deta4'
    execute_on = 'initial timestep_end'
  []
  [d2_eta2_eta3]
    type = MaterialRealAux
    variable = d2_eta2_eta3
    property = 'd^2f_loc/deta2deta3'
    execute_on = 'initial timestep_end'
  []
  [d2_eta2_eta4]
    type = MaterialRealAux
    variable = d2_eta2_eta4
    property = 'd^2f_loc/deta2deta4'
    execute_on = 'initial timestep_end'
  []
  [d2_eta3_eta4]
    type = MaterialRealAux
    variable = d2_eta3_eta4
    property = 'd^2f_loc/deta3deta4'
    execute_on = 'initial timestep_end'
  []
[]

[Materials]
  [local_energy]
    type = DerivativeParsedMaterial
    property_name = f_loc
    coupled_variables = 'c eta1 eta2 eta3 eta4'
    constant_names = 'rho_sq ca cb w alpha'
    constant_expressions = '2 0.3 0.7 1 5'
    derivative_order = 2
    expression = 'h:= eta1^3*(6*eta1^2 - 15*eta1 + 10) + eta2^3*(6*eta2^2 - 15*eta2 + 10) + eta3^3*(6*eta3^2 - 15*eta3 + 10) + eta4^3*(6*eta4^2 - 15*eta4 + 10);
                  g:= eta1^2*(1 - eta1)^2 + eta2^2*(1 - eta2)^2 + eta3^2*(1 - eta3)^2 + eta4^2*(1 - eta4)^2
                        + 2*alpha*(eta1^2*eta2^2 + eta1^2*eta3^2 + eta1^2*eta4^2 + eta2^2*eta3^2 + eta2^2*eta4^2 + eta3^2*eta4^2);
                  rho_sq*(c - ca)^2*(1 - h) + rho_sq*(cb - c)^2*h + w*g'
  []
[]

[Debug]
  show_material_props = true
[]

[Executioner]
  type = Steady
[]

[Outputs]
  csv = true
  execute_on = 'initial timestep_end'
  [exodus]
    type = Exodus
    execute_input_on = none
  []
[]
