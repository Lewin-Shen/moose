# CPU control of the sphere twin of port 05 (gates G0 and G6): the PFHub Benchmark 2 part (d) system of
# ../../flat/cpu_control/b2d_flat_cpu.i (the same kernels, materials, energy block, solver and
# tolerances, byte for byte) on the benchmark's sphere mesh, a copy of the PFHub sphere_mesh.e
# (6,146 nodes, 6,144 QUAD4 in 3D coordinates; SHA-256
# 4878bdbdf8f52e44315631a80ea7f675c4dbae25bf6eb39589bf1f181044860e), with dim = 2 and no refinement
# (the benchmark deck refines twice), the benchmark deck's five ParsedFunction initial conditions
# verbatim, ConstantDT 0.01 for three steps. The surface area is the integral of a constant aux
# variable of value 1 (not written to the Exodus file, so both twins carry the same names); the host
# assembly gives 125,598.09093269 on this mesh. The Kokkos twin is ../kokkos/b2d_sphere_kokkos.i.

[Mesh]
  [file]
    type = FileMeshGenerator
    file = ../sphere_mesh.e
  []
  dim = 2
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
  [ic_func_c]
    type = ParsedFunction
    symbol_names = 'c0 eps'
    symbol_values = '0.5 0.05'
    expression = 'c0 + eps*(cos(8*acos(z/sqrt(x^2 + y^2 + z^2)))*cos(15*atan2(y,x)) + (cos(12*acos(z/sqrt(x^2 + y^2 + z^2)))*cos(10*atan2(y,x)))^2 + cos(2.5*acos(z/sqrt(x^2 + y^2 + z^2)) - 1.5*atan2(y,x))*cos(7*acos(z/sqrt(x^2 + y^2 + z^2)) - 2*atan2(y,x)))'
  []
  [ic_func_eta1]
    type = ParsedFunction
    symbol_names = 'eps_eta psi'
    symbol_values = '0.1 1.5'
    expression = 'eps_eta*(cos(1*acos(z/sqrt(x^2 + y^2 + z^2)) - 4)*cos(1.7*atan2(y,x)) + cos(12*acos(z/sqrt(x^2 + y^2 + z^2)))*cos(12*atan2(y,x)) + psi*(cos(4.7*acos(z/sqrt(x^2 + y^2 + z^2)) + 4.15*atan2(y,x))*cos(3.2*acos(z/sqrt(x^2 + y^2 + z^2)) - 0.5*atan2(y,x)))^2)^2'
  []

  [ic_func_eta2]
    type = ParsedFunction
    symbol_names = 'eps_eta psi'
    symbol_values = '0.1 1.5'
    expression = 'eps_eta*(cos(2*acos(z/sqrt(x^2 + y^2 + z^2)) - 4)*cos(2.7*atan2(y,x)) + cos(13*acos(z/sqrt(x^2 + y^2 + z^2)))*cos(13*atan2(y,x)) + psi*(cos(4.8*acos(z/sqrt(x^2 + y^2 + z^2)) + 4.25*atan2(y,x))*cos(3.3*acos(z/sqrt(x^2 + y^2 + z^2)) - 0.6*atan2(y,x)))^2)^2'
  []

  [ic_func_eta3]
    type = ParsedFunction
    symbol_names = 'eps_eta psi'
    symbol_values = '0.1 1.5'
    expression = 'eps_eta*(cos(3*acos(z/sqrt(x^2 + y^2 + z^2)) - 4)*cos(3.7*atan2(y,x)) + cos(14*acos(z/sqrt(x^2 + y^2 + z^2)))*cos(14*atan2(y,x)) + psi*(cos(4.9*acos(z/sqrt(x^2 + y^2 + z^2)) + 4.35*atan2(y,x))*cos(3.4*acos(z/sqrt(x^2 + y^2 + z^2)) - 0.7*atan2(y,x)))^2)^2'
  []

  [ic_func_eta4]
    type = ParsedFunction
    symbol_names = 'eps_eta psi'
    symbol_values = '0.1 1.5'
    expression = 'eps_eta*(cos(4*acos(z/sqrt(x^2 + y^2 + z^2)) - 4)*cos(4.7*atan2(y,x)) + cos(15*acos(z/sqrt(x^2 + y^2 + z^2)))*cos(15*atan2(y,x)) + psi*(cos(5.0*acos(z/sqrt(x^2 + y^2 + z^2)) + 4.45*atan2(y,x))*cos(3.5*acos(z/sqrt(x^2 + y^2 + z^2)) - 0.8*atan2(y,x)))^2)^2'
  []
[]

[ICs]
  [concentrationIC]
    type = FunctionIC
    variable = c
    function = ic_func_c
  []
  [eta1_ic]
    type = FunctionIC
    variable = eta1
    function = ic_func_eta1
  []
  [eta2_ic]
    type = FunctionIC
    variable = eta2
    function = ic_func_eta2
  []
  [eta3_ic]
    type = FunctionIC
    variable = eta3
    function = ic_func_eta3
  []
  [eta4_ic]
    type = FunctionIC
    variable = eta4
    function = ic_func_eta4
  []
[]

[Kernels]
  # ConservedAction, FORWARD_SPLIT, on c and chem_pot_c
  [c_TimeDerivative]
    type = TimeDerivative
    variable = c
  []
  [c_MatDiffusion]
    type = MatDiffusion
    variable = c
    v = chem_pot_c
    diffusivity = M
  []
  [chem_pot_c_MatDiffusion]
    type = MatDiffusion
    variable = chem_pot_c
    v = c
    diffusivity = kappa_c
  []
  [chem_pot_c_CoupledMaterialDerivative]
    type = CoupledMaterialDerivative
    variable = chem_pot_c
    v = c
    f_name = f_loc
    coupled_variables = 'eta1 eta2 eta3 eta4'
  []
  [chem_pot_c_CoefReaction]
    type = Reaction
    variable = chem_pot_c
    rate = -1
  []

  # NonconservedAction on eta1 to eta4
  [eta1_TimeDerivative]
    type = TimeDerivative
    variable = eta1
  []
  [eta1_AllenCahn]
    type = AllenCahn
    variable = eta1
    mob_name = L
    f_name = f_loc
    coupled_variables = 'c eta2 eta3 eta4'
  []
  [eta1_ACInterface]
    type = ACInterface
    variable = eta1
    mob_name = L
    kappa_name = kappa_eta
    variable_L = true
    coupled_variables = 'c eta2 eta3 eta4'
  []

  [eta2_TimeDerivative]
    type = TimeDerivative
    variable = eta2
  []
  [eta2_AllenCahn]
    type = AllenCahn
    variable = eta2
    mob_name = L
    f_name = f_loc
    coupled_variables = 'c eta1 eta3 eta4'
  []
  [eta2_ACInterface]
    type = ACInterface
    variable = eta2
    mob_name = L
    kappa_name = kappa_eta
    variable_L = true
    coupled_variables = 'c eta1 eta3 eta4'
  []

  [eta3_TimeDerivative]
    type = TimeDerivative
    variable = eta3
  []
  [eta3_AllenCahn]
    type = AllenCahn
    variable = eta3
    mob_name = L
    f_name = f_loc
    coupled_variables = 'c eta1 eta2 eta4'
  []
  [eta3_ACInterface]
    type = ACInterface
    variable = eta3
    mob_name = L
    kappa_name = kappa_eta
    variable_L = true
    coupled_variables = 'c eta1 eta2 eta4'
  []

  [eta4_TimeDerivative]
    type = TimeDerivative
    variable = eta4
  []
  [eta4_AllenCahn]
    type = AllenCahn
    variable = eta4
    mob_name = L
    f_name = f_loc
    coupled_variables = 'c eta1 eta2 eta3'
  []
  [eta4_ACInterface]
    type = ACInterface
    variable = eta4
    mob_name = L
    kappa_name = kappa_eta
    variable_L = true
    coupled_variables = 'c eta1 eta2 eta3'
  []
[]

[AuxVariables]
  [f_density]
    order = CONSTANT
    family = MONOMIAL
  []
  [one]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1
    outputs = none
  []
[]

[AuxKernels]
  [f_density]
    type = TotalFreeEnergy
    variable = f_density
    f_name = 'f_loc'
    kappa_names = 'kappa_c kappa_eta kappa_eta kappa_eta kappa_eta'
    interfacial_vars = 'c eta1 eta2 eta3 eta4'
    execute_on = 'initial timestep_end'
  []
[]

[Materials]
  [constants_c]
    type = GenericConstantMaterial
    prop_names = 'M kappa_c'
    prop_values = '5 3'
  []
  [constants_eta]
    type = GenericConstantMaterial
    prop_names = 'L kappa_eta'
    prop_values = '5 3'
  []
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

[Postprocessors]
  [free_energy]
    type = ElementIntegralVariablePostprocessor
    variable = f_density
    execute_on = 'initial timestep_end'
  []
  [area]
    type = ElementIntegralVariablePostprocessor
    variable = one
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
  num_steps = 3
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
    execute_on = 'initial final'
    execute_input_on = none
  []
[]
