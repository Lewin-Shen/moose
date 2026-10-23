# CPU control of the Jacobian gate (G2) of port 05: the six-field PFHub Benchmark 2 part (d) system of
# b2d_flat_cpu.i on a 4 x 4 QUAD4 mesh of the unit square, one step, at a state where every block of
# the block Jacobian is nonzero: c = 0.4 + 0.05 cos(pi x) cos(pi y) in [0.35, 0.45] (off c = 0.5,
# where D = B - A = 0, and off ca = 0.3), eta_k = b_k + 0.05 sin(pi (x + k y)) with
# b = (0.2, 0.35, 0.6, 0.85) (distinct, never 0, 0.5 or 1, no pair mirrored about 0.5) and
# chem_pot_c = 0.1 sin(pi x) sin(pi y). The same finite-difference check on the CPU objects gives the
# tester's floor for the Kokkos twin ../kokkos/b2d_jac_kokkos.i, which differs by the Kokkos prefix,
# the free-energy block, the energy block and the absent coupled_variables on the interface kernel.
# Solver settings as b2d_flat_cpu.i (NEWTON, SMP full = true, asm + lu, automatic_scaling = false).

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
    expression = '0.4 + 0.05*cos(pi*x)*cos(pi*y)'
  []
  [chem_pot_c_ic]
    type = ParsedFunction
    expression = '0.1*sin(pi*x)*sin(pi*y)'
  []
  [eta1_ic]
    type = ParsedFunction
    expression = '0.2 + 0.05*sin(pi*(x + 1*y))'
  []
  [eta2_ic]
    type = ParsedFunction
    expression = '0.35 + 0.05*sin(pi*(x + 2*y))'
  []
  [eta3_ic]
    type = ParsedFunction
    expression = '0.6 + 0.05*sin(pi*(x + 3*y))'
  []
  [eta4_ic]
    type = ParsedFunction
    expression = '0.85 + 0.05*sin(pi*(x + 4*y))'
  []
[]

[ICs]
  [c_ic]
    type = FunctionIC
    variable = c
    function = c_ic
  []
  [chem_pot_c_ic]
    type = FunctionIC
    variable = chem_pot_c
    function = chem_pot_c_ic
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
  num_steps = 1
  [TimeStepper]
    type = ConstantDT
    dt = 0.01
  []
[]

[Outputs]
  exodus = true
  csv = true
  execute_on = 'initial timestep_end'
[]
