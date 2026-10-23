# Kokkos twin of ../cpu_control/b2d_flat_cpu.i (gates G3 and G5 of port 05): the PFHub Benchmark 2
# part (d) physics on a flat 20 x 20 QUAD4 mesh, every kernel and material on the device. It differs
# from the CPU control by the Kokkos prefix and three named blocks only:
#   - the free energy: KokkosTwoPhaseParabolicFreeEnergy (the closed form of the deck's
#     DerivativeParsedMaterial, the same 21 property names, its defaults the deck's constants);
#   - the energy: KokkosTotalFreeEnergyDensity declares f_density (TotalFreeEnergy's formula),
#     KokkosMaterialRealAux writes it to the f_density field and KokkosElementIntegralMaterialProperty
#     integrates it as free_energy (a host aux kernel or postprocessor over a Kokkos property reads a
#     silent zero);
#   - KokkosACInterface carries no coupled_variables (it accepts only an auxiliary T; every
#     off-diagonal term of ACInterface is zero at a constant L and kappa).
# The Exodus files of both twins carry the same names: the six nodal fields, f_density, free_energy.

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
[]

[Postprocessors]
  [free_energy]
    type = KokkosElementIntegralMaterialProperty
    mat_prop = f_density
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
