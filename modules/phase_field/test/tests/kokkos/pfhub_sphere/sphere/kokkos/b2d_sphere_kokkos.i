# Kokkos twin of ../cpu_control/b2d_sphere_cpu.i (gate G6 of port 05): the system of
# ../../flat/kokkos/b2d_flat_kokkos.i (byte for byte) on a copy of the PFHub sphere mesh (6,146 nodes,
# 6,144 QUAD4 in 3D coordinates; SHA-256 4878bdbdf8f52e44315631a80ea7f675c4dbae25bf6eb39589bf1f181044860e),
# dim = 2, no refinement, the benchmark deck's initial conditions verbatim, ConstantDT 0.01 for three
# steps. Needs the Kokkos element map for a mesh of dimension 2 in 3D coordinates (the surface
# metric). It differs from the CPU control by the Kokkos prefix, the free-energy block, the energy
# block, the absent coupled_variables on KokkosACInterface, and the area integral: the device
# integrates a constant material property of value 1 (KokkosGenericConstantMaterial 'one'), which must
# match the host assembly's area (125,598.09093269) within a relative 1e-12. The quads of the sphere are
# not affine, so the f_density field also exercises the JxW-weighted element average of the Kokkos
# auxiliary kernel.

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
  [area_one]
    type = KokkosGenericConstantMaterial
    prop_names = 'one'
    prop_values = '1'
  []
[]

[Postprocessors]
  [free_energy]
    type = KokkosElementIntegralMaterialProperty
    mat_prop = f_density
    execute_on = 'initial timestep_end'
  []
  [area]
    type = KokkosElementIntegralMaterialProperty
    mat_prop = one
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
