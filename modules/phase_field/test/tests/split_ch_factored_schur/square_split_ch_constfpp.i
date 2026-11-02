# The constant-f'' variant of square_split_ch.i: the double well replaced by the quadratic
# f = 5 (c - 0.5)^2, so f'' = 10 everywhere and the f'' block is 10 times the mass; the fit
# A01 = beta_fit A10 + gamma A00 is then exact (beta_fit = kappa/M = 0.4, gamma = -10). Nothing else
# differs from square_split_ch.i.

[Mesh]
  [square]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 8
    ny = 8
    xmax = 16
    ymax = 16
    elem_type = QUAD4
  []
[]

[Variables]
  [c]
  []
  [chem_pot_c]
  []
[]

[Functions]
  [ic_func]
    type = ParsedFunction
    expression = '0.5 + 0.15 * cos(2 * pi * x / 16) * cos(2 * pi * y / 16)'
  []
[]

[ICs]
  [c]
    type = FunctionIC
    variable = c
    function = ic_func
  []
[]

[Kernels]
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
  []
  [chem_pot_c_CoefReaction]
    type = CoefReaction
    variable = chem_pot_c
    coefficient = -1
  []
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'M kappa_c'
    prop_values = '5 2'
  []

  [local_energy]
    type = DerivativeParsedMaterial
    property_name = f_loc
    coupled_variables = c
    constant_names = 'rho_s cm'
    constant_expressions = '5 0.5'
    expression = 'rho_s*(c-cm)^2'
  []
[]

[Postprocessors]
  [average_c]
    type = ElementAverageValue
    variable = c
    execute_on = 'initial timestep_end'
  []
  [min_c_nodal]
    type = NodalExtremeValue
    variable = c
    value_type = min
    execute_on = 'initial timestep_end'
  []
  [max_c_nodal]
    type = NodalExtremeValue
    variable = c
    value_type = max
    execute_on = 'initial timestep_end'
  []
  [dt]
    type = TimestepSize
  []
  [nl_its]
    type = NumNonlinearIterations
  []
  [l_its]
    type = NumLinearIterations
  []
[]

[Preconditioning]
  active = 'fs_factored'
  [fs_factored]
    type = SplitCHFactoredSchurPreconditioner
    full = true # keep: full coupling puts every block pair in the pattern
    topsplit = wc
    mobility = 5
    kappa = 2
    schur_form = consistent_mass
    split_parameter = 1
    [wc]
      splitting = 'w c'
      splitting_type = schur
      schur_type = full
      schur_pre = A11
    []
    [w]
      vars = 'chem_pot_c'
      petsc_options_iname = '-ksp_type -pc_type -pc_jacobi_fixdiagonal'
      petsc_options_value = 'preonly jacobi false'
    []
    [c]
      vars = 'c'
      petsc_options_iname = '-ksp_type -pc_type
          -stiffness_factor_ksp_type -stiffness_factor_pc_type
          -hessian_factor_ksp_type -hessian_factor_pc_type'
      petsc_options_value = 'preonly shell
          preonly lu
          preonly lu'
    []
  []
  [fs_lu]
    type = FSP
    full = true
    topsplit = 'wc'
    [wc]
      splitting = 'w c'
      splitting_type = schur
      schur_type = full
      schur_pre = Sp
    []
    [w]
      vars = 'chem_pot_c'
      petsc_options_iname = '-ksp_type -pc_type -pc_jacobi_fixdiagonal'
      petsc_options_value = 'preonly jacobi false'
    []
    [c]
      vars = 'c'
      petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_type'
      petsc_options_value = 'preonly lu mumps'
    []
  []
  [smp_full]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  l_tol = 1e-4
  l_max_its = 1000
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-9
  nl_max_its = 10
  automatic_scaling = false
  start_time = 0
  end_time = 3
  [TimeStepper]
    type = TimeSequenceStepper
    time_sequence = '1 3'
  []
[]

[Outputs]
  csv = true
[]
