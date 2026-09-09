# Explicit CPU control of port/decks/test_fixed.i (P2 fixed-mesh circular-grain gate,
# derived from modules/phase_field/test/tests/grain_growth/test.i @ 6956830688).
# Changes vs the P2 deck (DESIGN-kokkos-grain-growth.md section 8.5.2): [Kernels/PolycrystalKernel] expanded
# into explicit kernels; the material is GBEvolution; hypre replaced by gamg with the assembled-matrix
# smoother (D-021). Mesh, ICs, aux, postprocessor, SMP, NEWTON, tolerances, dt and num_steps are unchanged.
# Non-periodic (natural boundaries), as the source deck.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  nz = 0
  xmin = 0
  xmax = 400
  ymin = 0
  ymax = 400
  zmin = 0
  zmax = 0
  elem_type = QUAD4
  uniform_refine = 2
[]

[GlobalParams]
  op_num = 2
  var_name_base = gr
[]

[Variables]
  [PolycrystalVariables]
  []
[]

[ICs]
  [PolycrystalICs]
    [BicrystalCircleGrainIC]
      radius = 300
      x = 400
      y = 0
    []
  []
[]

[AuxVariables]
  [bnds]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  [dt0]
    type = TimeDerivative
    variable = gr0
  []
  [bulk0]
    type = ACGrGrPoly
    variable = gr0
    v = 'gr1'
    mob_name = L
  []
  [int0]
    type = ACInterface
    variable = gr0
    mob_name = L
    kappa_name = kappa_op
    variable_L = true
  []
  [dt1]
    type = TimeDerivative
    variable = gr1
  []
  [bulk1]
    type = ACGrGrPoly
    variable = gr1
    v = 'gr0'
    mob_name = L
  []
  [int1]
    type = ACInterface
    variable = gr1
    mob_name = L
    kappa_name = kappa_op
    variable_L = true
  []
[]

[AuxKernels]
  [BndsCalc]
    type = BndsCalcAux
    variable = bnds
  []
[]

[Materials]
  [Copper]
    type = GBEvolution
    T = 500 # K
    wGB = 60 # nm
    GBmob0 = 2.5e-6 #m^4/(Js) from Schoenfelder 1997
    Q = 0.23 #Migration energy in eV
    GBenergy = 0.708 #GB energy in J/m^2
  []
[]

[Postprocessors]
  [gr1area]
    type = ElementIntegralVariablePostprocessor
    variable = gr1
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = 'NEWTON'

  petsc_options_iname = '-pc_type -pc_use_amat -ksp_gmres_restart'
  petsc_options_value = 'gamg false 31'

  l_tol = 1.0e-4
  l_max_its = 30
  nl_max_its = 20
  nl_rel_tol = 1.0e-9
  start_time = 0.0
  num_steps = 5
  dt = 80.0
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
  csv = true
[]
