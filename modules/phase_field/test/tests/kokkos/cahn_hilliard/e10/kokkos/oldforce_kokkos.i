# Port 02 tier 5, toy gate for KokkosOldCoupledForce (the E10 u_cold term): two diffusing fields u and v.
# u carries the source -coef * v_old, v at the previous time step, the Kokkos form of CoupledForce with
# implicit = false. CPU control: ../cpu_control/oldforce_cpu.i, identical apart from the Kokkos prefix and
# the implicit = false line that the Kokkos object replaces. Function initial conditions inside (0, 1),
# natural boundaries, NEWTON with SMP full = true, three steps of dt = 0.05.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 8
  ny = 8
  xmax = 1
  ymax = 1
  elem_type = QUAD4
[]

[Variables]
  [u]
  []
  [v]
  []
[]

[Functions]
  [u_ic]
    type = ParsedFunction
    expression = '0.5 + 0.2*cos(pi*x)*cos(pi*y)'
  []
  [v_ic]
    type = ParsedFunction
    expression = '0.4 + 0.3*cos(2*pi*x)*cos(pi*y)'
  []
[]

[ICs]
  [u]
    type = FunctionIC
    variable = u
    function = u_ic
  []
  [v]
    type = FunctionIC
    variable = v
    function = v_ic
  []
[]

[Kernels]
  [u_dt]
    type = KokkosTimeDerivative
    variable = u
  []
  [u_diff]
    type = KokkosDiffusion
    variable = u
  []
  [u_vold]
    type = KokkosOldCoupledForce
    variable = u
    v = v
    coef = 2.0
  []
  [v_dt]
    type = KokkosTimeDerivative
    variable = v
  []
  [v_diff]
    type = KokkosDiffusion
    variable = v
  []
  [v_react]
    type = KokkosReaction
    variable = v
    rate = 0.5
  []
[]

[Postprocessors]
  [u_int]
    type = KokkosElementIntegralVariablePostprocessor
    variable = u
  []
  [v_int]
    type = KokkosElementIntegralVariablePostprocessor
    variable = v
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
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -sub_pc_type'
  petsc_options_value = 'bjacobi lu'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  l_tol = 1e-12
  num_steps = 3
  dt = 0.05
[]

[Outputs]
  exodus = true
  csv = true
[]
