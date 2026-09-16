# Port 02 tier 5, dispatch gate for the KokkosMatCoupledForce configuration rule (Codex H04). The framework
# kernel creates its array of material-property handles only when material_properties is given
# (KokkosMatCoupledForce.K:57-65); an array of deep-copy element types that was never created is rejected
# by the dispatcher's copy of the kernel (KokkosArray.h:1381-1388, KokkosDispatcher.h:70-78). As written,
# with no material_properties on u_v and v_u, this deck must fail with
#   "Kokkos array error: cannot deep copy using constructor from array without host data"
# and with
#   Kernels/u_v/material_properties=one Kernels/v_u/material_properties=one
# (the E10 rule: a constant property of value 1 from KokkosGenericConstantMaterial, coef carrying the
# coefficient) it must run and pass the Jacobian gate. The two couplings mirror E10's g_s (coef = -18.3,
# residual +S M c) and c_u (coef = 1, residual -M u). Function initial conditions, natural boundaries,
# NEWTON with SMP full = true, two steps of dt = 0.05.

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

[Materials]
  [one]
    type = KokkosGenericConstantMaterial
    prop_names = 'one'
    prop_values = '1'
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
  [u_v]
    type = KokkosMatCoupledForce
    variable = u
    v = v
    coef = -18.3
  []
  [v_dt]
    type = KokkosTimeDerivative
    variable = v
  []
  [v_diff]
    type = KokkosDiffusion
    variable = v
  []
  [v_u]
    type = KokkosMatCoupledForce
    variable = v
    v = u
    coef = 1
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
  num_steps = 2
  dt = 0.05
[]

[Outputs]
  exodus = true
[]
