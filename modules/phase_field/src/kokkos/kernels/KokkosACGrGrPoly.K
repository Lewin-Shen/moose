//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosACGrGrPoly.h"

#include "MooseVariable.h"
#include "SystemBase.h"

#include <set>

registerKokkosResidualObject("PhaseFieldApp", KokkosACGrGrPoly);

InputParameters
KokkosACGrGrPoly::validParams()
{
  InputParameters params = Moose::Kokkos::KernelValue::validParams();
  params.addClassDescription("Grain-Boundary model poly-crystalline interface Allen-Cahn Kernel");
  params.addCoupledVar("coupled_variables",
                       "Vector of nonlinear variable arguments this object depends on");
  params.addParam<MaterialPropertyName>("mob_name", "L", "The mobility used with the kernel");
  params.addRequiredCoupledVar(
      "v", "Array of coupled order parameter names for other order parameters");
  return params;
}

KokkosACGrGrPoly::KokkosACGrGrPoly(const InputParameters & parameters)
  : Moose::Kokkos::KernelValue(parameters),
    _op_num(coupledComponents("v")),
    _v_var(coupledIndices("v")),
    _v(kokkosCoupledValues("v")),
    _L(getKokkosMaterialProperty<Real>("mob_name")),
    _mu(getKokkosMaterialPropertyByName<Real>("mu")),
    _gamma(getKokkosMaterialPropertyByName<Real>("gamma_asymm"))
{
  if (!_op_num)
    paramError("v", "At least one other order parameter is required by this specialization.");

  // Every entry of 'v' must be a nonlinear FE variable of this kernel's system, distinct from
  // the kernel variable and from each other (the device map is keyed by variable number).
  // The duplicate check uses a host-side set: Map::exists() is a device query and reads arrays
  // that only copyToDevice() fills (a null dereference on the host before that call).
  std::vector<VariableName> nonlinear_names = {_var.name()};
  std::set<unsigned int> seen;
  for (unsigned int c = 0; c < _op_num; ++c)
  {
    if (!isCoupled("v", c))
      paramError("v", "Entries must be nonlinear FE variables, not constants.");
    const auto * var = getVar("v", c);
    if (var->kind() != Moose::VAR_SOLVER || &var->sys() != &_var.sys())
      paramError("v", "Order parameters must belong to the kernel's nonlinear system.");
    if (_v_var[c] == _var.number())
      paramError("v", "The kernel variable must not appear in v.");
    if (!seen.insert(_v_var[c]).second)
      paramError("v", "Repeated order parameters are not supported.");
    _v_var_to_index[_v_var[c]] = c;
    nonlinear_names.push_back(var->name());
  }
  _v_var_to_index.copyToDevice();

  // 'coupled_variables' is accepted for action compatibility only: the prescribed auxiliary T
  const unsigned int n_args = coupledComponents("coupled_variables");
  if (n_args > 1)
    paramError("coupled_variables", "Only the prescribed auxiliary T is supported.");
  if (n_args)
  {
    if (!isCoupled("coupled_variables", 0))
      paramError("coupled_variables", "Use the auxiliary field T, or omit this list.");
    const auto * arg = getVar("coupled_variables", 0);
    if (arg->name() != "T" || arg->kind() != Moose::VAR_AUXILIARY)
      paramError("coupled_variables", "Only the prescribed auxiliary T is supported.");
  }

  // D-020 scope: L, mu and gamma_asymm must not depend on any order parameter. A declared
  // derivative property means a material outside this specialization is in use.
  const auto mobility = getMaterialPropertyName("mob_name");
  for (const auto & var : nonlinear_names)
    for (const auto & property : std::vector<MaterialPropertyName>{mobility, "mu", "gamma_asymm"})
    {
      const auto derivative = derivativePropertyNameFirst(property, var);
      if (hasKokkosMaterialProperty<Real>(derivative))
        mooseError("Unsupported order-parameter dependence in KokkosACGrGrPoly: ", derivative);
    }

  // A misspelled hook is silently replaced by the base default (KokkosDispatcher.h); fail loudly
  if (!Moose::Kokkos::hasUserJacobianHook<KokkosACGrGrPoly>() ||
      !Moose::Kokkos::hasUserOffDiagJacobianHook<KokkosACGrGrPoly>())
    mooseError(
        "KokkosACGrGrPoly Jacobian hooks were not registered with the expected signatures.");
}
