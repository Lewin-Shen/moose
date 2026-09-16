//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosCoupledMaterialDerivative.h"

#include "MooseVariable.h"
#include "SystemBase.h"

#include <algorithm>

registerKokkosResidualObject("PhaseFieldApp", KokkosCoupledMaterialDerivative);

InputParameters
KokkosCoupledMaterialDerivative::validParams()
{
  InputParameters params = Moose::Kokkos::KernelValue::validParams();
  params.addClassDescription("Kernel that implements the first derivative of a function material "
                             "property with respect to a coupled variable.");
  params.addRequiredCoupledVar("v", "Variable to take the derivative with respect to");
  params.addParam<MaterialPropertyName>("f_name",
                                        "F",
                                        "Function material to take the derivative of (should "
                                        "provide derivative properties - such as a "
                                        "DerivativeParsedMaterial)");
  params.addCoupledVar("coupled_variables",
                       "Vector of nonlinear variable arguments this object depends on");
  return params;
}

KokkosCoupledMaterialDerivative::KokkosCoupledMaterialDerivative(
    const InputParameters & parameters)
  : Moose::Kokkos::KernelValue(parameters),
    _v_var(coupled("v")),
    _n_args(coupledComponents("coupled_variables"))
{
  if (!isCoupled("v") || coupledComponents("v") != 1)
    paramError("v", "Exactly one variable is required (a literal value is not supported).");

  const auto F = getMaterialPropertyName("f_name");
  const auto vname = coupledName("v", 0);

  // The residual property is required: the CPU kernel would read a zero default and contribute
  // nothing, which is never what a deck means
  const auto dFdv_name = derivativePropertyNameFirst(F, vname);
  if (!hasKokkosMaterialPropertyByName<Real>(dFdv_name))
    mooseError("KokkosCoupledMaterialDerivative: no Kokkos material declares '",
               dFdv_name,
               "'; check f_name and v (the CPU kernel would silently use a zero property here).");
  _dFdv = getKokkosMaterialPropertyByName<Real>(dFdv_name);

  // Optional second derivatives under the CPU names: fetched only when declared, zero otherwise
  const auto fetch = [this](const MaterialPropertyName & derivative,
                            Moose::Kokkos::MaterialProperty<Real> & handle)
  {
    if (hasKokkosMaterialPropertyByName<Real>(derivative))
      handle = getKokkosMaterialPropertyByName<Real>(derivative);
  };
  fetch(derivativePropertyNameSecond(F, vname, _var.name()), _d2Fdvdu);
  fetch(derivativePropertyNameSecond(F, vname, vname), _d2Fdv2);

  // Always created, possibly empty: the dispatcher deep-copies the property array and aborts on
  // one that was never allocated (KokkosArray.h deepCopy)
  std::vector<VariableName> coupled_names = {vname};
  _arg_var.create(_n_args);
  _d2Fdvdarg.create(_n_args);
  for (unsigned int k = 0; k < _n_args; ++k)
  {
    if (!isCoupled("coupled_variables", k))
      paramError("coupled_variables", "Entries must be variables, not literal values.");
    const auto * arg = getVar("coupled_variables", k);
    // Only a nonlinear variable of this kernel's system has a Jacobian column
    _arg_var[k] = (arg->kind() == Moose::VAR_SOLVER && &arg->sys() == &_var.sys())
                      ? arg->number()
                      : libMesh::invalid_uint;
    fetch(derivativePropertyNameSecond(F, vname, arg->name()), _d2Fdvdarg[k]);
    coupled_names.push_back(arg->name());
  }
  if (_n_args)
    _arg_var.copyToDevice();

  // DerivativeMaterialInterface::validateNonlinearCoupling("f_name"), issued at construction: a
  // declared dF/dx for a nonlinear x that is neither the kernel variable nor coupled here means a
  // d^2F/dv dx block would be silently dropped from the Jacobian
  std::vector<VariableName> missing;
  for (const auto & x : _var.sys().getVariableNames())
  {
    if (x == _var.name())
      continue;
    if (std::find(coupled_names.begin(), coupled_names.end(), x) != coupled_names.end())
      continue;
    if (hasKokkosMaterialPropertyByName<Real>(derivativePropertyNameFirst(F, x)))
      missing.push_back(x);
  }
  if (!missing.empty())
    mooseWarning("Missing coupled variables {",
                 Moose::stringify(missing),
                 "} (add them to coupled_variables parameter of ",
                 name(),
                 ")");

  // A misspelled hook is silently replaced by the base default (KokkosDispatcher.h); fail loudly
  if (!Moose::Kokkos::hasUserJacobianHook<KokkosCoupledMaterialDerivative>() ||
      !Moose::Kokkos::hasUserOffDiagJacobianHook<KokkosCoupledMaterialDerivative>())
    mooseError("KokkosCoupledMaterialDerivative Jacobian hooks were not registered with the "
               "expected signatures.");
}
