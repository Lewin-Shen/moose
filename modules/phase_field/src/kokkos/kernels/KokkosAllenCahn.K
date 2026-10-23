//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosAllenCahn.h"

#include "MooseVariable.h"
#include "SystemBase.h"

#include <algorithm>

registerKokkosResidualObject("PhaseFieldApp", KokkosAllenCahn);

InputParameters
KokkosAllenCahn::validParams()
{
  InputParameters params = Moose::Kokkos::KernelValue::validParams();
  params.addClassDescription(
      "Allen-Cahn Kernel that uses a DerivativeMaterial Free Energy, for a constant mobility (the "
      "mobility derivative terms are omitted; a declared mobility derivative is an error)");
  params.addCoupledVar("coupled_variables",
                       "Vector of nonlinear variable arguments this object depends on");
  params.addParam<MaterialPropertyName>("mob_name", "L", "The mobility used with the kernel");
  params.addRequiredParam<MaterialPropertyName>(
      "f_name", "Base name of the free energy function F defined in a DerivativeParsedMaterial");
  return params;
}

KokkosAllenCahn::KokkosAllenCahn(const InputParameters & parameters)
  : Moose::Kokkos::KernelValue(parameters),
    _n_args(coupledComponents("coupled_variables")),
    _L(getKokkosMaterialProperty<Real>("mob_name"))
{
  const auto F = getMaterialPropertyName("f_name");
  const auto mobility = getMaterialPropertyName("mob_name");
  const auto uname = _var.name();

  // The residual property is required: the CPU kernel would read a zero default and contribute
  // nothing, which is never what a deck means
  const auto dFdu_name = derivativePropertyNameFirst(F, uname);
  if (!hasKokkosMaterialPropertyByName<Real>(dFdu_name))
    mooseError("KokkosAllenCahn: no Kokkos material declares '",
               dFdu_name,
               "'; check f_name and the kernel variable (the CPU kernel would silently use a zero "
               "property here).");
  _dFdu = getKokkosMaterialPropertyByName<Real>(dFdu_name);

  // Optional second derivatives under the CPU names: fetched only when declared, zero otherwise
  const auto fetch = [this](const MaterialPropertyName & derivative,
                            Moose::Kokkos::MaterialProperty<Real> & handle)
  {
    if (hasKokkosMaterialPropertyByName<Real>(derivative))
      handle = getKokkosMaterialPropertyByName<Real>(derivative);
  };
  fetch(derivativePropertyNameSecond(F, uname, uname), _d2Fdu2);

  // Constant mobility: ACBulk's dL/du and dL/darg product-rule terms are omitted by construction, so
  // a material that declares either derivative is outside this specialization
  const auto reject_mobility_derivative = [&](const VariableName & x)
  {
    const auto derivative = derivativePropertyNameFirst(mobility, x);
    if (hasKokkosMaterialPropertyByName<Real>(derivative))
      paramError("mob_name",
                 "KokkosAllenCahn supports a constant mobility only, but a material declares '",
                 derivative,
                 "'.");
  };
  reject_mobility_derivative(uname);

  // Always created, possibly empty: the dispatcher deep-copies the property array and aborts on
  // one that was never allocated (KokkosArray.h deepCopy)
  std::vector<VariableName> coupled_names;
  _arg_var.create(_n_args);
  _d2Fdudarg.create(_n_args);
  for (unsigned int k = 0; k < _n_args; ++k)
  {
    if (!isCoupled("coupled_variables", k))
      paramError("coupled_variables", "Entries must be variables, not literal values.");
    const auto * arg = getVar("coupled_variables", k);
    if (arg->name() == uname)
      paramError("coupled_variables", "The kernel variable must not appear in coupled_variables.");
    // Only a nonlinear variable of this kernel's system has a Jacobian column
    _arg_var[k] = (arg->kind() == Moose::VAR_SOLVER && &arg->sys() == &_var.sys())
                      ? arg->number()
                      : libMesh::invalid_uint;
    fetch(derivativePropertyNameSecond(F, uname, arg->name()), _d2Fdudarg[k]);
    reject_mobility_derivative(arg->name());
    coupled_names.push_back(arg->name());
  }
  if (_n_args)
    _arg_var.copyToDevice();

  // DerivativeMaterialInterface::validateNonlinearCoupling for "mob_name" (ACBulk::initialSetup)
  // and "f_name" (AllenCahn::initialSetup), issued at construction: a declared dF/dx or dL/dx for a
  // nonlinear x that is neither the kernel variable nor coupled here means a Jacobian block would
  // be silently dropped
  for (const auto & base : std::vector<MaterialPropertyName>{mobility, F})
  {
    std::vector<VariableName> missing;
    for (const auto & x : _var.sys().getVariableNames())
    {
      if (x == uname)
        continue;
      if (std::find(coupled_names.begin(), coupled_names.end(), x) != coupled_names.end())
        continue;
      if (hasKokkosMaterialPropertyByName<Real>(derivativePropertyNameFirst(base, x)))
        missing.push_back(x);
    }
    if (!missing.empty())
      mooseWarning("Missing coupled variables {",
                   Moose::stringify(missing),
                   "} (add them to coupled_variables parameter of ",
                   name(),
                   ")");
  }

  // A misspelled hook is silently replaced by the base default (KokkosDispatcher.h); fail loudly
  if (!Moose::Kokkos::hasUserJacobianHook<KokkosAllenCahn>() ||
      !Moose::Kokkos::hasUserOffDiagJacobianHook<KokkosAllenCahn>())
    mooseError("KokkosAllenCahn Jacobian hooks were not registered with the expected signatures.");
}
