//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosSplitCHParsed.h"

#include "MooseVariable.h"
#include "SystemBase.h"

#include <algorithm>
#include <set>

registerKokkosResidualObject("PhaseFieldApp", KokkosSplitCHParsed);

InputParameters
KokkosSplitCHParsed::validParams()
{
  InputParameters params = Moose::Kokkos::Kernel::validParams();
  params.addClassDescription(
      "Split formulation Cahn-Hilliard Kernel that uses a DerivativeMaterial Free Energy");
  params.addRequiredCoupledVar("w", "Chemical potential");
  params.addRequiredParam<MaterialPropertyName>("kappa_name", "The kappa used with the kernel");
  params.addRequiredParam<MaterialPropertyName>(
      "f_name", "Base name of the free energy function F defined in a DerivativeParsedMaterial");
  params.addCoupledVar("coupled_variables", "Vector of additional variable arguments to F");
  return params;
}

KokkosSplitCHParsed::KokkosSplitCHParsed(const InputParameters & parameters)
  : Moose::Kokkos::Kernel(parameters),
    _kappa(getKokkosMaterialProperty<Real>("kappa_name")),
    _w_var(coupled("w")),
    _w(kokkosCoupledValue("w")),
    _dFdc(getKokkosMaterialPropertyByName<Real>(
        derivativePropertyNameFirst(getMaterialPropertyName("f_name"), _var.name()))),
    _d2Fdc2(getKokkosMaterialPropertyByName<Real>(derivativePropertyNameSecond(
        getMaterialPropertyName("f_name"), _var.name(), _var.name()))),
    _n_args(coupledComponents("coupled_variables"))
{
  // A coupled w is restricted to a nonlinear variable of this kernel's system. This narrows the
  // host, which also runs an auxiliary or literal w (Coupleable::coupled() returns a non-aliasing
  // sentinel for either, Coupleable.C:112-116, 486-496); the aliasing hazard of trap T9 is the
  // host's JvarMap (JvarMapInterface.h:132-138), which the device argument map below avoids by
  // holding nonlinear variables of this system only. Kept as a deliberate restriction: no gate
  // deck and no ConservedAction deck couples w to anything else.
  // A literal default value has no variable object behind it (getVar() returns nullptr for it,
  // Coupleable.h:2181-2186), so it is rejected before the pointer is used
  if (!isCoupled("w"))
    paramError("w", "The chemical potential must be a variable, not a constant.");
  const auto * w = getVar("w", 0);
  if (w->kind() != Moose::VAR_SOLVER || &w->sys() != &_var.sys())
    paramError("w", "The chemical potential must be a nonlinear variable of the kernel's system.");
  if (_w_var == _var.number())
    paramError("w", "The chemical potential must differ from the kernel variable.");

  // The additional arguments of F: one d^2F/dc darg handle per entry (SplitCHParsed::_d2Fdcdarg).
  // The handle array is created even when the list is empty: a Kokkos array of property handles
  // is deep-copied by the dispatcher, and an array that was never allocated aborts that copy.
  const auto f = getMaterialPropertyName("f_name");
  _d2Fdcdarg.create(_n_args);
  std::set<unsigned int> seen;
  std::vector<VariableName> coupled_names;
  for (unsigned int i = 0; i < _n_args; ++i)
  {
    if (!isCoupled("coupled_variables", i))
      paramError("coupled_variables", "Entries must be variables, not constants.");
    const auto * arg = getVar("coupled_variables", i);
    const auto name = coupledName("coupled_variables", i);
    _arg_var.push_back(arg->number());
    coupled_names.push_back(name);

    // Only a nonlinear variable of this system has a Jacobian column (the host's JvarMap skips
    // auxiliary variables the same way); the duplicate check uses a host container because the
    // device map cannot be queried on the host before copyToDevice()
    if (arg->kind() == Moose::VAR_SOLVER && &arg->sys() == &_var.sys())
    {
      if (!seen.insert(arg->number()).second)
        paramError("coupled_variables", "Repeated arguments are not supported.");
      _arg_var_to_index[arg->number()] = i;
    }

    // An undeclared cross derivative stays an unset handle and contributes zero, as the host's
    // DerivativeMaterialInterface::getMaterialPropertyDerivative zero default does
    const auto derivative = derivativePropertyNameSecond(f, _var.name(), name);
    if (hasKokkosMaterialPropertyByName<Real>(derivative))
      _d2Fdcdarg[i] = getKokkosMaterialPropertyByName<Real>(derivative);
  }
  _arg_var_to_index.copyToDevice();

  // SplitCHParsed::initialSetup, validateNonlinearCoupling: every nonlinear variable whose cross
  // derivative with c the material declares must be coupled here (w counts as coupled on the
  // host, where it is in _coupled_standard_moose_vars); the host warns, so does this object
  std::vector<VariableName> missing;
  for (const auto & jname : _var.sys().getVariableNames())
  {
    if (jname == _var.name() || jname == w->name())
      continue;
    if (!hasKokkosMaterialPropertyByName<Real>(derivativePropertyNameSecond(f, _var.name(), jname)))
      continue;
    if (std::find(coupled_names.begin(), coupled_names.end(), jname) == coupled_names.end())
      missing.push_back(jname);
  }
  if (!missing.empty())
    mooseWarning("Missing coupled variables {",
                 Moose::stringify(missing),
                 "} (add them to coupled_variables parameter of ",
                 name(),
                 ")");

  // SplitCHParsed::initialSetup, validateDerivativeMaterialPropertyBase
  if (!hasKokkosMaterialPropertyByName<Real>(f))
    mooseWarning("The material property '",
                 f,
                 "' does not exist. The kernel '",
                 name(),
                 "' only needs its derivatives, but this may indicate a typo in the input file.");

  // A misspelled hook is silently replaced by the base default (KokkosDispatcher.h); fail loudly
  if (!Moose::Kokkos::hasUserJacobianHook<KokkosSplitCHParsed>() ||
      !Moose::Kokkos::hasUserOffDiagJacobianHook<KokkosSplitCHParsed>())
    mooseError(
        "KokkosSplitCHParsed Jacobian hooks were not registered with the expected signatures.");
}
