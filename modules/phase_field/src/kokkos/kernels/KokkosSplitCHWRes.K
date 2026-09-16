//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosSplitCHWRes.h"

#include "MooseVariable.h"
#include "SystemBase.h"

#include <set>

registerKokkosResidualObject("PhaseFieldApp", KokkosSplitCHWRes);

InputParameters
KokkosSplitCHWRes::validParams()
{
  InputParameters params = Moose::Kokkos::KernelGrad::validParams();
  params.addClassDescription("Split formulation Cahn-Hilliard Kernel for the chemical potential "
                             "variable with a scalar (isotropic) mobility");
  params.addParam<MaterialPropertyName>("mob_name", "mobtemp", "The mobility used with the kernel");
  params.addCoupledVar("coupled_variables", "Vector of variable arguments of the mobility");
  params.addCoupledVar(
      "w", "Coupled chemical potential (if not specified kernel variable will be used)");
  return params;
}

KokkosSplitCHWRes::KokkosSplitCHWRes(const InputParameters & parameters)
  : Moose::Kokkos::KernelGrad(parameters),
    _mob(getKokkosMaterialProperty<Real>("mob_name")),
    _is_coupled(isCoupled("w")),
    _w_var(_is_coupled ? coupled("w") : _var.number()),
    _w_is_variable(!(_is_coupled && _w_var != _var.number())),
    _grad_w(_is_coupled ? kokkosCoupledGradient("w") : _grad_u),
    _n_args(coupledComponents("coupled_variables"))
{
  // A coupled w is restricted to a nonlinear variable of this kernel's system. This narrows the
  // host, which also runs an auxiliary or literal w (Coupleable::coupled() returns a non-aliasing
  // sentinel for either, Coupleable.C:112-116, 486-496); the aliasing hazard of trap T9 is the
  // host's JvarMap (JvarMapInterface.h:132-138), which the device argument map below avoids by
  // holding nonlinear variables of this system only. Kept as a deliberate restriction: no gate
  // deck and no ConservedAction deck couples w to anything else.
  if (_is_coupled)
  {
    const auto * w = getVar("w", 0);
    if (w->kind() != Moose::VAR_SOLVER || &w->sys() != &_var.sys())
      paramError("w", "The chemical potential must be a nonlinear variable of the kernel's system.");
  }

  // The mobility's arguments: one dM/darg handle per entry (SplitCHWResBase::_dmobdarg). The
  // handle array is created even when the list is empty: a Kokkos array of property handles is
  // deep-copied by the dispatcher, and an array that was never allocated aborts that copy.
  const auto mobility = getMaterialPropertyName("mob_name");
  _dmob_darg.create(_n_args);
  std::set<unsigned int> seen;
  for (unsigned int i = 0; i < _n_args; ++i)
  {
    if (!isCoupled("coupled_variables", i))
      paramError("coupled_variables", "Entries must be variables, not constants.");
    const auto * arg = getVar("coupled_variables", i);
    const auto name = coupledName("coupled_variables", i);
    _arg_var.push_back(arg->number());

    // Only a nonlinear variable of this system has a Jacobian column (the host's JvarMap skips
    // auxiliary variables the same way); the duplicate check uses a host container because the
    // device map cannot be queried on the host before copyToDevice()
    if (arg->kind() == Moose::VAR_SOLVER && &arg->sys() == &_var.sys())
    {
      if (!seen.insert(arg->number()).second)
        paramError("coupled_variables", "Repeated arguments are not supported.");
      _arg_var_to_index[arg->number()] = i;
    }

    // An undeclared derivative stays an unset handle and contributes zero, as the host's
    // DerivativeMaterialInterface::getMaterialPropertyDerivative zero default does
    const auto derivative = derivativePropertyNameFirst(mobility, name);
    if (hasKokkosMaterialPropertyByName<Real>(derivative))
      _dmob_darg[i] = getKokkosMaterialPropertyByName<Real>(derivative);
  }
  _arg_var_to_index.copyToDevice();

  // A misspelled hook is silently replaced by the base default (KokkosDispatcher.h); fail loudly
  if (!Moose::Kokkos::hasUserJacobianHook<KokkosSplitCHWRes>() ||
      !Moose::Kokkos::hasUserOffDiagJacobianHook<KokkosSplitCHWRes>())
    mooseError(
        "KokkosSplitCHWRes Jacobian hooks were not registered with the expected signatures.");
}
