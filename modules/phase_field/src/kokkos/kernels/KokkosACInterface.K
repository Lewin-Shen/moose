//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosACInterface.h"

#include "MooseVariable.h"

registerKokkosResidualObject("PhaseFieldApp", KokkosACInterface);

InputParameters
KokkosACInterface::validParams()
{
  InputParameters params = Moose::Kokkos::Kernel::validParams();
  params.addClassDescription("Gradient energy Allen-Cahn Kernel");
  params.addCoupledVar("coupled_variables",
                       "Vector of nonlinear variable arguments this object depends on");
  params.addParam<MaterialPropertyName>("mob_name", "L", "The mobility used with the kernel");
  params.addParam<MaterialPropertyName>("kappa_name", "kappa_op", "The kappa used with the kernel");
  params.addParam<bool>("variable_L",
                        true,
                        "The mobility is a function of any MOOSE variable (if "
                        "this is set to false L must be constant over the "
                        "entire domain!)");
  return params;
}

KokkosACInterface::KokkosACInterface(const InputParameters & parameters)
  : Moose::Kokkos::Kernel(parameters),
    _L(getKokkosMaterialProperty<Real>("mob_name")),
    _kappa(getKokkosMaterialProperty<Real>("kappa_name")),
    _variable_L(getParam<bool>("variable_L")),
    _n_args(coupledComponents("coupled_variables"))
{
  if (_n_args > 1)
    paramError("coupled_variables", "Only the prescribed auxiliary T is supported.");

  // Resolve the base property names first, then reject every derivative the general CPU
  // formula would use: this class must not silently acquire nonlinear-temperature physics
  const auto mobility = getMaterialPropertyName("mob_name");
  const auto kappa = getMaterialPropertyName("kappa_name");
  const auto reject = [this](const MaterialPropertyName & derivative)
  {
    if (hasKokkosMaterialProperty<Real>(derivative))
      mooseError("Unsupported material derivative in KokkosACInterface: ", derivative);
  };
  reject(derivativePropertyNameFirst(mobility, _var.name()));
  reject(derivativePropertyNameFirst(kappa, _var.name()));
  reject(derivativePropertyNameSecond(mobility, _var.name(), _var.name()));

  // A material that declares dL/dT (KokkosGBEvolution with an auxiliary T) carries a gradient term this
  // kernel would silently drop if T is not coupled here; the CPU kernel validates that coupling
  if (_variable_L && !_n_args)
  {
    const auto dLdT = derivativePropertyNameFirst(mobility, "T");
    if (hasKokkosMaterialProperty<Real>(dLdT))
      mooseError("KokkosACInterface: the mobility material declares ",
                 dLdT,
                 "; couple T through coupled_variables, or set variable_L = false only if the "
                 "mobility is spatially constant");
  }

  if (_n_args)
  {
    if (!isCoupled("coupled_variables", 0))
      paramError("coupled_variables", "Use the auxiliary field T, or omit this list.");
    const auto * arg = getVar("coupled_variables", 0);
    if (arg->name() != "T" || arg->kind() != Moose::VAR_AUXILIARY)
      paramError("coupled_variables", "Only the prescribed auxiliary T is supported.");

    const auto name = coupledName("coupled_variables", 0);
    reject(derivativePropertyNameFirst(kappa, name));
    reject(derivativePropertyNameSecond(mobility, name, _var.name()));
    reject(derivativePropertyNameSecond(mobility, name, name));

    // A missing first derivative leaves a zero contribution (the handle stays unset)
    _grad_T = kokkosCoupledGradient("coupled_variables", 0);
    const auto derivative = derivativePropertyNameFirst(mobility, name);
    if (hasKokkosMaterialProperty<Real>(derivative))
      _dL_dT = getKokkosMaterialPropertyByName<Real>(derivative);
  }

  if (!Moose::Kokkos::hasUserJacobianHook<KokkosACInterface>() ||
      !Moose::Kokkos::hasUserOffDiagJacobianHook<KokkosACInterface>())
    mooseError(
        "KokkosACInterface Jacobian hooks were not registered with the expected signatures.");
}
