//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosOldCoupledForce.h"

#include "MooseVariable.h"

registerKokkosResidualObject("PhaseFieldApp", KokkosOldCoupledForce);

InputParameters
KokkosOldCoupledForce::validParams()
{
  InputParameters params = Moose::Kokkos::KernelValue::validParams();

  params.addClassDescription(
      "Implements a source term proportional to the value of a coupled variable at the previous "
      "time step, the Kokkos form of CoupledForce with implicit = false. Weak form: $(\\psi_i, "
      "-\\sigma v_{old})$. Contributes nothing to the Jacobian.");
  params.addRequiredCoupledVar("v", "The coupled variable which provides the force");
  params.addParam<Real>(
      "coef", 1.0, "Coefficent ($\\sigma$) multiplier for the coupled force term.");

  // The old-state read is built in; a Kokkos kernel would silently ignore this parameter
  params.suppressParameter<bool>("implicit");

  return params;
}

KokkosOldCoupledForce::KokkosOldCoupledForce(const InputParameters & parameters)
  : Moose::Kokkos::KernelValue(parameters),
    _v_var(coupled("v")),
    _v_old(kokkosCoupledValueOld("v")),
    _coef(getParam<Real>("coef"))
{
  if (_var.number() == _v_var)
    paramError("v",
               "Coupled variable 'v' needs to be different from 'variable' with "
               "KokkosOldCoupledForce.");

  // The residual depends on the old solution only: no Jacobian hook may ever be registered
  if (Moose::Kokkos::hasUserJacobianHook<KokkosOldCoupledForce>() ||
      Moose::Kokkos::hasUserOffDiagJacobianHook<KokkosOldCoupledForce>())
    mooseError("KokkosOldCoupledForce must not define a Jacobian hook.");
}
