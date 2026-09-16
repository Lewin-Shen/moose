//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosKernelValue.h"

/**
 * Kokkos specialization of CoupledForce with implicit = false: a source term proportional to the
 * value of a coupled variable at the previous time step. Weak form: (psi_i, -sigma v_old).
 *
 * The CPU kernel reads the old solution through Coupleable::coupledValue() when its 'implicit'
 * parameter is false (Coupleable.C:540) and every Jacobian loop skips it
 * (ComputeFullJacobianThread.C:62, 94, 117). A Kokkos kernel does not honor 'implicit': the
 * parameter is inherited from TransientInterface through ResidualObject::validParams()
 * (ResidualObject.C:18) but nothing under framework/{include,src}/kokkos reads it, so this object
 * reads the old state explicitly with kokkosCoupledValueOld() and suppresses the parameter.
 *
 * No Jacobian hook is defined, on purpose: the residual depends on the old solution only, so the
 * diagonal and every off-diagonal block are structurally zero. The Kokkos dispatcher detects hooks
 * by function-pointer comparison (KokkosDispatcher.h:441-462), so an undefined hook and a misspelled
 * one look the same; the constructor asserts that none is registered.
 */
class KokkosOldCoupledForce : public Moose::Kokkos::KernelValue
{
public:
  static InputParameters validParams();

  KokkosOldCoupledForce(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION Real precomputeQpResidual(const unsigned int qp, AssemblyDatum & datum) const;

private:
  /// Coupled variable number
  const unsigned int _v_var;
  /// Coupled variable value at the previous time step (OLD_SOLUTION_TAG)
  const Moose::Kokkos::VariableValue _v_old;
  /// Multiplier for the coupled force term
  const Real _coef;
};

template <typename Derived>
KOKKOS_FUNCTION Real
KokkosOldCoupledForce::precomputeQpResidual(const unsigned int qp, AssemblyDatum & datum) const
{
  // CoupledForce::computeQpResidual, -_coef * _v[_qp], with _v the old solution (implicit = false)
  return -_coef * _v_old(datum, qp);
}
