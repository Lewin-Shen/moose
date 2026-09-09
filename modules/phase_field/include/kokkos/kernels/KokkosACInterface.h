//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosKernel.h"
#include "DerivativeMaterialPropertyNameInterface.h"

/**
 * Kokkos specialization of ACInterface (the gradient-energy Allen-Cahn term):
 *
 *   R_i = grad(u) . kappa (L grad(psi_i) + psi_i grad(L))
 *
 * where grad(L) = sum_c dL/darg_c grad(arg_c) over an optional prescribed auxiliary
 * temperature (variable_L = true). L and kappa must not depend on the nonlinear order
 * parameters: dL/dop, dkappa/dop and the second derivatives ACInterface carries are all zero in
 * this specialization, so the diagonal Jacobian reduces to grad(phi_j) . kappa grad(L psi_i) and
 * every off-diagonal block is zero (the auxiliary T has no nonlinear column). Declared
 * derivatives outside this scope are rejected at construction, and a material that declares dL/dT
 * while T is not coupled here is rejected too (the gradient term would otherwise be dropped).
 */
class KokkosACInterface : public Moose::Kokkos::Kernel,
                          public DerivativeMaterialPropertyNameInterface
{
  using Real3 = Moose::Kokkos::Real3;

public:
  static InputParameters validParams();

  KokkosACInterface(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION Real computeQpResidual(const unsigned int i,
                                         const unsigned int qp,
                                         AssemblyDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION Real computeQpJacobian(const unsigned int i,
                                         const unsigned int j,
                                         const unsigned int qp,
                                         AssemblyDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION Real computeQpOffDiagJacobian(const unsigned int i,
                                                const unsigned int j,
                                                const unsigned int jvar,
                                                const unsigned int qp,
                                                AssemblyDatum & datum) const;

private:
  /// kappa * (L grad(psi_i) + psi_i grad(L)), ACInterface::kappaNablaLPsi
  KOKKOS_FUNCTION Real3 kappaNablaLPsi(const unsigned int i,
                                       const unsigned int qp,
                                       AssemblyDatum & datum) const
  {
    Real3 sum = _L(datum, qp) * _grad_test(datum, i, qp);
    // grad(L) = dL/dT grad(T) for the prescribed auxiliary T (zero when no derivative is declared)
    if (_variable_L && _dL_dT)
      sum += _grad_T(datum, qp) * (_dL_dT(datum, qp) * _test(datum, i, qp));
    return _kappa(datum, qp) * sum;
  }

  /// Mobility L and gradient-energy coefficient kappa
  const Moose::Kokkos::MaterialProperty<Real> _L;
  const Moose::Kokkos::MaterialProperty<Real> _kappa;
  /// Whether grad(L) is carried (ACInterface's variable_L)
  const bool _variable_L;
  /// Number of coupled_variables entries (0 or 1: the prescribed auxiliary T)
  const unsigned int _n_args;
  /// Gradient of the prescribed auxiliary T (unset when coupled_variables is empty)
  Moose::Kokkos::VariableGradient _grad_T;
  /// dL/dT handle; unset (false) when the material declares no such derivative
  Moose::Kokkos::MaterialProperty<Real> _dL_dT;
};

template <typename Derived>
KOKKOS_FUNCTION Real
KokkosACInterface::computeQpResidual(const unsigned int i,
                                     const unsigned int qp,
                                     AssemblyDatum & datum) const
{
  return _grad_u(datum, qp) * kappaNablaLPsi(i, qp, datum);
}

template <typename Derived>
KOKKOS_FUNCTION Real
KokkosACInterface::computeQpJacobian(const unsigned int i,
                                     const unsigned int j,
                                     const unsigned int qp,
                                     AssemblyDatum & datum) const
{
  // ACInterface::computeQpJacobian with dkappa/dop = dL/dop = d2L/dop2 = d2L/darg dop = 0
  return _grad_phi(datum, j, qp) * kappaNablaLPsi(i, qp, datum);
}

template <typename Derived>
KOKKOS_FUNCTION Real
KokkosACInterface::computeQpOffDiagJacobian(const unsigned int /* i */,
                                            const unsigned int /* j */,
                                            const unsigned int /* jvar */,
                                            const unsigned int /* qp */,
                                            AssemblyDatum & /* datum */) const
{
  // No nonlinear coupled argument exists in this specialization (T is auxiliary)
  return 0.0;
}
