//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosKernelValue.h"
#include "KokkosMap.h"
#include "DerivativeMaterialPropertyNameInterface.h"

/**
 * Kokkos specialization of ACGrGrPoly (the polycrystalline Allen-Cahn bulk term of the
 * GBEvolution grain-growth model):
 *
 *   R_i = L mu (op^3 - op + 2 gamma op sum_j op_j^2) psi_i
 *
 * with hand-coded diagonal and off-diagonal (other order parameter) Jacobians. The physics is
 * transcribed from ACGrGrPoly::computeDFDOP / computeQpOffDiagJacobian and ACBulk's product
 * rule. L, mu and gamma_asymm must be independent of every nonlinear order parameter; a
 * declared dependence (a d<prop>/d<op> material property) is rejected at construction, so the
 * dL/dop product-rule terms ACBulk carries are omitted by construction, not by accident.
 * An optional auxiliary temperature 'T' may be listed in coupled_variables for action
 * compatibility; it contributes no bulk term (the mobility reads T through its material).
 */
class KokkosACGrGrPoly : public Moose::Kokkos::KernelValue,
                         public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  KokkosACGrGrPoly(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION Real precomputeQpResidual(const unsigned int qp, AssemblyDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION Real precomputeQpJacobian(const unsigned int j,
                                            const unsigned int qp,
                                            AssemblyDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION Real precomputeQpOffDiagJacobian(const unsigned int j,
                                                   const unsigned int jvar,
                                                   const unsigned int qp,
                                                   AssemblyDatum & datum) const;

private:
  /// sum_j op_j^2 over the OTHER order parameters (ACGrGrPoly's SumOPj)
  KOKKOS_FUNCTION Real sumOtherOpsSquared(const unsigned int qp, AssemblyDatum & datum) const
  {
    Real sum = 0.0;
    for (unsigned int c = 0; c < _op_num; ++c)
    {
      const Real other = _v(datum, qp, c);
      sum += other * other;
    }
    return sum;
  }

  /// Number of OTHER order parameters coupled through 'v' (N - 1)
  const unsigned int _op_num;
  /// Host-only metadata: variable numbers of the coupled order parameters, in 'v' order
  std::vector<unsigned int> _v_var;
  /// Component-indexed values of the coupled order parameters
  const Moose::Kokkos::VariableValue _v;
  /// Variable number -> component index in 'v' (device-usable map)
  Moose::Kokkos::Map<unsigned int, unsigned int> _v_var_to_index;
  /// Mobility L, and the GBEvolution properties mu and gamma_asymm
  const Moose::Kokkos::MaterialProperty<Real> _L;
  const Moose::Kokkos::MaterialProperty<Real> _mu;
  const Moose::Kokkos::MaterialProperty<Real> _gamma;
};

template <typename Derived>
KOKKOS_FUNCTION Real
KokkosACGrGrPoly::precomputeQpResidual(const unsigned int qp, AssemblyDatum & datum) const
{
  // ACBulk::precomputeQpResidual = L * dF/dop, with ACGrGrPoly::computeDFDOP(Residual)
  const Real op = _u(datum, qp);
  const Real sum = sumOtherOpsSquared(qp, datum);
  return _L(datum, qp) * _mu(datum, qp) *
         (op * op * op - op + 2.0 * _gamma(datum, qp) * op * sum);
}

template <typename Derived>
KOKKOS_FUNCTION Real
KokkosACGrGrPoly::precomputeQpJacobian(const unsigned int j,
                                       const unsigned int qp,
                                       AssemblyDatum & datum) const
{
  // ACBulk::precomputeQpJacobian = L * dF/dop(Jacobian); the dL/dop * phi_j * dF/dop term is
  // zero for every material this specialization accepts (see the constructor guards)
  const Real op = _u(datum, qp);
  const Real sum = sumOtherOpsSquared(qp, datum);
  return _L(datum, qp) * _mu(datum, qp) * _phi(datum, j, qp) *
         (3.0 * op * op - 1.0 + 2.0 * _gamma(datum, qp) * sum);
}

template <typename Derived>
KOKKOS_FUNCTION Real
KokkosACGrGrPoly::precomputeQpOffDiagJacobian(const unsigned int j,
                                              const unsigned int jvar,
                                              const unsigned int qp,
                                              AssemblyDatum & datum) const
{
  // ACGrGrPoly::computeQpOffDiagJacobian divided by the test function (KernelValue applies it):
  // L * mu * 2 gamma op * dSumOPj, dSumOPj = 2 op_c phi_j for the column's order parameter
  if (!_v_var_to_index.exists(jvar))
    return 0.0;
  const unsigned int c = _v_var_to_index[jvar];
  const Real dsum = 2.0 * _v(datum, qp, c) * _phi(datum, j, qp);
  return _L(datum, qp) * _mu(datum, qp) * 2.0 * _gamma(datum, qp) * _u(datum, qp) * dsum;
}
