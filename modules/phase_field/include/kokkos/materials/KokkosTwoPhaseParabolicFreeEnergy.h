//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosMaterial.h"
#include "DerivativeMaterialPropertyNameInterface.h"

/**
 * Kokkos material for the two-phase parabolic free energy of a concentration c and n order
 * parameters eta_1..eta_n, with its first and second derivatives in closed form:
 *
 *   f = rho_sq (c - ca)^2 (1 - h) + rho_sq (cb - c)^2 h + w g
 *   h = sum_k p(eta_k),  p(e) = e^3 (6 e^2 - 15 e + 10)
 *   g = sum_k eta_k^2 (1 - eta_k)^2 + 2 alpha sum_{k<m} eta_k^2 eta_m^2
 *
 * Two parabolic phase energies interpolated by the sum of the smoothstep functions of the order
 * parameters, plus the multi-well with the pairwise penalty; w is the barrier height (not a
 * chemical potential). This is NOT a general expression material: the form and its derivatives are
 * hand-written and compiled ahead of time (the Kokkos function parser evaluates values only and
 * cannot differentiate). The defaults are the constants of PFHub Benchmark 2.
 *
 * Declared properties, (n + 2)(n + 3) / 2 per quadrature point, under the names a
 * DerivativeParsedMaterial with derivative_order = 2 declares (DerivativeMaterialPropertyNameInterface
 * on the host, sorted denominators, so with ten or more order parameters eta10 sorts before eta2):
 *   f, df/dc, d^2f/dc^2 (the constant 2 rho_sq, declared and written at every point),
 *   df/deta_k, d^2f/dcdeta_k, d^2f/deta_k^2 for every k, and
 *   d^2f/deta_kdeta_m for every pair k < m, stored at the packed upper-triangle index
 *   k n - k (k + 1) / 2 + m - k - 1 (0-based positions in 'etas').
 * The pair terms are always declared, also where alpha or w is zero (the parsed material omits a
 * derivative that simplifies to zero; here it is declared and written as zero, which a kernel reads
 * the same way). No third derivatives are declared.
 */
class KokkosTwoPhaseParabolicFreeEnergy : public Moose::Kokkos::Material,
                                          public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  KokkosTwoPhaseParabolicFreeEnergy(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION void computeQpProperties(const unsigned int qp, Datum & datum) const;

private:
  /// Packed upper-triangle index of the pair (k, m), 0 <= k < m < n
  KOKKOS_FUNCTION unsigned int pairIndex(const unsigned int k, const unsigned int m) const
  {
    return k * _n - k * (k + 1) / 2 + m - k - 1;
  }

  /// Number of order parameters
  const unsigned int _n;

  /// rho_sq (enters linearly), the two phase compositions, the barrier height and the penalty
  const Real _rho_sq;
  const Real _ca;
  const Real _cb;
  const Real _w;
  const Real _alpha;

  /// Concentration and the order parameters (component k is eta_k); set after the guards
  Moose::Kokkos::VariableValue _c;
  Moose::Kokkos::VariableValue _etas;

  /// f, df/dc, d^2f/dc^2
  Moose::Kokkos::MaterialProperty<Real> _f;
  Moose::Kokkos::MaterialProperty<Real> _dfdc;
  Moose::Kokkos::MaterialProperty<Real> _d2fdc2;

  /// Per order parameter: df/deta_k, d^2f/dcdeta_k, d^2f/deta_k^2 (each created with n)
  Moose::Kokkos::Array<Moose::Kokkos::MaterialProperty<Real>> _dfdeta;
  Moose::Kokkos::Array<Moose::Kokkos::MaterialProperty<Real>> _d2fdcdeta;
  Moose::Kokkos::Array<Moose::Kokkos::MaterialProperty<Real>> _d2fdeta2;
  /// Per pair k < m: d^2f/deta_kdeta_m (created with n (n - 1) / 2, empty when n = 1)
  Moose::Kokkos::Array<Moose::Kokkos::MaterialProperty<Real>> _d2fdetadeta;
};

template <typename Derived>
KOKKOS_FUNCTION void
KokkosTwoPhaseParabolicFreeEnergy::computeQpProperties(const unsigned int qp, Datum & datum) const
{
  // A = (c - ca)^2, B = (cb - c)^2, D = B - A; one pass over the order parameters for
  // h = sum p(e), Q = sum e^2 (1 - e)^2, S = sum e^2, S4 = sum e^4, then P = (S^2 - S4) / 2
  const Real c = _c(datum, qp);
  const Real A = (c - _ca) * (c - _ca), B = (_cb - c) * (_cb - c), D = B - A, dcab = _cb - _ca;
  Real h = 0, S = 0, S4 = 0, Q = 0;
  for (unsigned int k = 0; k < _n; ++k)
  {
    const Real e = _etas(datum, qp, k), e2 = e * e, ome = 1 - e;
    h += e2 * e * (6 * e2 - 15 * e + 10);
    Q += e2 * ome * ome;
    S += e2;
    S4 += e2 * e2;
  }
  const Real P = 0.5 * (S * S - S4);
  _f(datum, qp) = _rho_sq * (A * (1 - h) + B * h) + _w * (Q + 2 * _alpha * P);
  _dfdc(datum, qp) = 2 * _rho_sq * ((c - _ca) - h * dcab);
  _d2fdc2(datum, qp) = 2 * _rho_sq;

  // Per order parameter, with S_k = S - e^2, p' = 30 e^2 (1 - e)^2, p'' = 60 e (1 - e)(1 - 2e),
  // q' = 2 e (1 - e)(1 - 2e), q'' = 2 (1 - 6e + 6e^2); the pair loop gives 8 w alpha eta_k eta_m
  for (unsigned int k = 0; k < _n; ++k)
  {
    const Real e = _etas(datum, qp, k), ome = 1 - e, Sk = S - e * e;
    const Real dp = 30 * e * e * ome * ome, d2p = 60 * e * ome * (1 - 2 * e);
    const Real dq = 2 * e * ome * (1 - 2 * e), d2q = 2 * (1 - 6 * e + 6 * e * e);
    _dfdeta[k](datum, qp) = _rho_sq * D * dp + _w * (dq + 4 * _alpha * e * Sk);
    _d2fdcdeta[k](datum, qp) = -2 * _rho_sq * dcab * dp;
    _d2fdeta2[k](datum, qp) = _rho_sq * D * d2p + _w * (d2q + 4 * _alpha * Sk);
    for (unsigned int m = k + 1; m < _n; ++m)
      _d2fdetadeta[pairIndex(k, m)](datum, qp) = 8 * _w * _alpha * e * _etas(datum, qp, m);
  }
}
