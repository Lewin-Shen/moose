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
 * Test-only two-variable free energy for exercising the cross-derivative branches of
 * KokkosMatDiffusion and KokkosCoupledMaterialDerivative (implementation decision I-10):
 *
 *   f = W (c - ca)^2 (c - cb)^2 + a c eta^2
 *
 * Declared properties (DerivativeMaterialPropertyNameInterface names, sorted denominators):
 * f, df/dc, d^2f/dc^2, df/deta, d^2f/deta^2, d^2f/dcdeta. The CPU twin is a
 * DerivativeParsedMaterial with the same expression and coupled_variables = 'c eta'.
 */
class KokkosCHFreeEnergyCoupledTest : public Moose::Kokkos::Material,
                                      public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  KokkosCHFreeEnergyCoupledTest(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION void computeQpProperties(const unsigned int qp, Datum & datum) const;

private:
  const Real _W;
  const Real _ca;
  const Real _cb;
  const Real _a;

  const Moose::Kokkos::VariableValue _c;
  const Moose::Kokkos::VariableValue _eta;

  Moose::Kokkos::MaterialProperty<Real> _f;
  Moose::Kokkos::MaterialProperty<Real> _df_dc;
  Moose::Kokkos::MaterialProperty<Real> _d2f_dc2;
  Moose::Kokkos::MaterialProperty<Real> _df_deta;
  Moose::Kokkos::MaterialProperty<Real> _d2f_deta2;
  Moose::Kokkos::MaterialProperty<Real> _d2f_dcdeta;
};

template <typename Derived>
KOKKOS_FUNCTION void
KokkosCHFreeEnergyCoupledTest::computeQpProperties(const unsigned int qp, Datum & datum) const
{
  const Real c = _c(datum, qp);
  const Real eta = _eta(datum, qp);
  const Real p = c - _ca;
  const Real q = c - _cb;

  _f(datum, qp) = _W * p * p * q * q + _a * c * eta * eta;
  _df_dc(datum, qp) = 2.0 * _W * p * q * (p + q) + _a * eta * eta;
  _d2f_dc2(datum, qp) = 2.0 * _W * (p * p + q * q + 4.0 * p * q);
  _df_deta(datum, qp) = 2.0 * _a * c * eta;
  _d2f_deta2(datum, qp) = 2.0 * _a * c;
  _d2f_dcdeta(datum, qp) = 2.0 * _a * eta;
}
