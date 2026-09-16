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

/**
 * Total free energy density of a single conserved variable c, as a Kokkos material property:
 *
 *   f_total = f(c) + kappa / 2 |grad(c)|^2
 *
 * the formula of the TotalFreeEnergy aux kernel for one interfacial variable, declared as a material
 * property so that KokkosElementIntegralMaterialProperty integrates it on the device. A host aux
 * kernel or postprocessor over a Kokkos material property reads a silent zero, so the total free
 * energy of a Kokkos deck must be formed and integrated on the device. The host twin is CHEnergyDensity.
 */
class KokkosCHEnergyDensity : public Moose::Kokkos::Material
{
public:
  static InputParameters validParams();

  KokkosCHEnergyDensity(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION void computeQpProperties(const unsigned int qp, Datum & datum) const;

private:
  /// Bulk free energy density f(c)
  const Moose::Kokkos::MaterialProperty<Real> _f;
  /// Gradient energy coefficient kappa
  const Moose::Kokkos::MaterialProperty<Real> _kappa;
  /// grad(c)
  const Moose::Kokkos::VariableGradient _grad_c;
  /// f_total
  Moose::Kokkos::MaterialProperty<Real> _f_total;
};

template <typename Derived>
KOKKOS_FUNCTION void
KokkosCHEnergyDensity::computeQpProperties(const unsigned int qp, Datum & datum) const
{
  // TotalFreeEnergy::computeValue for one interfacial variable: F + kappa / 2 |grad(c)|^2
  const Real f = _f(datum, qp);
  const Real kappa = _kappa(datum, qp);
  const Moose::Kokkos::Real3 grad_c = _grad_c(datum, qp);
  _f_total(datum, qp) = f + kappa / 2.0 * (grad_c * grad_c);
}
