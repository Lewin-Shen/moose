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
 * Total free energy density (bulk plus gradient parts) over several interfacial variables, as a
 * Kokkos material property:
 *
 *   f_density = F + additional_free_energy + sum_i kappa_i / 2 |grad(v_i)|^2
 *
 * the formula of the TotalFreeEnergy aux kernel, declared as a material property so that
 * KokkosElementIntegralMaterialProperty integrates it and KokkosMaterialRealAux writes it to an
 * elemental field on the device. A host aux kernel or postprocessor over a Kokkos material property
 * reads a silent zero, so the total free energy of a Kokkos deck must be formed and integrated on
 * the device. kappa_names[i] is paired with interfacial_vars[i] by position, as on the CPU; each
 * entry is a material property name (a plain number is read as a constant). An uncoupled
 * additional_free_energy reads its default value, 0, as on the CPU.
 */
class KokkosTotalFreeEnergyDensity : public Moose::Kokkos::Material
{
public:
  static InputParameters validParams();

  KokkosTotalFreeEnergyDensity(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION void computeQpProperties(const unsigned int qp, Datum & datum) const;

private:
  /// Number of interfacial variables (equal to the number of kappa names)
  const unsigned int _nvars;
  /// Bulk free energy density F
  const Moose::Kokkos::MaterialProperty<Real> _F;
  /// Gradients of the interfacial variables (component i is v_i)
  const Moose::Kokkos::VariableGradient _grad_vars;
  /// Additional free energy contribution (a constant default value when uncoupled)
  const Moose::Kokkos::VariableValue _additional_free_energy;
  /// Gradient energy coefficient of each interfacial variable; always created, possibly empty,
  /// because the dispatcher deep-copies this array and aborts on an unallocated one
  Moose::Kokkos::Array<Moose::Kokkos::MaterialProperty<Real>> _kappas;
  /// The declared total free energy density
  Moose::Kokkos::MaterialProperty<Real> _total_energy;
};

template <typename Derived>
KOKKOS_FUNCTION void
KokkosTotalFreeEnergyDensity::computeQpProperties(const unsigned int qp, Datum & datum) const
{
  // TotalFreeEnergy::computeValue: bulk energy contribution and additional contributions
  const Real F = _F(datum, qp);
  const Real additional = _additional_free_energy(datum, qp);
  Real total_energy = F + additional;

  // Interfacial energy of each variable: kappa_i / 2 |grad(v_i)|^2
  for (unsigned int i = 0; i < _nvars; ++i)
  {
    const Real kappa = _kappas[i](datum, qp);
    const Moose::Kokkos::Real3 grad = _grad_vars(datum, qp, i);
    total_energy += kappa / 2.0 * (grad * grad);
  }

  _total_energy(datum, qp) = total_energy;
}
