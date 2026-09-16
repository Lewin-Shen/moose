//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

/**
 * Total free energy density of a single conserved variable c, as a material property:
 *
 *   f_total = f(c) + kappa / 2 |grad(c)|^2
 *
 * the formula of the TotalFreeEnergy aux kernel for one interfacial variable, declared as a material
 * property so that ElementIntegralMaterialProperty integrates it at the quadrature points. Host twin of
 * KokkosCHEnergyDensity, so that a CPU control deck forms and integrates the energy exactly as the
 * Kokkos deck does on the device. A test-app object (the CPU controls are tests; the product object is
 * the Kokkos one), so decks using it need --allow-test-objects.
 */
class CHEnergyDensity : public Material
{
public:
  static InputParameters validParams();

  CHEnergyDensity(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Bulk free energy density f(c)
  const MaterialProperty<Real> & _f;
  /// Gradient energy coefficient kappa
  const MaterialProperty<Real> & _kappa;
  /// grad(c)
  const VariableGradient & _grad_c;
  /// f_total
  MaterialProperty<Real> & _f_total;
};
