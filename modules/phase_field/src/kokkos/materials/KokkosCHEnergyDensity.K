//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosCHEnergyDensity.h"

#include "MooseVariable.h"

registerKokkosMaterial("PhaseFieldApp", KokkosCHEnergyDensity);

InputParameters
KokkosCHEnergyDensity::validParams()
{
  InputParameters params = Moose::Kokkos::Material::validParams();
  params.addClassDescription("Total free energy density f + kappa / 2 |grad(c)|^2 of one conserved "
                             "variable as a material property, for integration on the device with "
                             "KokkosElementIntegralMaterialProperty.");
  params.addParam<MaterialPropertyName>("f_name", "F", "Bulk free energy density property");
  params.addParam<MaterialPropertyName>(
      "kappa_name", "kappa_c", "Gradient energy coefficient property (or a literal value)");
  params.addRequiredCoupledVar("c", "Conserved variable whose gradient carries the interfacial energy");
  params.addParam<MaterialPropertyName>(
      "property_name", "f_total", "Name of the declared total free energy density property");
  return params;
}

KokkosCHEnergyDensity::KokkosCHEnergyDensity(const InputParameters & parameters)
  : Moose::Kokkos::Material(parameters),
    _f(getKokkosMaterialProperty<Real>("f_name")),
    _kappa(getKokkosMaterialProperty<Real>("kappa_name")),
    _grad_c(kokkosCoupledGradient("c")),
    _f_total(declareKokkosProperty<Real>("property_name"))
{
  if (!isCoupled("c") || coupledComponents("c") != 1)
    paramError("c", "Exactly one variable is required.");
}
