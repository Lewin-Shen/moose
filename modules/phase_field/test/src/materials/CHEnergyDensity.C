//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CHEnergyDensity.h"

registerMooseObject("PhaseFieldTestApp", CHEnergyDensity);

InputParameters
CHEnergyDensity::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Total free energy density f + kappa / 2 |grad(c)|^2 of one conserved "
                             "variable as a material property, for integration with "
                             "ElementIntegralMaterialProperty (host twin of KokkosCHEnergyDensity).");
  params.addParam<MaterialPropertyName>("f_name", "F", "Bulk free energy density property");
  params.addParam<MaterialPropertyName>(
      "kappa_name", "kappa_c", "Gradient energy coefficient property (or a literal value)");
  params.addRequiredCoupledVar("c", "Conserved variable whose gradient carries the interfacial energy");
  params.addParam<MaterialPropertyName>(
      "property_name", "f_total", "Name of the declared total free energy density property");
  return params;
}

CHEnergyDensity::CHEnergyDensity(const InputParameters & parameters)
  : Material(parameters),
    _f(getMaterialProperty<Real>("f_name")),
    _kappa(getMaterialProperty<Real>("kappa_name")),
    _grad_c(coupledGradient("c")),
    _f_total(declareProperty<Real>(getParam<MaterialPropertyName>("property_name")))
{
}

void
CHEnergyDensity::computeQpProperties()
{
  // TotalFreeEnergy::computeValue for one interfacial variable: F + kappa / 2 |grad(c)|^2
  _f_total[_qp] = _f[_qp] + _kappa[_qp] / 2.0 * _grad_c[_qp].norm_sq();
}
