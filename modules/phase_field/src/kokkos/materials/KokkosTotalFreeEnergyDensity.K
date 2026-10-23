//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosTotalFreeEnergyDensity.h"

#include "MooseVariable.h"

registerKokkosMaterial("PhaseFieldApp", KokkosTotalFreeEnergyDensity);

InputParameters
KokkosTotalFreeEnergyDensity::validParams()
{
  InputParameters params = Moose::Kokkos::Material::validParams();
  params.addClassDescription(
      "Total free energy density (both the bulk and gradient parts) as a material property, where "
      "the bulk free energy has been defined in a material, for integration on the device with "
      "KokkosElementIntegralMaterialProperty");
  params.addCoupledVar("interfacial_vars", "Variable names that contribute to interfacial energy");
  params.addCoupledVar(
      "additional_free_energy",
      0.0,
      "Coupled variable holding additional free energy contributions to be summed up");
  params.addParam<MaterialPropertyName>("f_name", "F", " Base name of the free energy function");
  params.addParam<std::vector<MaterialPropertyName>>("kappa_names",
                                                     std::vector<MaterialPropertyName>(),
                                                     "Vector of kappa names corresponding to "
                                                     "each variable name in interfacial_vars "
                                                     "in the same order.");
  params.addParam<MaterialPropertyName>(
      "property_name", "f_density", "Name of the declared total free energy density property");
  return params;
}

KokkosTotalFreeEnergyDensity::KokkosTotalFreeEnergyDensity(const InputParameters & parameters)
  : Moose::Kokkos::Material(parameters),
    _nvars(coupledComponents("interfacial_vars")),
    _F(getKokkosMaterialProperty<Real>("f_name")),
    _grad_vars(kokkosCoupledGradients("interfacial_vars")),
    _additional_free_energy(kokkosCoupledValue("additional_free_energy")),
    _total_energy(declareKokkosProperty<Real>("property_name"))
{
  // Error check to ensure size of interfacial_vars is the same as kappa_names (TotalFreeEnergy)
  const auto & kappa_names = getParam<std::vector<MaterialPropertyName>>("kappa_names");
  if (_nvars != kappa_names.size())
    paramError("kappa_names",
               "Size of interfacial_vars is not equal to the size of kappa_names in "
               "KokkosTotalFreeEnergyDensity");

  // Assign kappa values by property name, paired with interfacial_vars by position
  _kappas.create(_nvars);
  for (unsigned int i = 0; i < _nvars; ++i)
    _kappas[i] = getKokkosMaterialPropertyByName<Real>(kappa_names[i]);
}
