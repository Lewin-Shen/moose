//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosCHMobility.h"

#include "MooseVariable.h"

registerKokkosMaterial("PhaseFieldApp", KokkosCHMobility);

InputParameters
KokkosCHMobility::validParams()
{
  InputParameters params = Moose::Kokkos::Material::validParams();
  params.addClassDescription("Cahn-Hilliard mobility M and its derivative dM/dc in closed form: "
                             "CONSTANT (M = M0) or DEGENERATE (M = M0 c (1 - c)).");
  params.addParam<MaterialPropertyName>(
      "property_name", "M", "Name of the mobility property (the derivative name derives from it)");
  params.addRequiredCoupledVar("c", "Concentration variable");

  MooseEnum form("CONSTANT DEGENERATE", "CONSTANT");
  params.addParam<MooseEnum>("form", form, "Closed form of the mobility");

  // Default: the tutorial's mobility (spinodal.i, [Materials/mobility] prop_values = 2.2e-5)
  params.addParam<Real>("M0", 2.2e-5, "Mobility prefactor M0");
  return params;
}

KokkosCHMobility::KokkosCHMobility(const InputParameters & parameters)
  : Moose::Kokkos::Material(parameters),
    _form(getParam<MooseEnum>("form") == "DEGENERATE" ? DEGENERATE : CONSTANT),
    _M0(getParam<Real>("M0")),
    _c(kokkosCoupledValue("c")),
    _M(declareKokkosProperty<Real>("property_name"))
{
  if (!isCoupled("c") || coupledComponents("c") != 1)
    paramError("c", "Exactly one variable is required (a literal value is not supported).");

  const MaterialPropertyName M = getParam<MaterialPropertyName>("property_name");
  _dMdc = declareKokkosPropertyByName<Real>(derivativePropertyNameFirst(M, coupledName("c", 0)));
}
