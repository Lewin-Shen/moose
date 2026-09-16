//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosCHFreeEnergy.h"

#include "MooseVariable.h"

registerKokkosMaterial("PhaseFieldApp", KokkosCHFreeEnergy);

InputParameters
KokkosCHFreeEnergy::validParams()
{
  InputParameters params = Moose::Kokkos::Material::validParams();
  params.addClassDescription(
      "Cahn-Hilliard local free energy f(c) and its first two derivatives in closed form. Two "
      "hand-coded forms, not a general expression material: DOUBLE_WELL, W (c - ca)^2 (c - cb)^2, "
      "or REDLICH_KISTER, the regular-solution expression with logarithms of the tutorial "
      "(eVpJ/Vm*length_scale^3*(A*c+B*(1-c)+C*c*log(c)+D*(1-c)*log(1-c)+E*c*(1-c)"
      "+F*c*(1-c)*(2*c-1)+G*c*(1-c)*(2*c-1)^2)). Declares f, df/dc and d^2f/dc^2 under the "
      "DerivativeParsedMaterial property names.");
  params.addParam<MaterialPropertyName>(
      "property_name", "f_loc", "Name of the free energy property (derivative names derive from it)");
  params.addRequiredCoupledVar("c", "Concentration variable");

  MooseEnum form("DOUBLE_WELL REDLICH_KISTER", "DOUBLE_WELL");
  params.addParam<MooseEnum>("form", form, "Closed form of the free energy");

  // DOUBLE_WELL defaults: spinodal_polynomial.i (W = 16, ca = 0.27, cb = 0.83)
  params.addParam<Real>("W", 16.0, "DOUBLE_WELL: well height W in f = W (c - ca)^2 (c - cb)^2");
  params.addParam<Real>("ca", 0.27, "DOUBLE_WELL: first well position");
  params.addParam<Real>("cb", 0.83, "DOUBLE_WELL: second well position");

  // REDLICH_KISTER defaults: spinodal.i (the 500 K Fe-Cr regular-solution coefficients)
  params.addParam<Real>("A", -2.45e4, "REDLICH_KISTER: coefficient of c");
  params.addParam<Real>("B", -2.83e4, "REDLICH_KISTER: coefficient of (1 - c)");
  params.addParam<Real>("C", 4.17e3, "REDLICH_KISTER: coefficient of c ln(c)");
  params.addParam<Real>("D", 7.05e3, "REDLICH_KISTER: coefficient of (1 - c) ln(1 - c)");
  params.addParam<Real>("E", 1.21e4, "REDLICH_KISTER: coefficient of c (1 - c)");
  params.addParam<Real>("F", 2.57e3, "REDLICH_KISTER: coefficient of c (1 - c) (2c - 1)");
  params.addParam<Real>("G", -2.35e3, "REDLICH_KISTER: coefficient of c (1 - c) (2c - 1)^2");
  params.addParam<Real>("length_scale", 1e-9, "REDLICH_KISTER: length scale in m (prefactor)");
  params.addParam<Real>("eVpJ", 6.24150934e18, "REDLICH_KISTER: eV per J (prefactor)");
  params.addParam<Real>("Vm", 7.1e-6, "REDLICH_KISTER: molar volume in m^3/mol (prefactor)");
  return params;
}

KokkosCHFreeEnergy::KokkosCHFreeEnergy(const InputParameters & parameters)
  : Moose::Kokkos::Material(parameters),
    _form(getParam<MooseEnum>("form") == "REDLICH_KISTER" ? REDLICH_KISTER : DOUBLE_WELL),
    _W(getParam<Real>("W")),
    _ca(getParam<Real>("ca")),
    _cb(getParam<Real>("cb")),
    _A(getParam<Real>("A")),
    _B(getParam<Real>("B")),
    _C(getParam<Real>("C")),
    _D(getParam<Real>("D")),
    _E(getParam<Real>("E")),
    _F(getParam<Real>("F")),
    _G(getParam<Real>("G")),
    _p(getParam<Real>("eVpJ") / getParam<Real>("Vm") *
       (getParam<Real>("length_scale") * getParam<Real>("length_scale") *
        getParam<Real>("length_scale"))),
    _c(kokkosCoupledValue("c")),
    _f(declareKokkosProperty<Real>("property_name"))
{
  // The derivative property names are built from the variable name, so 'c' must be a variable
  // (auxiliary or nonlinear), not a literal default value
  if (!isCoupled("c") || coupledComponents("c") != 1)
    paramError("c", "Exactly one variable is required (a literal value is not supported).");

  const MaterialPropertyName f = getParam<MaterialPropertyName>("property_name");
  const auto cname = coupledName("c", 0);

  // DerivativeParsedMaterial declares 'd<f>/d<c>' and 'd^2<f>/d<c>^2' through
  // DerivativeMaterialPropertyNameInterface; the same helper gives the same strings here
  _df = declareKokkosPropertyByName<Real>(derivativePropertyNameFirst(f, cname));
  _d2f = declareKokkosPropertyByName<Real>(derivativePropertyNameSecond(f, cname, cname));
}
