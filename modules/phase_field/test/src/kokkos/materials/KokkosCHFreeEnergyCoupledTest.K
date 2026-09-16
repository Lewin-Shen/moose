//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosCHFreeEnergyCoupledTest.h"

#include "MooseVariable.h"

registerKokkosMaterial("PhaseFieldTestApp", KokkosCHFreeEnergyCoupledTest);

InputParameters
KokkosCHFreeEnergyCoupledTest::validParams()
{
  InputParameters params = Moose::Kokkos::Material::validParams();
  params.addClassDescription("Test-only two-variable free energy f = W (c - ca)^2 (c - cb)^2 + "
                             "a c eta^2 with all first and second derivatives, for the "
                             "cross-derivative Jacobian gates of the Kokkos CH kernels.");
  params.addParam<MaterialPropertyName>("property_name", "f", "Name of the free energy property");
  params.addRequiredCoupledVar("c", "Concentration variable");
  params.addRequiredCoupledVar("eta", "Second variable");
  params.addParam<Real>("W", 16.0, "Well height");
  params.addParam<Real>("ca", 0.27, "First well position");
  params.addParam<Real>("cb", 0.83, "Second well position");
  params.addParam<Real>("a", 1.0, "Coupling coefficient of the a c eta^2 term");
  return params;
}

KokkosCHFreeEnergyCoupledTest::KokkosCHFreeEnergyCoupledTest(const InputParameters & parameters)
  : Moose::Kokkos::Material(parameters),
    _W(getParam<Real>("W")),
    _ca(getParam<Real>("ca")),
    _cb(getParam<Real>("cb")),
    _a(getParam<Real>("a")),
    _c(kokkosCoupledValue("c")),
    _eta(kokkosCoupledValue("eta")),
    _f(declareKokkosProperty<Real>("property_name"))
{
  if (!isCoupled("c") || coupledComponents("c") != 1)
    paramError("c", "Exactly one variable is required.");
  if (!isCoupled("eta") || coupledComponents("eta") != 1)
    paramError("eta", "Exactly one variable is required.");

  const MaterialPropertyName f = getParam<MaterialPropertyName>("property_name");
  const auto cname = coupledName("c", 0);
  const auto ename = coupledName("eta", 0);

  _df_dc = declareKokkosPropertyByName<Real>(derivativePropertyNameFirst(f, cname));
  _d2f_dc2 = declareKokkosPropertyByName<Real>(derivativePropertyNameSecond(f, cname, cname));
  _df_deta = declareKokkosPropertyByName<Real>(derivativePropertyNameFirst(f, ename));
  _d2f_deta2 = declareKokkosPropertyByName<Real>(derivativePropertyNameSecond(f, ename, ename));
  // Sorted denominator: d^2f/dcdeta whichever order is asked for
  _d2f_dcdeta = declareKokkosPropertyByName<Real>(derivativePropertyNameSecond(f, cname, ename));
}
