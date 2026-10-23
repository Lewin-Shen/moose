//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosTwoPhaseParabolicFreeEnergy.h"

#include "MooseVariable.h"

#include <set>

registerKokkosMaterial("PhaseFieldApp", KokkosTwoPhaseParabolicFreeEnergy);

InputParameters
KokkosTwoPhaseParabolicFreeEnergy::validParams()
{
  InputParameters params = Moose::Kokkos::Material::validParams();
  params.addClassDescription(
      "Two-phase parabolic free energy of a concentration c and n order parameters with its first "
      "and second derivatives in closed form: f = rho_sq (c - ca)^2 (1 - h) + rho_sq (cb - c)^2 h + "
      "w g, the two parabolic phase energies interpolated by the sum of the smoothsteps "
      "h = sum_k eta_k^3 (6 eta_k^2 - 15 eta_k + 10), plus the multi-well with the pairwise penalty "
      "g = sum_k eta_k^2 (1 - eta_k)^2 + 2 alpha sum_{k<m} eta_k^2 eta_m^2, where w is the barrier "
      "height. Declares f and all its first and second derivatives under the "
      "DerivativeParsedMaterial property names; the pair terms are always declared.");
  params.addParam<MaterialPropertyName>(
      "property_name", "f_loc", "Name of the free energy property (derivative names derive from it)");
  params.addRequiredCoupledVar("c", "Concentration variable");
  params.addRequiredCoupledVar("etas", "Order parameter variables (one or more)");
  params.addParam<Real>("rho_sq", 2.0, "Curvature of the two parabolic phase energies");
  params.addParam<Real>("ca", 0.3, "Composition at the minimum of the first phase energy");
  params.addParam<Real>("cb", 0.7, "Composition at the minimum of the second phase energy");
  params.addParam<Real>("w", 1.0, "Barrier height of the multi-well (not a chemical potential)");
  params.addParam<Real>("alpha", 5.0, "Pairwise order parameter penalty coefficient");
  return params;
}

KokkosTwoPhaseParabolicFreeEnergy::KokkosTwoPhaseParabolicFreeEnergy(
    const InputParameters & parameters)
  : Moose::Kokkos::Material(parameters),
    _n(coupledComponents("etas")),
    _rho_sq(getParam<Real>("rho_sq")),
    _ca(getParam<Real>("ca")),
    _cb(getParam<Real>("cb")),
    _w(getParam<Real>("w")),
    _alpha(getParam<Real>("alpha")),
    _f(declareKokkosProperty<Real>("property_name"))
{
  // The derivative property names are built from the variable names, so 'c' and every entry of
  // 'etas' must be a variable (auxiliary or nonlinear), not a literal default value
  if (!isCoupled("c") || coupledComponents("c") != 1)
    paramError("c", "Exactly one variable is required (a literal value is not supported).");
  if (!_n)
    paramError("etas", "At least one order parameter variable is required.");

  const auto cname = coupledName("c", 0);
  std::vector<VariableName> eta_names(_n);
  std::set<VariableName> seen;
  for (unsigned int k = 0; k < _n; ++k)
  {
    if (!isCoupled("etas", k))
      paramError("etas", "Entries must be variables, not literal values.");
    eta_names[k] = coupledName("etas", k);
    if (eta_names[k] == cname)
      paramError("etas", "The concentration variable 'c' must not appear in 'etas'.");
    if (!seen.insert(eta_names[k]).second)
      paramError("etas", "Repeated order parameters are not supported.");
  }

  _c = kokkosCoupledValue("c");
  _etas = kokkosCoupledValues("etas");

  // Every derivative name in one block, through DerivativeMaterialPropertyNameInterface (the same
  // helper DerivativeParsedMaterial uses, so the strings are the CPU names, sorted). Every array is
  // created before it is filled, the pair array with n (n - 1) / 2 entries even when that is zero:
  // the dispatcher deep-copies these arrays and aborts on one that was never allocated.
  const MaterialPropertyName f = getParam<MaterialPropertyName>("property_name");
  _dfdc = declareKokkosPropertyByName<Real>(derivativePropertyNameFirst(f, cname));
  _d2fdc2 = declareKokkosPropertyByName<Real>(derivativePropertyNameSecond(f, cname, cname));
  _dfdeta.create(_n);
  _d2fdcdeta.create(_n);
  _d2fdeta2.create(_n);
  _d2fdetadeta.create(_n * (_n - 1) / 2);
  for (unsigned int k = 0; k < _n; ++k)
  {
    _dfdeta[k] = declareKokkosPropertyByName<Real>(derivativePropertyNameFirst(f, eta_names[k]));
    _d2fdcdeta[k] =
        declareKokkosPropertyByName<Real>(derivativePropertyNameSecond(f, cname, eta_names[k]));
    _d2fdeta2[k] = declareKokkosPropertyByName<Real>(
        derivativePropertyNameSecond(f, eta_names[k], eta_names[k]));
    for (unsigned int m = k + 1; m < _n; ++m)
      _d2fdetadeta[pairIndex(k, m)] = declareKokkosPropertyByName<Real>(
          derivativePropertyNameSecond(f, eta_names[k], eta_names[m]));
  }
}
