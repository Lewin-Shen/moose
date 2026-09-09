//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosGBEvolution.h"

#include "MooseVariable.h"

registerKokkosMaterial("PhaseFieldApp", KokkosGBEvolution);

InputParameters
KokkosGBEvolution::validParams()
{
  InputParameters params = Moose::Kokkos::Material::validParams();
  params.addClassDescription(
      "Computes necessary material properties for the isotropic grain growth model");
  params.addRequiredCoupledVar("T", "Temperature in Kelvin");
  params.addParam<Real>("f0s", 0.125, "The GB energy constant ");
  params.addRequiredParam<Real>("wGB", "Diffuse GB width in the length scale of the model");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale in m, where default is nm");
  params.addParam<Real>("time_scale", 1.0e-9, "Time scale in s, where default is ns");
  params.addParam<Real>(
      "GBMobility",
      -1,
      "GB mobility input in m^4/(J*s), that overrides the temperature dependent calculation");
  params.addParam<Real>("GBmob0", 0, "Grain boundary mobility prefactor in m^4/(J*s)");
  params.addParam<Real>("Q", 0, "Grain boundary migration activation energy in eV");
  params.addParam<Real>("molar_volume",
                        24.62e-6,
                        "Molar volume in m^3/mol, needed for temperature gradient driving force");
  params.addRequiredParam<Real>("GBenergy", "Grain boundary energy in J/m^2");
  return params;
}

KokkosGBEvolution::KokkosGBEvolution(const InputParameters & parameters)
  : Moose::Kokkos::Material(parameters),
    _f0s(getParam<Real>("f0s")),
    _wGB(getParam<Real>("wGB")),
    _length_scale(getParam<Real>("length_scale")),
    _time_scale(getParam<Real>("time_scale")),
    _GBmob0(getParam<Real>("GBmob0")),
    _Q(getParam<Real>("Q")),
    _GBMobility(getParam<Real>("GBMobility")),
    _molar_vol(getParam<Real>("molar_volume")),
    _GBEnergy(getParam<Real>("GBenergy")),
    _T(kokkosCoupledValue("T")),
    _has_T(!parameters.hasDefaultCoupledValue("T")),
    _sigma(declareKokkosProperty<Real>("sigma")),
    _M_GB(declareKokkosProperty<Real>("M_GB")),
    _kappa(declareKokkosProperty<Real>("kappa_op")),
    _gamma(declareKokkosProperty<Real>("gamma_asymm")),
    _L(declareKokkosProperty<Real>("L")),
    _l_GB(declareKokkosProperty<Real>("l_GB")),
    _mu(declareKokkosProperty<Real>("mu")),
    _entropy_diff(declareKokkosProperty<Real>("entropy_diff")),
    _molar_volume(declareKokkosProperty<Real>("molar_volume")),
    _act_wGB(declareKokkosProperty<Real>("act_wGB"))
{
  // GBEvolutionBase's sentinel check, verbatim
  if (_GBMobility == -1 && _GBmob0 == 0)
    mooseError("Either a value for GBMobility or for GBmob0 and Q must be provided");

  if (coupledComponents("T") != 1)
    paramError("T", "Exactly one prescribed temperature value or auxiliary field is supported.");

  if (_has_T)
  {
    // D-020 scope: a prescribed (auxiliary) field sampled at the quadrature points
    if (getVar("T", 0)->kind() != Moose::VAR_AUXILIARY)
      paramError("T", "Nonlinear temperature is outside the Kokkos grain-growth specialization.");
    if (getParam<MooseEnum>("constant_on") != "NONE")
      paramError("constant_on", "Use NONE for a quadrature-point auxiliary temperature.");
    _dLdT = declareKokkosProperty<Real>(derivativePropertyNameFirst("L", coupledName("T", 0)));
  }
}
