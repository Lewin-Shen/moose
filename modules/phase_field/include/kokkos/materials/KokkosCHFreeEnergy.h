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
#include "DerivativeMaterialPropertyNameInterface.h"

/**
 * Kokkos replacement for the DerivativeParsedMaterial that supplies the Cahn-Hilliard local free
 * energy f(c) together with its first two derivatives. Two closed forms, selected by 'form':
 *
 *   DOUBLE_WELL:     f = W (c - ca)^2 (c - cb)^2
 *   REDLICH_KISTER:  f = p [A c + B (1 - c) + C c ln(c) + D (1 - c) ln(1 - c) + E c (1 - c)
 *                          + F c (1 - c) (2c - 1) + G c (1 - c) (2c - 1)^2],   p = eVpJ / Vm * length_scale^3
 *
 * This is NOT a general expression material: the two energies and their derivatives are
 * hand-written closed forms compiled ahead of time (the Kokkos function parser evaluates values
 * only and cannot differentiate). Every coefficient is an input parameter whose default reproduces
 * the tutorial decks (spinodal_polynomial.i for DOUBLE_WELL, spinodal.i for REDLICH_KISTER).
 *
 * The logarithms are evaluated as written, without clamping, so that the material matches the CPU
 * DerivativeParsedMaterial exactly; outside (0, 1) the REDLICH_KISTER form returns NaN on both
 * backends.
 *
 * Declared properties: <property_name>, d<property_name>/d<c>, d^2<property_name>/d<c>^2, with the
 * names built by DerivativeMaterialPropertyNameInterface on the host, so they are exactly the names
 * the CPU objects (DerivativeParsedMaterial, MatDiffusion, CoupledMaterialDerivative, SplitCHParsed)
 * declare and request.
 */
class KokkosCHFreeEnergy : public Moose::Kokkos::Material,
                           public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  KokkosCHFreeEnergy(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION void computeQpProperties(const unsigned int qp, Datum & datum) const;

private:
  enum Form : unsigned int
  {
    DOUBLE_WELL = 0,
    REDLICH_KISTER = 1
  };

  /// Selected closed form (the MooseEnum resolved on the host)
  const unsigned int _form;

  /// DOUBLE_WELL: well height and the two well positions
  const Real _W;
  const Real _ca;
  const Real _cb;

  /// REDLICH_KISTER: the seven coefficients and the prefactor p = eVpJ / Vm * length_scale^3
  const Real _A;
  const Real _B;
  const Real _C;
  const Real _D;
  const Real _E;
  const Real _F;
  const Real _G;
  const Real _p;

  /// Concentration
  const Moose::Kokkos::VariableValue _c;

  /// f, df/dc, d^2f/dc^2
  Moose::Kokkos::MaterialProperty<Real> _f;
  Moose::Kokkos::MaterialProperty<Real> _df;
  Moose::Kokkos::MaterialProperty<Real> _d2f;
};

template <typename Derived>
KOKKOS_FUNCTION void
KokkosCHFreeEnergy::computeQpProperties(const unsigned int qp, Datum & datum) const
{
  const Real c = _c(datum, qp);

  if (_form == DOUBLE_WELL)
  {
    // f = W p^2 q^2 with p = c - ca, q = c - cb; f' = 2 W p q (p + q); f'' = 2 W (p^2 + q^2 + 4 p q)
    const Real p = c - _ca;
    const Real q = c - _cb;
    _f(datum, qp) = _W * p * p * q * q;
    _df(datum, qp) = 2.0 * _W * p * q * (p + q);
    _d2f(datum, qp) = 2.0 * _W * (p * p + q * q + 4.0 * p * q);
  }
  else
  {
    // The tutorial expression as written (no clamp), with
    //   c (1 - c) (2c - 1)   = -2 c^3 + 3 c^2 - c
    //   c (1 - c) (2c - 1)^2 = -4 c^4 + 8 c^3 - 5 c^2 + c
    // differentiated term by term.
    const Real omc = 1.0 - c;
    const Real lnc = ::Kokkos::log(c);
    const Real lnomc = ::Kokkos::log(omc);
    const Real c2 = c * c;
    const Real c3 = c2 * c;
    const Real twocm1 = 2.0 * c - 1.0;

    _f(datum, qp) = _p * (_A * c + _B * omc + _C * c * lnc + _D * omc * lnomc + _E * c * omc +
                          _F * c * omc * twocm1 + _G * c * omc * twocm1 * twocm1);
    _df(datum, qp) = _p * ((_A - _B) + _C * (lnc + 1.0) - _D * (lnomc + 1.0) +
                           _E * (1.0 - 2.0 * c) + _F * (-6.0 * c2 + 6.0 * c - 1.0) +
                           _G * (-16.0 * c3 + 24.0 * c2 - 10.0 * c + 1.0));
    _d2f(datum, qp) = _p * (_C / c + _D / omc - 2.0 * _E + _F * (6.0 - 12.0 * c) +
                            _G * (-48.0 * c2 + 48.0 * c - 10.0));
  }
}
