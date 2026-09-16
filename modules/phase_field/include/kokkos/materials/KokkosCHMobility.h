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
 * Cahn-Hilliard mobility M(c) with its derivative dM/dc, in two closed forms selected by 'form':
 *
 *   CONSTANT:    M = M0,            dM/dc = 0
 *   DEGENERATE:  M = M0 c (1 - c),  dM/dc = M0 (1 - 2c)
 *
 * Declared properties: <property_name> and d<property_name>/d<c> (both forms declare both, so a deck
 * switches between them by the 'form' parameter alone). The derivative name is built by
 * DerivativeMaterialPropertyNameInterface on the host, so KokkosMatDiffusion and KokkosSplitCHWRes
 * pick it up under exactly the name the CPU kernels request from a DerivativeParsedMaterial.
 */
class KokkosCHMobility : public Moose::Kokkos::Material,
                         public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  KokkosCHMobility(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION void computeQpProperties(const unsigned int qp, Datum & datum) const;

private:
  enum Form : unsigned int
  {
    CONSTANT = 0,
    DEGENERATE = 1
  };

  /// Selected closed form (the MooseEnum resolved on the host)
  const unsigned int _form;
  /// Mobility prefactor
  const Real _M0;
  /// Concentration
  const Moose::Kokkos::VariableValue _c;
  /// M and dM/dc
  Moose::Kokkos::MaterialProperty<Real> _M;
  Moose::Kokkos::MaterialProperty<Real> _dMdc;
};

template <typename Derived>
KOKKOS_FUNCTION void
KokkosCHMobility::computeQpProperties(const unsigned int qp, Datum & datum) const
{
  if (_form == CONSTANT)
  {
    _M(datum, qp) = _M0;
    _dMdc(datum, qp) = 0.0;
  }
  else
  {
    const Real c = _c(datum, qp);
    _M(datum, qp) = _M0 * c * (1.0 - c);
    _dMdc(datum, qp) = _M0 * (1.0 - 2.0 * c);
  }
}
