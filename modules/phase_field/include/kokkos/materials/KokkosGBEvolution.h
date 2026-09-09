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
 * Kokkos fusion of GBEvolution and GBEvolutionBase: the ten fixed-name properties of the
 * isotropic grain-growth model (sigma, M_GB, kappa_op, gamma_asymm, L, l_GB, mu, entropy_diff,
 * molar_volume, act_wGB) for a literal or a prescribed auxiliary temperature. dL/dT is declared
 * only when T is a coupled field (as GBEvolutionBase does), including the constant-mobility
 * branch where its value is zero. Nonlinear T is outside this specialization.
 */
class KokkosGBEvolution : public Moose::Kokkos::Material,
                          public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  KokkosGBEvolution(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION void computeQpProperties(const unsigned int qp, Datum & datum) const;

private:
  const Real _f0s;
  const Real _wGB;
  const Real _length_scale;
  const Real _time_scale;
  const Real _GBmob0;
  const Real _Q;
  const Real _GBMobility;
  const Real _molar_vol;
  const Real _GBEnergy;

  /// Temperature: a coupled auxiliary field or the literal default value
  const Moose::Kokkos::VariableValue _T;
  /// Whether T is a coupled field (then dL/dT is declared and written)
  const bool _has_T;

  Moose::Kokkos::MaterialProperty<Real> _sigma;
  Moose::Kokkos::MaterialProperty<Real> _M_GB;
  Moose::Kokkos::MaterialProperty<Real> _kappa;
  Moose::Kokkos::MaterialProperty<Real> _gamma;
  Moose::Kokkos::MaterialProperty<Real> _L;
  Moose::Kokkos::MaterialProperty<Real> _l_GB;
  Moose::Kokkos::MaterialProperty<Real> _mu;
  Moose::Kokkos::MaterialProperty<Real> _entropy_diff;
  Moose::Kokkos::MaterialProperty<Real> _molar_volume;
  Moose::Kokkos::MaterialProperty<Real> _act_wGB;
  Moose::Kokkos::MaterialProperty<Real> _dLdT;

  /// Boltzmann constant in eV/K and the Joule-to-eV conversion, as in GBEvolutionBase
  const Real _kb = 8.617343e-5;
  const Real _JtoeV = 6.24150974e18;
};

template <typename Derived>
KOKKOS_FUNCTION void
KokkosGBEvolution::computeQpProperties(const unsigned int qp, Datum & datum) const
{
  const Real temperature = _T(datum, qp);
  const Real length_scale4 = _length_scale * _length_scale * _length_scale * _length_scale;

  // GBEvolution::computeQpProperties: eV/nm^2
  _sigma(datum, qp) = _GBEnergy * _JtoeV * (_length_scale * _length_scale);

  // GBEvolutionBase::computeQpProperties
  Real mobility;
  Real dM_GBdT;
  if (_GBMobility < 0)
  {
    // Convert to lengthscale^4/(eV*timescale)
    const Real M0 = _GBmob0 * _time_scale / (_JtoeV * length_scale4);
    mobility = M0 * ::Kokkos::exp(-_Q / (_kb * temperature));
    dM_GBdT = mobility * _Q / (_kb * temperature * temperature);
  }
  else
  {
    mobility = _GBMobility * _time_scale / (_JtoeV * length_scale4);
    dM_GBdT = 0.0;
  }
  _M_GB(datum, qp) = mobility;

  // in the length scale of the system
  _l_GB(datum, qp) = _wGB;

  _L(datum, qp) = 4.0 / 3.0 * mobility / _wGB;
  if (_has_T)
    _dLdT(datum, qp) = 4.0 / 3.0 * dM_GBdT / _wGB;
  _kappa(datum, qp) = 3.0 / 4.0 * _sigma(datum, qp) * _wGB;
  _gamma(datum, qp) = 1.5;
  _mu(datum, qp) = 3.0 / 4.0 * 1.0 / _f0s * _sigma(datum, qp) / _wGB;

  // J/(K mol) converted to eV/(K mol)
  _entropy_diff(datum, qp) = 8.0e3 * _JtoeV;

  // m^3/mol converted to ls^3/mol
  _molar_volume(datum, qp) = _molar_vol / (_length_scale * _length_scale * _length_scale);
  _act_wGB(datum, qp) = 0.5e-9 / _length_scale;
}
