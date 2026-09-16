//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosKernel.h"
#include "KokkosMap.h"
#include "DerivativeMaterialPropertyNameInterface.h"

/**
 * Kokkos fusion of SplitCHParsed, SplitCHCRes and SplitCHBase: the concentration row of the split
 * Cahn-Hilliard equation with a free energy F(c, args...) supplied by a material that declares
 * dF/dc, d^2F/dc^2 and the cross derivatives d^2F/dc darg under the names of
 * DerivativeMaterialPropertyNameInterface:
 *
 *   R_i = (dF/dc - w) psi_i + kappa grad(c) . grad(psi_i)
 *
 * Hand-coded Jacobians transcribed from SplitCHCRes.C and SplitCHParsed.C: the diagonal is
 * d^2F/dc^2 phi_j psi_i + kappa grad(phi_j) . grad(psi_i); the off-diagonal against w is
 * -phi_j psi_i; the off-diagonal against each nonlinear argument in coupled_variables is
 * d^2F/dc darg phi_j psi_i, zero when the material declares no such derivative (the host's zero
 * default). The residual mixes a value term and a gradient term, so neither precompute base
 * applies and the test function is applied here. SplitCHBase's computeDEDC hook is zero in every
 * class of the host chain and is not carried.
 */
class KokkosSplitCHParsed : public Moose::Kokkos::Kernel,
                            public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  KokkosSplitCHParsed(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION Real computeQpResidual(const unsigned int i,
                                         const unsigned int qp,
                                         AssemblyDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION Real computeQpJacobian(const unsigned int i,
                                         const unsigned int j,
                                         const unsigned int qp,
                                         AssemblyDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION Real computeQpOffDiagJacobian(const unsigned int i,
                                                const unsigned int j,
                                                const unsigned int jvar,
                                                const unsigned int qp,
                                                AssemblyDatum & datum) const;

private:
  /// Gradient energy coefficient kappa
  const Moose::Kokkos::MaterialProperty<Real> _kappa;
  /// Variable number of the chemical potential w
  const unsigned int _w_var;
  /// Chemical potential value
  const Moose::Kokkos::VariableValue _w;
  /// dF/dc
  const Moose::Kokkos::MaterialProperty<Real> _dFdc;
  /// d^2F/dc^2
  const Moose::Kokkos::MaterialProperty<Real> _d2Fdc2;
  /// Number of coupled_variables entries (the additional arguments of F)
  const unsigned int _n_args;
  /// Host-only metadata: variable numbers of the arguments, in coupled_variables order
  std::vector<unsigned int> _arg_var;
  /// Variable number -> argument index, for the nonlinear arguments of this kernel's system
  Moose::Kokkos::Map<unsigned int, unsigned int> _arg_var_to_index;
  /// d^2F/dc darg handles indexed by argument; an unset handle means no declared derivative (zero)
  Moose::Kokkos::Array<Moose::Kokkos::MaterialProperty<Real>> _d2Fdcdarg;
};

template <typename Derived>
KOKKOS_FUNCTION Real
KokkosSplitCHParsed::computeQpResidual(const unsigned int i,
                                       const unsigned int qp,
                                       AssemblyDatum & datum) const
{
  // SplitCHBase::computeQpResidual with SplitCHParsed::computeDFDC(Residual) = dF/dc and
  // computeDEDC = 0, plus SplitCHCRes's -w psi_i and kappa grad(c) . grad(psi_i)
  const Real dFdc = _dFdc(datum, qp);
  const Real kappa = _kappa(datum, qp);
  return (dFdc - _w(datum, qp)) * _test(datum, i, qp) +
         kappa * (_grad_u(datum, qp) * _grad_test(datum, i, qp));
}

template <typename Derived>
KOKKOS_FUNCTION Real
KokkosSplitCHParsed::computeQpJacobian(const unsigned int i,
                                       const unsigned int j,
                                       const unsigned int qp,
                                       AssemblyDatum & datum) const
{
  // SplitCHBase::computeQpJacobian with SplitCHParsed::computeDFDC(Jacobian) = d^2F/dc^2 phi_j,
  // plus SplitCHCRes's kappa grad(phi_j) . grad(psi_i)
  const Real d2Fdc2 = _d2Fdc2(datum, qp);
  const Real kappa = _kappa(datum, qp);
  return d2Fdc2 * _phi(datum, j, qp) * _test(datum, i, qp) +
         kappa * (_grad_phi(datum, j, qp) * _grad_test(datum, i, qp));
}

template <typename Derived>
KOKKOS_FUNCTION Real
KokkosSplitCHParsed::computeQpOffDiagJacobian(const unsigned int i,
                                              const unsigned int j,
                                              const unsigned int jvar,
                                              const unsigned int qp,
                                              AssemblyDatum & datum) const
{
  // SplitCHCRes::computeQpOffDiagJacobian: the w column
  if (jvar == _w_var)
    return -_phi(datum, j, qp) * _test(datum, i, qp);

  // SplitCHParsed::computeQpOffDiagJacobian: d^2F/dc darg phi_j psi_i for the nonlinear
  // arguments of F; every other column is zero (the host skips them in
  // JvarMapKernelInterface::computeOffDiagJacobian)
  if (!_n_args || !_arg_var_to_index.exists(jvar))
    return 0.0;

  const auto & d2F = _d2Fdcdarg[_arg_var_to_index[jvar]];
  if (!d2F)
    return 0.0;

  const Real d = d2F(datum, qp);
  return d * _phi(datum, j, qp) * _test(datum, i, qp);
}
