//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosKernelGrad.h"
#include "KokkosMap.h"
#include "DerivativeMaterialPropertyNameInterface.h"

/**
 * Kokkos specialization of SplitCHWRes (SplitCHWResBase<Real>), the chemical-potential row of the
 * split Cahn-Hilliard equation with a scalar (isotropic) mobility M:
 *
 *   R_i = M grad(w) . grad(psi_i)
 *
 * where w is the kernel variable unless 'w' couples another variable. Hand-coded Jacobians
 * transcribed from SplitCHWResBase.h: the diagonal is M grad(phi_j) . grad(psi_i) when w is the
 * kernel variable and zero when 'w' couples a different variable; the off-diagonal against a coupled
 * w is M grad(phi_j) . grad(psi_i); the off-diagonal against every nonlinear argument of the
 * mobility listed in coupled_variables is dM/darg phi_j grad(w) . grad(psi_i). A derivative
 * property the mobility material does not declare contributes zero, which is what the host's
 * DerivativeMaterialInterface zero default does. The test-function gradient is factored out by
 * KernelGrad, so every hook returns the vector it multiplies.
 */
class KokkosSplitCHWRes : public Moose::Kokkos::KernelGrad,
                          public DerivativeMaterialPropertyNameInterface
{
  using Real3 = Moose::Kokkos::Real3;

public:
  static InputParameters validParams();

  KokkosSplitCHWRes(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION Real3 precomputeQpResidual(const unsigned int qp, AssemblyDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION Real3 precomputeQpJacobian(const unsigned int j,
                                             const unsigned int qp,
                                             AssemblyDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION Real3 precomputeQpOffDiagJacobian(const unsigned int j,
                                                    const unsigned int jvar,
                                                    const unsigned int qp,
                                                    AssemblyDatum & datum) const;

private:
  /// SplitCHWResBase::computeQpWJacobian without the test-function gradient: M grad(phi_j)
  KOKKOS_FUNCTION Real3 wJacobian(const unsigned int j,
                                  const unsigned int qp,
                                  AssemblyDatum & datum) const
  {
    const Real mob = _mob(datum, qp);
    return mob * _grad_phi(datum, j, qp);
  }

  /// Mobility M
  const Moose::Kokkos::MaterialProperty<Real> _mob;
  /// Whether 'w' couples a variable (SplitCHWResBase::_is_coupled)
  const bool _is_coupled;
  /// Variable number of the chemical potential (the kernel variable when 'w' is not coupled)
  const unsigned int _w_var;
  /// Whether the chemical potential is the kernel variable, so the diagonal carries M grad(phi_j)
  const bool _w_is_variable;
  /// Gradient of the chemical potential (the kernel variable's gradient when 'w' is not coupled)
  const Moose::Kokkos::VariableGradient _grad_w;
  /// Number of coupled_variables entries (the arguments of the mobility)
  const unsigned int _n_args;
  /// Host-only metadata: variable numbers of the arguments, in coupled_variables order
  std::vector<unsigned int> _arg_var;
  /// Variable number -> argument index, for the nonlinear arguments of this kernel's system
  Moose::Kokkos::Map<unsigned int, unsigned int> _arg_var_to_index;
  /// dM/darg handles indexed by argument; an unset handle means no declared derivative (zero)
  Moose::Kokkos::Array<Moose::Kokkos::MaterialProperty<Real>> _dmob_darg;
};

template <typename Derived>
KOKKOS_FUNCTION Moose::Kokkos::Real3
KokkosSplitCHWRes::precomputeQpResidual(const unsigned int qp, AssemblyDatum & datum) const
{
  // SplitCHWResBase::computeQpResidual: M grad(w) . grad(psi_i)
  const Real mob = _mob(datum, qp);
  return mob * _grad_w(datum, qp);
}

template <typename Derived>
KOKKOS_FUNCTION Moose::Kokkos::Real3
KokkosSplitCHWRes::precomputeQpJacobian(const unsigned int j,
                                        const unsigned int qp,
                                        AssemblyDatum & datum) const
{
  // SplitCHWResBase::computeQpJacobian: zero when 'w' couples a variable other than the kernel
  // variable, computeQpWJacobian otherwise
  if (!_w_is_variable)
    return Real3(0);
  return wJacobian(j, qp, datum);
}

template <typename Derived>
KOKKOS_FUNCTION Moose::Kokkos::Real3
KokkosSplitCHWRes::precomputeQpOffDiagJacobian(const unsigned int j,
                                               const unsigned int jvar,
                                               const unsigned int qp,
                                               AssemblyDatum & datum) const
{
  // SplitCHWResBase::computeQpOffDiagJacobian: the w column (reached only when 'w' couples a
  // variable other than the kernel variable) carries computeQpWJacobian
  if (jvar == _w_var)
    return wJacobian(j, qp, datum);

  // dM/darg phi_j grad(w) for the nonlinear arguments of the mobility; every other column is zero
  // (the host skips them in JvarMapKernelInterface::computeOffDiagJacobian)
  if (!_n_args || !_arg_var_to_index.exists(jvar))
    return Real3(0);

  const auto & dmob = _dmob_darg[_arg_var_to_index[jvar]];
  if (!dmob)
    return Real3(0);

  const Real d = dmob(datum, qp);
  return (d * _phi(datum, j, qp)) * _grad_w(datum, qp);
}
