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
#include "DerivativeMaterialPropertyNameInterface.h"

/**
 * Kokkos specialization of MatDiffusion (MatDiffusionBase<Real>): the isotropic diffusion term
 *
 *   R_i = D grad(v) . grad(psi_i)
 *
 * with D a material property that may depend on the kernel variable u and on the nonlinear
 * arguments listed in 'args' (MatDiffusionBase's parameter name), and v either the coupled
 * variable 'v' or, when 'v' is not given, the kernel variable itself. The Jacobian is
 * MatDiffusionBase's, transcribed exactly:
 *
 *   diagonal     : dD/du phi_j grad(v)   + [ D grad(phi_j) only when v is NOT coupled ]
 *   off-diagonal : dD/darg phi_j grad(v) + [ D grad(phi_j) when jvar is v ]
 *
 * both dotted with grad(psi_i) by the KernelGrad base. So with 'v' set the Laplacian lives in the
 * (u, v) off-diagonal block only, and with 'v' unset it lives on the diagonal only; a twin that
 * always or never adds it double-counts or drops the Laplacian.
 *
 * The derivative properties dD/du, dD/dv and dD/darg are optional: the CPU kernel gets zero
 * defaults from DerivativeMaterialInterface, the Kokkos side has no such default, so each handle is
 * fetched only when a material declares it and reads as zero otherwise. The names come from
 * DerivativeMaterialPropertyNameInterface on the host, so a material declaring dM/dc is picked up
 * under exactly the name the CPU kernel requests. Auxiliary entries of 'args' are accepted but
 * never matched against jvar (they have no nonlinear column).
 */
class KokkosMatDiffusion : public Moose::Kokkos::KernelGrad,
                           public DerivativeMaterialPropertyNameInterface
{
  using Real3 = Moose::Kokkos::Real3;

public:
  static InputParameters validParams();

  KokkosMatDiffusion(const InputParameters & parameters);

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
  /// MatDiffusionBase::precomputeQpCJacobian: D grad(phi_j)
  KOKKOS_FUNCTION Real3 cJacobian(const unsigned int j,
                                  const unsigned int qp,
                                  AssemblyDatum & datum) const
  {
    const Real D = _diffusivity(datum, qp);
    return D * _grad_phi(datum, j, qp);
  }

  /// Diffusivity D
  const Moose::Kokkos::MaterialProperty<Real> _diffusivity;
  /// Whether the kernel operates on a coupled variable v (MatDiffusionBase::_is_coupled)
  const bool _is_coupled;
  /// Variable number of v, or of the kernel variable when v is not coupled (MatDiffusionBase::_v_var)
  const unsigned int _v_var;
  /// grad(v), or grad(u) when v is not coupled
  const Moose::Kokkos::VariableGradient _grad_v;
  /// Number of 'args' entries
  const unsigned int _n_args;
  /// dD/du (unset, reads as zero, when no material declares it)
  Moose::Kokkos::MaterialProperty<Real> _dD_du;
  /// dD/dv (unset when no material declares it, or when v is not coupled)
  Moose::Kokkos::MaterialProperty<Real> _dD_dv;
  /// Nonlinear variable numbers of the 'args' entries (invalid_uint for auxiliary ones)
  Moose::Kokkos::Array<unsigned int> _arg_var;
  /// dD/darg per 'args' entry (each unset when not declared); always created, possibly empty,
  /// because the dispatcher deep-copies this array and aborts on an unallocated one
  Moose::Kokkos::Array<Moose::Kokkos::MaterialProperty<Real>> _dD_darg;
};

template <typename Derived>
KOKKOS_FUNCTION Moose::Kokkos::Real3
KokkosMatDiffusion::precomputeQpResidual(const unsigned int qp, AssemblyDatum & datum) const
{
  // MatDiffusionBaseTempl::precomputeQpResidual
  const Real D = _diffusivity(datum, qp);
  return D * _grad_v(datum, qp);
}

template <typename Derived>
KOKKOS_FUNCTION Moose::Kokkos::Real3
KokkosMatDiffusion::precomputeQpJacobian(const unsigned int j,
                                         const unsigned int qp,
                                         AssemblyDatum & datum) const
{
  // MatDiffusionBase::precomputeQpJacobian
  Real3 sum(0);
  if (_dD_du)
  {
    const Real dDdu = _dD_du(datum, qp);
    sum += (dDdu * _phi(datum, j, qp)) * _grad_v(datum, qp);
  }
  // The Laplacian enters the diagonal only when v is not coupled (MatDiffusionBase.C:86-87)
  if (!_is_coupled)
    sum += cJacobian(j, qp, datum);
  return sum;
}

template <typename Derived>
KOKKOS_FUNCTION Moose::Kokkos::Real3
KokkosMatDiffusion::precomputeQpOffDiagJacobian(const unsigned int j,
                                                const unsigned int jvar,
                                                const unsigned int qp,
                                                AssemblyDatum & datum) const
{
  // MatDiffusionBase::computeQpOffDiagJacobian, divided by grad(psi_i) (KernelGrad applies it).
  // jvar is compared against fetched variable numbers, never positions.
  Real3 sum(0);
  if (jvar == _v_var)
  {
    if (_dD_dv)
    {
      const Real dDdv = _dD_dv(datum, qp);
      sum += (dDdv * _phi(datum, j, qp)) * _grad_v(datum, qp);
    }
    sum += cJacobian(j, qp, datum);
    return sum;
  }
  for (unsigned int k = 0; k < _n_args; ++k)
    if (jvar == _arg_var[k])
    {
      if (_dD_darg[k])
      {
        const Real dDdarg = _dD_darg[k](datum, qp);
        sum += (dDdarg * _phi(datum, j, qp)) * _grad_v(datum, qp);
      }
      return sum;
    }
  // Not coupled to jvar: JvarMapKernelInterface skips the block on the CPU
  return sum;
}
