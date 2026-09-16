//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosKernelValue.h"
#include "DerivativeMaterialPropertyNameInterface.h"

/**
 * Kokkos specialization of CoupledMaterialDerivative: the first derivative of a function material
 * property F with respect to a coupled variable v, tested against psi_i:
 *
 *   R_i = dF/dv psi_i
 *
 * with the Jacobian of CoupledMaterialDerivative transcribed exactly (KernelValue applies psi_i):
 *
 *   diagonal     : d^2F/dv du phi_j
 *   off-diagonal : d^2F/dv dv phi_j when jvar is v; d^2F/dv darg phi_j for an argument of
 *                  'coupled_variables'
 *
 * dF/dv must be declared by a Kokkos material (a missing derivative is a zero property on the CPU,
 * which makes the kernel a no-op, so it is an error here). The second derivatives are optional:
 * fetched only when declared, zero otherwise. All names come from
 * DerivativeMaterialPropertyNameInterface on the host (sorted denominators, d^2F/dcdeta), so they
 * match what a DerivativeParsedMaterial or KokkosCHFreeEnergy declares. Auxiliary entries of
 * 'coupled_variables' are accepted but never matched against jvar.
 */
class KokkosCoupledMaterialDerivative : public Moose::Kokkos::KernelValue,
                                        public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  KokkosCoupledMaterialDerivative(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION Real precomputeQpResidual(const unsigned int qp, AssemblyDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION Real precomputeQpJacobian(const unsigned int j,
                                            const unsigned int qp,
                                            AssemblyDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION Real precomputeQpOffDiagJacobian(const unsigned int j,
                                                   const unsigned int jvar,
                                                   const unsigned int qp,
                                                   AssemblyDatum & datum) const;

private:
  /// Variable number of v
  const unsigned int _v_var;
  /// Number of coupled_variables entries
  const unsigned int _n_args;
  /// dF/dv (required)
  Moose::Kokkos::MaterialProperty<Real> _dFdv;
  /// d^2F/dv du (unset, reads as zero, when no material declares it)
  Moose::Kokkos::MaterialProperty<Real> _d2Fdvdu;
  /// d^2F/dv dv (unset when no material declares it)
  Moose::Kokkos::MaterialProperty<Real> _d2Fdv2;
  /// Nonlinear variable numbers of the coupled_variables entries (invalid_uint for auxiliary ones)
  Moose::Kokkos::Array<unsigned int> _arg_var;
  /// d^2F/dv darg per coupled_variables entry (each unset when not declared); always created,
  /// possibly empty, because the dispatcher deep-copies this array and aborts on an unallocated one
  Moose::Kokkos::Array<Moose::Kokkos::MaterialProperty<Real>> _d2Fdvdarg;
};

template <typename Derived>
KOKKOS_FUNCTION Real
KokkosCoupledMaterialDerivative::precomputeQpResidual(const unsigned int qp,
                                                      AssemblyDatum & datum) const
{
  // CoupledMaterialDerivative::computeQpResidual without the test function
  const Real dFdv = _dFdv(datum, qp);
  return dFdv;
}

template <typename Derived>
KOKKOS_FUNCTION Real
KokkosCoupledMaterialDerivative::precomputeQpJacobian(const unsigned int j,
                                                      const unsigned int qp,
                                                      AssemblyDatum & datum) const
{
  // CoupledMaterialDerivative::computeQpJacobian: d^2F/dv du phi_j (psi_i applied by the base)
  if (!_d2Fdvdu)
    return 0.0;
  const Real d2Fdvdu = _d2Fdvdu(datum, qp);
  return d2Fdvdu * _phi(datum, j, qp);
}

template <typename Derived>
KOKKOS_FUNCTION Real
KokkosCoupledMaterialDerivative::precomputeQpOffDiagJacobian(const unsigned int j,
                                                             const unsigned int jvar,
                                                             const unsigned int qp,
                                                             AssemblyDatum & datum) const
{
  // CoupledMaterialDerivative::computeQpOffDiagJacobian: (*_d2Fdvdarg[cvar]) phi_j, where cvar runs
  // over the coupled variables including v itself. jvar is compared against fetched variable
  // numbers, never positions.
  if (jvar == _v_var)
  {
    if (!_d2Fdv2)
      return 0.0;
    const Real d2Fdv2 = _d2Fdv2(datum, qp);
    return d2Fdv2 * _phi(datum, j, qp);
  }
  for (unsigned int k = 0; k < _n_args; ++k)
    if (jvar == _arg_var[k])
    {
      if (!_d2Fdvdarg[k])
        return 0.0;
      const Real d2Fdvdarg = _d2Fdvdarg[k](datum, qp);
      return d2Fdvdarg * _phi(datum, j, qp);
    }
  // Not coupled to jvar: JvarMapKernelInterface skips the block on the CPU
  return 0.0;
}
