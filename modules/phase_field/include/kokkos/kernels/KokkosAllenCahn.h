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
 * Kokkos specialization of AllenCahn (the Allen-Cahn bulk term of ACBulk<Real>) for a constant
 * mobility L:
 *
 *   R_i = L dF/du psi_i
 *
 * with the Jacobian of AllenCahn transcribed exactly for constant L (KernelValue applies psi_i;
 * phi_j in the off-diagonal hook is the column variable's shape function, as the Kokkos base supplies
 * it, where the host uses the kernel variable's: the same function while every variable has the same
 * finite element family and order, which the gate decks and the benchmark deck have):
 *
 *   diagonal     : L d^2F/du^2 phi_j
 *   off-diagonal : L d^2F/du darg phi_j for an argument of 'coupled_variables'
 *
 * Constant mobility: the ACBulk product-rule terms dL/du phi_j dF/du (diagonal) and
 * dL/darg phi_j dF/du (off-diagonal) are omitted by construction; a declared dL/du or dL/darg is
 * rejected at construction, so the omission can never drop a term a material provides.
 *
 * dF/du must be declared by a Kokkos material (a missing derivative is a zero property on the CPU,
 * which makes the kernel a no-op, so it is an error here). The second derivatives are optional:
 * fetched only when declared, zero otherwise. All names come from
 * DerivativeMaterialPropertyNameInterface on the host (sorted denominators), so they match what a
 * DerivativeParsedMaterial or KokkosTwoPhaseParabolicFreeEnergy declares. Auxiliary entries of
 * 'coupled_variables' are accepted but never matched against jvar.
 */
class KokkosAllenCahn : public Moose::Kokkos::KernelValue,
                        public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  KokkosAllenCahn(const InputParameters & parameters);

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
  /// Number of coupled_variables entries
  const unsigned int _n_args;
  /// Mobility L (constant: no declared derivative is accepted)
  const Moose::Kokkos::MaterialProperty<Real> _L;
  /// dF/du (required)
  Moose::Kokkos::MaterialProperty<Real> _dFdu;
  /// d^2F/du^2 (unset, reads as zero, when no material declares it)
  Moose::Kokkos::MaterialProperty<Real> _d2Fdu2;
  /// Nonlinear variable numbers of the coupled_variables entries (invalid_uint for auxiliary ones)
  Moose::Kokkos::Array<unsigned int> _arg_var;
  /// d^2F/du darg per coupled_variables entry (each unset when not declared); always created,
  /// possibly empty, because the dispatcher deep-copies this array and aborts on an unallocated one
  Moose::Kokkos::Array<Moose::Kokkos::MaterialProperty<Real>> _d2Fdudarg;
};

template <typename Derived>
KOKKOS_FUNCTION Real
KokkosAllenCahn::precomputeQpResidual(const unsigned int qp, AssemblyDatum & datum) const
{
  // ACBulk::precomputeQpResidual = L * dF/dop, with AllenCahn::computeDFDOP(Residual) = dF/du
  const Real dFdop = _dFdu(datum, qp);
  return _L(datum, qp) * dFdop;
}

template <typename Derived>
KOKKOS_FUNCTION Real
KokkosAllenCahn::precomputeQpJacobian(const unsigned int j,
                                      const unsigned int qp,
                                      AssemblyDatum & datum) const
{
  // ACBulk::precomputeQpJacobian = L * JdFdop + dL/dop phi_j dFdop, with
  // AllenCahn::computeDFDOP(Jacobian) = d^2F/du^2 phi_j; the dL/dop term is zero for every material
  // this specialization accepts (see the constructor guards)
  if (!_d2Fdu2)
    return 0.0;
  const Real JdFdop = _d2Fdu2(datum, qp) * _phi(datum, j, qp);
  return _L(datum, qp) * JdFdop;
}

template <typename Derived>
KOKKOS_FUNCTION Real
KokkosAllenCahn::precomputeQpOffDiagJacobian(const unsigned int j,
                                             const unsigned int jvar,
                                             const unsigned int qp,
                                             AssemblyDatum & datum) const
{
  // AllenCahn::computeQpOffDiagJacobian divided by the test function (KernelValue applies it):
  // L (*_d2FdEtadarg[cvar]) phi_j, plus ACBulk's dL/darg phi_j dF/du, zero here by the guards.
  // jvar is compared against fetched variable numbers, never positions.
  for (unsigned int k = 0; k < _n_args; ++k)
    if (jvar == _arg_var[k])
    {
      if (!_d2Fdudarg[k])
        return 0.0;
      const Real d2Fdudarg = _d2Fdudarg[k](datum, qp);
      return _L(datum, qp) * d2Fdudarg * _phi(datum, j, qp);
    }
  // Not coupled to jvar: JvarMapKernelInterface skips the block on the CPU
  return 0.0;
}
