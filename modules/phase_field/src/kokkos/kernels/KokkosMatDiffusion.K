//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosMatDiffusion.h"

#include "MooseVariable.h"
#include "SystemBase.h"

#include <algorithm>

registerKokkosResidualObject("PhaseFieldApp", KokkosMatDiffusion);

InputParameters
KokkosMatDiffusion::validParams()
{
  InputParameters params = Moose::Kokkos::KernelGrad::validParams();
  params.addClassDescription(
      "Diffusion equation Kernel that takes an isotropic Diffusivity from a material property");
  params.addParam<MaterialPropertyName>(
      "diffusivity", "D", "The diffusivity value or material property");
  params.addCoupledVar("v",
                       "Coupled concentration variable for kernel to operate on; if this "
                       "is not specified, the kernel's nonlinear variable will be used as "
                       "usual");
  // MatDiffusionBase<T>::validParams adds 'args' (not JvarMapInterface's 'coupled_variables')
  params.addCoupledVar("args",
                       "Optional vector of arguments for the diffusivity. If provided and "
                       "diffusivity is a derivative parsed material, Jacobian contributions from "
                       "the diffusivity will be automatically computed");
  return params;
}

KokkosMatDiffusion::KokkosMatDiffusion(const InputParameters & parameters)
  : Moose::Kokkos::KernelGrad(parameters),
    _diffusivity(getKokkosMaterialProperty<Real>("diffusivity")),
    _is_coupled(isCoupled("v")),
    _v_var(_is_coupled ? coupled("v") : _var.number()),
    _grad_v(_is_coupled ? kokkosCoupledGradient("v") : _grad_u),
    _n_args(coupledComponents("args"))
{
  if (_is_coupled && coupledComponents("v") != 1)
    paramError("v", "Exactly one coupled variable is supported.");
  // With v equal to the kernel variable MatDiffusionBase drops the Laplacian from every Jacobian
  // block (the diagonal skips it because v is coupled, and jvar never equals the row variable in
  // the off-diagonal loop); omit v to diffuse the kernel variable itself
  if (_is_coupled && _v_var == _var.number())
    paramError("v", "Omit v to diffuse the kernel variable itself.");

  // Optional derivative handles under the CPU names (DerivativeMaterialPropertyNameInterface):
  // fetched only when declared, otherwise left unset and read as zero on the device
  const auto D = getMaterialPropertyName("diffusivity");
  const auto fetch = [this](const MaterialPropertyName & derivative,
                            Moose::Kokkos::MaterialProperty<Real> & handle)
  {
    if (hasKokkosMaterialPropertyByName<Real>(derivative))
      handle = getKokkosMaterialPropertyByName<Real>(derivative);
  };

  std::vector<VariableName> coupled_names;
  fetch(derivativePropertyNameFirst(D, _var.name()), _dD_du);
  if (_is_coupled)
  {
    const auto vname = coupledName("v", 0);
    fetch(derivativePropertyNameFirst(D, vname), _dD_dv);
    coupled_names.push_back(vname);
  }

  // Always created, possibly empty: the dispatcher deep-copies the property array and aborts on
  // one that was never allocated (KokkosArray.h deepCopy)
  _arg_var.create(_n_args);
  _dD_darg.create(_n_args);
  for (unsigned int k = 0; k < _n_args; ++k)
  {
    if (!isCoupled("args", k))
      paramError("args", "Entries must be variables, not literal values.");
    const auto * arg = getVar("args", k);
    // Only a nonlinear variable of this kernel's system has a Jacobian column; an auxiliary
    // argument keeps its derivative handle unused (as on the CPU) and never matches a jvar
    _arg_var[k] = (arg->kind() == Moose::VAR_SOLVER && &arg->sys() == &_var.sys())
                      ? arg->number()
                      : libMesh::invalid_uint;
    fetch(derivativePropertyNameFirst(D, arg->name()), _dD_darg[k]);
    coupled_names.push_back(arg->name());
  }
  if (_n_args)
    _arg_var.copyToDevice();

  // DerivativeMaterialInterface::validateNonlinearCoupling("diffusivity"), issued at construction
  // (Kokkos materials are created before kernels): a declared dD/dx for a nonlinear x that is
  // neither the kernel variable nor coupled here would be silently dropped from the Jacobian
  std::vector<VariableName> missing;
  for (const auto & x : _var.sys().getVariableNames())
  {
    if (x == _var.name())
      continue;
    if (std::find(coupled_names.begin(), coupled_names.end(), x) != coupled_names.end())
      continue;
    if (hasKokkosMaterialPropertyByName<Real>(derivativePropertyNameFirst(D, x)))
      missing.push_back(x);
  }
  if (!missing.empty())
    mooseWarning("Missing coupled variables {",
                 Moose::stringify(missing),
                 "} (add them to args parameter of ",
                 name(),
                 ")");

  // A misspelled hook is silently replaced by the base default (KokkosDispatcher.h); fail loudly
  if (!Moose::Kokkos::hasUserJacobianHook<KokkosMatDiffusion>() ||
      !Moose::Kokkos::hasUserOffDiagJacobianHook<KokkosMatDiffusion>())
    mooseError("KokkosMatDiffusion Jacobian hooks were not registered with the expected signatures.");
}
