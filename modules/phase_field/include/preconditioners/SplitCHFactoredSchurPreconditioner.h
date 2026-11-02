//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FieldSplitPreconditioner.h"
#include "SplitCHFactoredSchur.h"

#include "libmesh/petsc_macro.h"

#include <memory>

/**
 * A FieldSplitPreconditioner for the split Cahn-Hilliard system whose Schur complement sub-solver
 * is a PCSHELL applying SplitCHFactoredSchur. The split tree is the base class's: two splits under
 * the top split, the first holding chem_pot_c only and the second c only, Schur FULL with
 * schur_pre = A11, and '-ksp_type preonly -pc_type shell' under the concentration split. After the
 * base class's setupDM() the override installs a post-setup callback on the FieldSplit; the
 * callback installs the shell's callbacks after every FieldSplit rebuild and marks its context
 * stale. The counters are printed when the object is destroyed, at the end of the run, and with
 * print_counters also at the start of every solve after the first. Requires PETSc 3.20 or later.
 */
class SplitCHFactoredSchurPreconditioner : public FieldSplitPreconditioner
{
public:
  static InputParameters validParams();

  SplitCHFactoredSchurPreconditioner(const InputParameters & parameters);

#if PETSC_RELEASE_GREATER_EQUALS(3, 20, 0)
  virtual ~SplitCHFactoredSchurPreconditioner();

  virtual void initialSetup() override;
  virtual void setupDM() override;

protected:
  /// Creates the shell's context; called at the first setupDM(), so a derived class's override
  /// is the one called
  virtual std::unique_ptr<SplitCHFactoredSchur>
  buildContext(const SplitCHFactoredSchur::Settings & settings);

private:
  static PetscErrorCode postSetUpCallback(PC pc);
  PetscErrorCode postSetUp(PC pc);

  /// The split tree the forms assume: the chemical potential split first and the concentration
  /// split second under a Schur FULL top split with schur_pre = A11, no split restricted beyond its
  /// variable, and the two variables of one finite element type on the same blocks, so their
  /// degrees of freedom pair node by node
  void checkSplitTree();
  /// Implicit Euler on c's equation, because the factors take the mass as dt A11
  void checkTimeIntegrator() const;
  /// '-ksp_type preonly -pc_type shell' under the concentration split and the full Schur
  /// factorization under the top split, read back from the options database
  void checkShellOptions() const;

  SplitCHFactoredSchur::Settings _settings;
  /// The name of the concentration split (the second under the top split)
  std::string _c_split;
  std::unique_ptr<SplitCHFactoredSchur> _context;
  /// Whether this rank prints the counters
  const bool _print_counters;
  /// Whether the counters are also printed at the start of every solve after the first
  const bool _print_solve_counters;
  /// The number of solves started so far (setupDM() calls)
  unsigned int _n_solves = 0;
#endif
};
