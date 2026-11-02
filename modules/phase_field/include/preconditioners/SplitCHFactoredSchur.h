//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseTypes.h"

#include "libmesh/petsc_macro.h"

#include <petscksp.h>

#include <array>
#include <limits>
#include <string>

class FEProblemBase;
class ConsoleStream;

#if PETSC_RELEASE_GREATER_EQUALS(3, 20, 0)
/**
 * The context of the PCSHELL that preconditions the Schur complement of the split Cahn-Hilliard
 * system (first split the chemical potential w, second split the concentration c, Schur FULL with
 * schur_pre = A11) with a product of two shifted operators. With A00 = -Mm, A01 = H = kappa K + F,
 * A10 = M K and A11 = Mm / dt (implicit Euler, unscaled rows), tau = dt M, eps = sqrt(tau kappa),
 * a = beta eps and g = tau / a:
 *
 *  - consistent_mass: Op1 = dt A11 + (a / M) A10, Fc = A01 - (kappa / M) A10,
 *    Op2 = dt A11 + ((eps^2 / a) / M) A10 + (tau / a) Fc, and y = dt^2 Op2^-1 A11 Op1^-1 x;
 *  - diagonal_mass: D = -diag(A00), Op1 = (a / M) A10 + D, Op2 = g A01 + D, and
 *    y = dt Op2^-1 D Op1^-1 x.
 *
 * Each inverse is one solve of a preonly factor KSP (prefixes <shell prefix>stiffness_factor_ and
 * <shell prefix>hessian_factor_; another KSP type is refused), whose PC must be a fixed linear
 * operator (one V-cycle or an LU). The dt is captured with the blocks at the build and never read by
 * the apply. The object that owns this context (SplitCHFactoredSchurPreconditioner) installs the
 * callbacks from the FieldSplit's post-setup and destroys the context; no PCShellSetDestroy is set.
 * Requires PETSc 3.20 or later.
 */
class SplitCHFactoredSchur
{
public:
  enum class Form
  {
    ConsistentMass,
    DiagonalMass
  };
  enum class BlockCheck
  {
    FirstBuild,
    EveryBuild,
    None
  };

  /// The owner's parameters, fixed for the run
  struct Settings
  {
    std::string name;
    Real mobility;
    Real kappa;
    Real split_parameter;
    Form form;
    BlockCheck block_check;
    Real block_check_tolerance;
    bool verify_apply;
  };

  SplitCHFactoredSchur(FEProblemBase & problem,
                       const ConsoleStream & console,
                       const Settings & settings);
  virtual ~SplitCHFactoredSchur();

  /// Installs the callbacks on the Schur complement sub-solver's PC and marks the context stale
  PetscErrorCode install(PC shell);

  /// The PETSc callbacks
  static PetscErrorCode setUpCallback(PC pc);
  static PetscErrorCode applyCallback(PC pc, Vec x, Vec y);
  static PetscErrorCode viewCallback(PC pc, PetscViewer viewer);

  /// The counters (post-setup calls, first builds, rebuilds, refills, guard builds) in one line,
  /// "<name> <label>: post-setup calls ..."
  std::string counters(const std::string & label = "counters") const;

protected:
  /// Whether the setup or the apply's guard must build (the stale mark of the post-setup)
  virtual bool buildRequested() const { return _stale; }

  /// Whether the shell's setup builds on a request; when it does not, the apply's guard does every
  /// build (a test-application control turns it off to run every build through the guard)
  virtual bool setUpBuilds() const { return true; }

  /// Assembles Op2 (and Fc in the consistent_mass form) at a first build or a rebuild
  /// (refill = false) or in place on a refill (refill = true)
  virtual PetscErrorCode assembleHessianFactor(bool refill);

  /// The dt the apply scales by: the dt captured at the build
  virtual Real applyTimeStep() const { return _dt_build; }

  FEProblemBase & _problem;
  const ConsoleStream & _console;
  const Settings _settings;
  const MPI_Comm _comm;

  /// The blocks held by reference since the last first build or rebuild
  Mat _A00 = nullptr;
  Mat _A01 = nullptr;
  Mat _A10 = nullptr;
  Mat _A11 = nullptr;

  /// The factor matrices, Fc (consistent_mass) and D (diagonal_mass)
  Mat _op1 = nullptr;
  Mat _op2 = nullptr;
  Mat _fc = nullptr;
  Vec _d = nullptr;

  /// The constants of the last build
  Real _dt_build = 0;
  Real _tau = 0;
  Real _eps = 0;
  Real _a = 0;
  Real _g = 0;

  /// Whether factors exist
  bool _built = false;

private:
  enum class Branch
  {
    None = 0,
    FirstBuild = 1,
    Refill = 2,
    Rebuild = 3,
    Guard = 4
  };

  PetscErrorCode setUp(PC pc);
  PetscErrorCode apply(PC pc, Vec x, Vec y);
  PetscErrorCode view(PetscViewer viewer) const;

  /// One build: the checks of every build, the block check, then a first build, refill or rebuild
  PetscErrorCode build(PC pc, bool guard);
  /// The shell's operator without the guard, the counters or the prints; pc may be null (dumps,
  /// where a failed first factor solve stops the run)
  PetscErrorCode applyOperator(PC pc, Vec x, Vec y);

  /// At the first build: A00 and A11 own the same number of rows on every process
  PetscErrorCode checkLocalSizes(Mat A00, Mat A11);
  PetscErrorCode checkTimeStep(Mat A00, Mat A11, Real dt, unsigned int build_index);
  PetscErrorCode blockCheck(Mat A00, Mat A01, Mat A10, unsigned int build_index);
  PetscErrorCode holdBlocks(const std::array<Mat, 4> & blocks);
  PetscErrorCode assembleStiffnessFactor(bool refill);
  PetscErrorCode checkPatterns();
  PetscErrorCode createFactorSolver(PC pc, Mat op, const char * suffix, KSP * ksp) const;
  PetscErrorCode destroyBuild();
  PetscErrorCode blockFingerprint(const std::array<Mat, 4> & blocks,
                                  std::array<Real, 4> & norms,
                                  std::array<Real, 4> & states) const;

  /// verify_apply: the per-apply print, the per-rebuild dump and the reuse dump
  PetscErrorCode verifyApply(PC pc);
  PetscErrorCode writeDump(PC pc, bool reuse, unsigned int index);

  std::string formName() const;

  KSP _k1 = nullptr;
  KSP _k2 = nullptr;
  Vec _t1 = nullptr;
  Vec _t2 = nullptr;

  std::array<PetscObjectId, 4> _ids{};
  std::array<PetscObjectState, 4> _nonzero_states{};
  bool _refillable = false;
  bool _stale = false;

  /// The counters; post-setup calls = first builds + rebuilds + refills + guard builds
  unsigned int _n_post_setup = 0;
  unsigned int _n_first = 0;
  unsigned int _n_rebuild = 0;
  unsigned int _n_refill = 0;
  unsigned int _n_guard = 0;
  unsigned int _n_builds = 0;
  unsigned long _n_apply = 0;
  unsigned long _n_factor_fail = 0;

  /// The last build: the post-setup count it answered, how it was counted and the path it took
  unsigned int _generation = 0;
  Branch _branch = Branch::None;
  Branch _path = Branch::None;

  /// The checks of the last build
  Real _r = 0;
  Real _dt_deviation = 0;
  bool _fit_ran = false;
  Real _beta_fit = 0;
  Real _gamma_fit = 0;
  Real _fpp_min = 0;
  Real _fpp_max = 0;
  Real _fpp_negative_fraction = 0;

  /// verify_apply's records
  std::array<Real, 4> _norm_before{};
  std::array<Real, 4> _norm_after{};
  std::array<Real, 4> _state_before{};
  std::array<Real, 4> _state_after{};
  unsigned int _dumped_post_setup = 0;
  unsigned int _reuse_dumped_generation = std::numeric_limits<unsigned int>::max();
};
#endif
