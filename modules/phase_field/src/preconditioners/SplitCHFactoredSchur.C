//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SplitCHFactoredSchur.h"

#include "ConsoleStream.h"
#include "FEProblemBase.h"
#include "MooseApp.h"
#include "MooseError.h"
#include "MooseUtils.h"

#include <petscviewer.h>

#include <cmath>
#include <iomanip>
#include <sstream>

#if PETSC_RELEASE_GREATER_EQUALS(3, 20, 0)
namespace
{
/// The tolerance of the r and per-row dt ratio checks; fixed, not a parameter, because unscaled
/// rows under implicit Euler meet both to rounding
const Real split_ch_schur_dt_tolerance = 1e-10;
/// verify_apply refuses a problem above this many c rows
const PetscInt split_ch_schur_dump_max_rows = 2000;
/// The dump header: its layout version and its length (the layout is in verify_apply.py)
const Real split_ch_schur_dump_version = 1;
const std::size_t split_ch_schur_header_length = 64;
/// The seeds of the dump's two random vectors (offset by 7919 per rank)
const unsigned long split_ch_schur_seed_x = 20261031;
const unsigned long split_ch_schur_seed_y = 20261101;
}

SplitCHFactoredSchur::SplitCHFactoredSchur(FEProblemBase & problem,
                                           const ConsoleStream & console,
                                           const Settings & settings)
  : _problem(problem), _console(console), _settings(settings), _comm(problem.comm().get())
{
}

SplitCHFactoredSchur::~SplitCHFactoredSchur()
{
  // A destructor cannot throw: abort on a PETSc error, as NavierStokesProblem's destructor does
  PetscCallAbort(_comm, destroyBuild());
}

PetscErrorCode
SplitCHFactoredSchur::install(PC shell)
{
  PetscFunctionBeginUser;
  const std::string shell_name = "SplitCHFactoredSchur " + formName();
  PetscCall(PCShellSetContext(shell, this));
  PetscCall(PCShellSetSetUp(shell, &SplitCHFactoredSchur::setUpCallback));
  PetscCall(PCShellSetApply(shell, &SplitCHFactoredSchur::applyCallback));
  PetscCall(PCShellSetName(shell, shell_name.c_str()));
  PetscCall(PCShellSetView(shell, &SplitCHFactoredSchur::viewCallback));
  // The FieldSplit rebuilt its blocks: the next setup (or the apply's guard) builds the factors
  _stale = true;
  ++_n_post_setup;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::setUpCallback(PC pc)
{
  PetscFunctionBeginUser;
  void * ctx = nullptr;
  PetscCall(PCShellGetContext(pc, &ctx));
  PetscCall(static_cast<SplitCHFactoredSchur *>(ctx)->setUp(pc));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::applyCallback(PC pc, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  void * ctx = nullptr;
  PetscCall(PCShellGetContext(pc, &ctx));
  PetscCall(static_cast<SplitCHFactoredSchur *>(ctx)->apply(pc, x, y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::viewCallback(PC pc, PetscViewer viewer)
{
  PetscFunctionBeginUser;
  void * ctx = nullptr;
  PetscCall(PCShellGetContext(pc, &ctx));
  PetscCall(static_cast<const SplitCHFactoredSchur *>(ctx)->view(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::setUp(PC pc)
{
  PetscFunctionBeginUser;
  PetscCall(PCSetFailedReason(pc, PC_NOERROR));
  if (setUpBuilds() && buildRequested())
    PetscCall(build(pc, /*guard=*/false));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::apply(PC pc, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  // A PC failure reason is sticky: this apply reports its own factor solves only
  PetscCall(PCSetFailedReason(pc, PC_NOERROR));
  // The guard: the setup did not run, or did not build, after the post-setup's stale mark
  if (buildRequested())
    PetscCall(build(pc, /*guard=*/true));
  ++_n_apply;
  if (_settings.verify_apply)
    PetscCall(verifyApply(pc));
  PetscCall(applyOperator(pc, x, y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::applyOperator(PC pc, Vec x, Vec y)
{
  PetscFunctionBeginUser;
  const Real s = applyTimeStep();
  KSPConvergedReason reason;

  // _t1 = Op1^-1 x. The factor solvers are preonly (checked at their creation), and a preonly
  // solve reduces its PC's failure over the processes before it sets its reason: the reason is the
  // same on every process, so every process takes the same branch below
  PetscCall(KSPSolve(_k1, x, _t1));
  PetscCall(KSPGetConvergedReason(_k1, &reason));
  if (pc)
  {
    if (reason < 0)
      ++_n_factor_fail;
    PetscCall(KSPCheckSolve(_k1, pc, _t1));
    // A failed first factor ends the apply with y flagged, before the second factor sees its output
    if (reason < 0)
    {
#if PETSC_RELEASE_GREATER_EQUALS(3, 22, 0)
      PetscCall(VecFlag(y, 1));
#else
      PetscCall(VecSetInf(y));
#endif
      PetscFunctionReturn(PETSC_SUCCESS);
    }
  }
  // Without an outer PC (the verify_apply dump) nothing carries the failure on: the run stops
  else if (reason < 0)
    mooseError(_settings.name,
               ": the first factor solve (the stiffness factor, Op1) of the verify_apply "
               "diagnostic dump failed at build ",
               _n_builds,
               " (KSP reason ",
               KSPConvergedReasons[reason],
               ")");

  // _t2 = A11 _t1 (consistent_mass) or D _t1 (diagonal_mass)
  if (_settings.form == Form::ConsistentMass)
    PetscCall(MatMult(_A11, _t1, _t2));
  else
    PetscCall(VecPointwiseMult(_t2, _d, _t1));

  // y = Op2^-1 _t2
  PetscCall(KSPSolve(_k2, _t2, y));
  if (pc)
  {
    PetscCall(KSPGetConvergedReason(_k2, &reason));
    if (reason < 0)
      ++_n_factor_fail;
    PetscCall(KSPCheckSolve(_k2, pc, y));
  }

  // y = dt^2 y (consistent_mass) or dt y (diagonal_mass), with the build's dt
  PetscCall(VecScale(y, _settings.form == Form::ConsistentMass ? s * s : s));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::build(PC pc, bool guard)
{
  PetscFunctionBeginUser;
  Mat S, A00, A01, A10, A11;
  PetscCall(PCGetOperators(pc, &S, nullptr));
  PetscCall(MatSchurComplementGetSubMatrices(S, &A00, nullptr, &A01, &A10, &A11));
  const std::array<Mat, 4> blocks{{A00, A01, A10, A11}};
  const unsigned int build_index = _n_builds + 1;

  // The dt of the Jacobian these blocks were extracted from: the only read of the problem's dt in
  // the build; the apply uses the captured value
  const Real dt = _problem.dt();

  // At the run's first build, before any operation that combines the two variables' blocks: the
  // blocks own the same number of rows on every process (what the pairing of the two variables'
  // degrees of freedom rests on is stated in SplitCHFactoredSchurPreconditioner::checkSplitTree)
  if (build_index == 1)
    PetscCall(checkLocalSizes(A00, A11));

  if (_settings.verify_apply)
    PetscCall(blockFingerprint(blocks, _norm_before, _state_before));

  // The diagonal checks, at every build: the signs of diag(A00) and diag(A11), r and the per-row
  // dt ratio
  PetscCall(checkTimeStep(A00, A11, dt, build_index));

  // The fit of A01 against A10 and A00 and the nodal f'' statistics, per block_check
  _fit_ran = false;
  if (_settings.block_check == BlockCheck::EveryBuild ||
      (_settings.block_check == BlockCheck::FirstBuild && build_index == 1))
    PetscCall(blockCheck(A00, A01, A10, build_index));

  // The refill keys: every block's id and nonzero state equal the stored ones
  bool same_blocks = _built;
  for (std::size_t b = 0; same_blocks && b < blocks.size(); ++b)
  {
    PetscObjectId id;
    PetscObjectState nonzero_state;
    PetscCall(PetscObjectGetId((PetscObject)blocks[b], &id));
    PetscCall(MatGetNonzeroState(blocks[b], &nonzero_state));
    same_blocks = (id == _ids[b] && nonzero_state == _nonzero_states[b]);
  }

  _dt_build = dt;
  _tau = dt * _settings.mobility;
  _eps = std::sqrt(_tau * _settings.kappa);
  _a = _settings.split_parameter * _eps;
  _g = _tau / _a;

  if (same_blocks && _refillable)
  {
    // Refill in place: the first build proved the patterns equal, and the blocks are the same
    // objects with the same patterns
    _path = Branch::Refill;
    PetscCall(assembleStiffnessFactor(/*refill=*/true));
    PetscCall(assembleHessianFactor(/*refill=*/true));
    // KSPSetUp returns at once after a solve unless the operators are set again
    PetscCall(KSPSetOperators(_k1, _op1, _op1));
    PetscCall(KSPSetUp(_k1));
    PetscCall(KSPSetOperators(_k2, _op2, _op2));
    PetscCall(KSPSetUp(_k2));
  }
  else
  {
    // A first build, or a rebuild from scratch (new block objects or a changed pattern): release
    // the old references, factor matrices and KSPs
    _path = _built ? Branch::Rebuild : Branch::FirstBuild;
    PetscCall(destroyBuild());
    PetscCall(holdBlocks(blocks));
    PetscCall(assembleStiffnessFactor(/*refill=*/false));
    PetscCall(assembleHessianFactor(/*refill=*/false));
    PetscCall(checkPatterns());
    PetscCall(createFactorSolver(pc, _op1, "stiffness_factor_", &_k1));
    PetscCall(createFactorSolver(pc, _op2, "hessian_factor_", &_k2));
    PetscCall(MatCreateVecs(_A11, &_t1, &_t2));
  }

  _branch = guard ? Branch::Guard : _path;
  switch (_branch)
  {
    case Branch::FirstBuild:
      ++_n_first;
      break;
    case Branch::Refill:
      ++_n_refill;
      break;
    case Branch::Rebuild:
      ++_n_rebuild;
      break;
    case Branch::Guard:
      ++_n_guard;
      break;
    case Branch::None:
      break;
  }
  _built = true;
  _stale = false;
  _generation = _n_post_setup;
  _n_builds = build_index;

  if (_settings.verify_apply)
  {
    PetscCall(blockFingerprint(blocks, _norm_after, _state_after));
    const char * names[] = {"none", "first build", "refill", "rebuild", "guard build"};
    std::ostringstream os;
    os << std::setprecision(17) << _settings.name << " build " << build_index << " (post-setup "
       << _n_post_setup << "): " << names[static_cast<int>(_branch)] << ", path "
       << names[static_cast<int>(_path)] << ", dt " << _dt_build << ", r " << _r
       << ", dt ratio deviation " << _dt_deviation << ", refillable "
       << (_refillable ? "yes" : "no");
    _console << os.str() << std::endl;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::checkLocalSizes(Mat A00, Mat A11)
{
  PetscFunctionBeginUser;
  // The factors add A10 and A01 into matrices of A11's rows and the diagonal of A00 onto them, so
  // the two blocks must own the same number of rows on every process. The flag is reduced before
  // any error, so every process takes the same branch
  PetscInt m00, m11;
  PetscCall(MatGetLocalSize(A00, &m00, nullptr));
  PetscCall(MatGetLocalSize(A11, &m11, nullptr));
  int equal = (m00 == m11) ? 1 : 0;
  PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, &equal, 1, MPI_INT, MPI_MIN, _comm));
  if (!equal)
  {
    PetscInt n00, n11;
    PetscCall(MatGetSize(A00, &n00, nullptr));
    PetscCall(MatGetSize(A11, &n11, nullptr));
    mooseError(_settings.name,
               ": the blocks of 'chem_pot_c' (A00, ",
               n00,
               " rows) and of 'c' (A11, ",
               n11,
               " rows) must own the same number of rows on every process, because the factors add "
               "them entry by entry; this process owns ",
               m00,
               " rows of A00 and ",
               m11,
               " rows of A11");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::checkTimeStep(Mat A00, Mat A11, Real dt, unsigned int build_index)
{
  PetscFunctionBeginUser;
  Vec d00, d11, q;
  PetscCall(MatCreateVecs(A00, nullptr, &d00));
  PetscCall(MatCreateVecs(A11, nullptr, &d11));
  PetscCall(MatGetDiagonal(A00, d00));
  PetscCall(MatGetDiagonal(A11, d11));

  // The signs, first and at any dt: A00 = -Mm has a negative diagonal and A11 = Mm / dt a positive
  // one. Swapped splits give A00 = Mm / dt and A11 = -Mm, which pass the ratios below at dt = 1
  PetscReal d00_max, d11_min;
  PetscCall(VecMax(d00, nullptr, &d00_max));
  PetscCall(VecMin(d11, nullptr, &d11_min));
  // Written so that a NaN fails
  if (!(d00_max < 0) || !(d11_min > 0))
  {
    PetscCall(VecDestroy(&d11));
    PetscCall(VecDestroy(&d00));
    std::ostringstream os;
    os << std::setprecision(12) << _settings.name << ": the blocks at build " << build_index
       << " have max_i diag(A00)_i = " << d00_max
       << " (must be negative: A00 is the chemical potential's block, -Mm) and min_i diag(A11)_i = "
       << d11_min
       << " (must be positive: A11 is the concentration's block, Mm / dt). The splits look "
          "swapped (the first split must hold 'chem_pot_c' and the second 'c'; a "
          "'dm_moose_fieldsplit_names' option on the command line reorders them), the sign of the "
          "equation of 'c' or 'chem_pot_c' is reversed, or a Dirichlet or other row-replacing "
          "condition sits on 'chem_pot_c'.";
    mooseError(os.str());
  }

  // The ratio r = -(d00 . d11) / (dt d11 . d11)
  PetscScalar num, den;
  PetscCall(VecDot(d00, d11, &num));
  PetscCall(VecDot(d11, d11, &den));
  _r = -PetscRealPart(num) / (dt * PetscRealPart(den));

  // The per-row dt ratio deviation, max_i |-d00_i / d11_i - dt| / dt
  PetscReal deviation;
  PetscCall(VecDuplicate(d11, &q));
  PetscCall(VecPointwiseDivide(q, d00, d11));
  PetscCall(VecScale(q, -1.0 / dt));
  PetscCall(VecShift(q, -1.0));
  PetscCall(VecNorm(q, NORM_INFINITY, &deviation));
  _dt_deviation = deviation;

  PetscCall(VecDestroy(&q));
  PetscCall(VecDestroy(&d11));
  PetscCall(VecDestroy(&d00));

  // Written so that a NaN fails
  if (!(std::abs(_r - 1.0) <= split_ch_schur_dt_tolerance) ||
      !(_dt_deviation <= split_ch_schur_dt_tolerance))
  {
    std::ostringstream os;
    os << std::setprecision(12) << _settings.name << ": the blocks at build " << build_index
       << " are not those of unscaled rows under implicit Euler at the problem's dt = " << dt
       << ": r = -(diag(A00) . diag(A11)) / (dt diag(A11) . diag(A11)) = " << _r
       << " (must equal 1 within 1e-10) and max_i |-diag(A00)_i / diag(A11)_i - dt| / dt = "
       << _dt_deviation
       << " (must be at most 1e-10). The object requires unscaled rows (no variable 'scaling' "
          "and 'automatic_scaling = false'), 'Executioner/scheme = implicit-euler', and no "
          "Dirichlet or other row-replacing condition on 'c' or 'chem_pot_c' (a nodal boundary "
          "condition replaces its rows' diagonals).";
    mooseError(os.str());
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::blockCheck(Mat A00, Mat A01, Mat A10, unsigned int build_index)
{
  PetscFunctionBeginUser;
  // The least-squares fit A01 = beta_fit A10 + gamma_fit A00 over A01's stored entries, by its
  // 2 x 2 normal equations: s = (A10.A10, A10.A00, A00.A00, A10.A01, A00.A01) over those entries
  PetscInt rstart, rend, rstart10, rend10, rstart00, rend00;
  PetscCall(MatGetOwnershipRange(A01, &rstart, &rend));
  PetscCall(MatGetOwnershipRange(A10, &rstart10, &rend10));
  PetscCall(MatGetOwnershipRange(A00, &rstart00, &rend00));
  if (rstart10 != rstart || rend10 != rend || rstart00 != rstart || rend00 != rend)
    mooseError(_settings.name,
               ": the block check needs A00, A01 and A10 with the same row ownership; found [",
               rstart,
               ", ",
               rend,
               "), [",
               rstart10,
               ", ",
               rend10,
               ") and [",
               rstart00,
               ", ",
               rend00,
               ")");

  double s[5] = {0, 0, 0, 0, 0};
  for (PetscInt i = rstart; i < rend; ++i)
  {
    PetscInt n01, n10, n00;
    const PetscInt *c01, *c10, *c00;
    const PetscScalar *v01, *v10, *v00;
    PetscCall(MatGetRow(A01, i, &n01, &c01, &v01));
    PetscCall(MatGetRow(A10, i, &n10, &c10, &v10));
    PetscCall(MatGetRow(A00, i, &n00, &c00, &v00));
    // The rows' column indices are sorted: merge A10's and A00's entries onto A01's
    PetscInt p10 = 0, p00 = 0;
    for (PetscInt k = 0; k < n01; ++k)
    {
      const PetscInt j = c01[k];
      while (p10 < n10 && c10[p10] < j)
        ++p10;
      while (p00 < n00 && c00[p00] < j)
        ++p00;
      const double x10 = (p10 < n10 && c10[p10] == j) ? PetscRealPart(v10[p10]) : 0.0;
      const double x00 = (p00 < n00 && c00[p00] == j) ? PetscRealPart(v00[p00]) : 0.0;
      const double y = PetscRealPart(v01[k]);
      s[0] += x10 * x10;
      s[1] += x10 * x00;
      s[2] += x00 * x00;
      s[3] += x10 * y;
      s[4] += x00 * y;
    }
    PetscCall(MatRestoreRow(A00, i, &n00, &c00, &v00));
    PetscCall(MatRestoreRow(A10, i, &n10, &c10, &v10));
    PetscCall(MatRestoreRow(A01, i, &n01, &c01, &v01));
  }
  PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, s, 5, MPI_DOUBLE, MPI_SUM, _comm));
  const double det = s[0] * s[2] - s[1] * s[1];
  _beta_fit = (s[3] * s[2] - s[1] * s[4]) / det;
  _gamma_fit = (s[0] * s[4] - s[1] * s[3]) / det;

  // The nodal f'' = rowsum(A01) / (-rowsum(A00)): its range and negative fraction
  Vec f, m;
  PetscCall(MatCreateVecs(A01, nullptr, &f));
  PetscCall(MatCreateVecs(A00, nullptr, &m));
  PetscCall(MatGetRowSum(A01, f));
  PetscCall(MatGetRowSum(A00, m));
  PetscCall(VecPointwiseDivide(f, f, m));
  PetscCall(VecScale(f, -1.0));
  PetscReal fmin, fmax;
  PetscCall(VecMin(f, nullptr, &fmin));
  PetscCall(VecMax(f, nullptr, &fmax));
  PetscInt n_local, n_global;
  PetscCall(VecGetLocalSize(f, &n_local));
  PetscCall(VecGetSize(f, &n_global));
  const PetscScalar * fa;
  double negative = 0;
  PetscCall(VecGetArrayRead(f, &fa));
  for (PetscInt i = 0; i < n_local; ++i)
    if (PetscRealPart(fa[i]) < 0)
      negative += 1;
  PetscCall(VecRestoreArrayRead(f, &fa));
  PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, &negative, 1, MPI_DOUBLE, MPI_SUM, _comm));
  PetscCall(VecDestroy(&m));
  PetscCall(VecDestroy(&f));
  _fpp_min = fmin;
  _fpp_max = fmax;
  _fpp_negative_fraction = negative / n_global;
  _fit_ran = true;

  const Real expected = _settings.kappa / _settings.mobility;
  const Real deviation = std::abs(_beta_fit - expected) / expected;
  {
    std::ostringstream os;
    os << std::setprecision(17) << _settings.name << " block check at build " << build_index
       << ": beta_fit = " << _beta_fit << " against kappa/M = " << expected
       << ", relative deviation " << deviation << " (block_check_tolerance "
       << _settings.block_check_tolerance << "); gamma_fit = " << _gamma_fit
       << "; nodal f'' = rowsum(A01) / (-rowsum(A00)) in [" << _fpp_min << ", " << _fpp_max
       << "], negative fraction " << _fpp_negative_fraction;
    _console << os.str() << std::endl;
  }
  // Written so that a NaN fails
  if (!(deviation <= _settings.block_check_tolerance))
  {
    std::ostringstream os;
    os << std::setprecision(12) << _settings.name << ": block check failed at build "
       << build_index << ": beta_fit = " << _beta_fit << " against kappa/M = " << expected
       << " (relative deviation " << deviation << ", block_check_tolerance "
       << _settings.block_check_tolerance
       << "): the parameters 'mobility' and 'kappa' do not match the assembled blocks";
    mooseError(os.str());
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::holdBlocks(const std::array<Mat, 4> & blocks)
{
  PetscFunctionBeginUser;
  for (std::size_t b = 0; b < blocks.size(); ++b)
  {
    PetscCall(PetscObjectReference((PetscObject)blocks[b]));
    PetscCall(PetscObjectGetId((PetscObject)blocks[b], &_ids[b]));
    PetscCall(MatGetNonzeroState(blocks[b], &_nonzero_states[b]));
  }
  _A00 = blocks[0];
  _A01 = blocks[1];
  _A10 = blocks[2];
  _A11 = blocks[3];
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::assembleStiffnessFactor(bool refill)
{
  PetscFunctionBeginUser;
  const Real M = _settings.mobility;
  if (_settings.form == Form::ConsistentMass)
  {
    // Op1 = dt A11 + (a / M) A10
    if (refill)
      PetscCall(MatCopy(_A11, _op1, SAME_NONZERO_PATTERN));
    else
      PetscCall(MatDuplicate(_A11, MAT_COPY_VALUES, &_op1));
    PetscCall(MatScale(_op1, _dt_build));
    PetscCall(MatAXPY(
        _op1, _a / M, _A10, refill ? SAME_NONZERO_PATTERN : UNKNOWN_NONZERO_PATTERN));
  }
  else
  {
    // D = -diag(A00); Op1 = (a / M) A10 + D
    if (!refill)
      PetscCall(MatCreateVecs(_A00, nullptr, &_d));
    PetscCall(MatGetDiagonal(_A00, _d));
    PetscCall(VecScale(_d, -1.0));
    if (refill)
      PetscCall(MatCopy(_A10, _op1, SAME_NONZERO_PATTERN));
    else
      PetscCall(MatDuplicate(_A10, MAT_COPY_VALUES, &_op1));
    PetscCall(MatScale(_op1, _a / M));
    PetscCall(MatDiagonalSet(_op1, _d, ADD_VALUES));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::assembleHessianFactor(bool refill)
{
  PetscFunctionBeginUser;
  const Real M = _settings.mobility;
  const MatStructure pattern = refill ? SAME_NONZERO_PATTERN : UNKNOWN_NONZERO_PATTERN;
  if (_settings.form == Form::ConsistentMass)
  {
    // Op2 = dt A11 + ((eps^2 / a) / M) A10
    if (refill)
      PetscCall(MatCopy(_A11, _op2, SAME_NONZERO_PATTERN));
    else
      PetscCall(MatDuplicate(_A11, MAT_COPY_VALUES, &_op2));
    PetscCall(MatScale(_op2, _dt_build));
    PetscCall(MatAXPY(_op2, (_eps * _eps / _a) / M, _A10, pattern));
    // Fc = A01 - (kappa / M) A10, the f'' term (r = 1, checked at every build)
    if (refill)
      PetscCall(MatCopy(_A01, _fc, SAME_NONZERO_PATTERN));
    else
      PetscCall(MatDuplicate(_A01, MAT_COPY_VALUES, &_fc));
    PetscCall(MatAXPY(_fc, -_settings.kappa / M, _A10, pattern));
    // Op2 = Op2 + (tau / a) Fc
    PetscCall(MatAXPY(_op2, _tau / _a, _fc, pattern));
  }
  else
  {
    // Op2 = g A01 + D
    if (refill)
      PetscCall(MatCopy(_A01, _op2, SAME_NONZERO_PATTERN));
    else
      PetscCall(MatDuplicate(_A01, MAT_COPY_VALUES, &_op2));
    PetscCall(MatScale(_op2, _g));
    PetscCall(MatDiagonalSet(_op2, _d, ADD_VALUES));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::checkPatterns()
{
  PetscFunctionBeginUser;
  // Equal local nonzero counts prove equal patterns here, because each factor's pattern is the
  // union of the blocks' it is made from; a refill is allowed only then
  auto nonzeros = [](Mat A, PetscLogDouble & nz) -> PetscErrorCode
  {
    MatInfo info;
    PetscErrorCode ierr = MatGetInfo(A, MAT_LOCAL, &info);
    nz = info.nz_used;
    return ierr;
  };
  PetscLogDouble n11, n10, n01, n1, n2, nf = 0;
  PetscCall(nonzeros(_A11, n11));
  PetscCall(nonzeros(_A10, n10));
  PetscCall(nonzeros(_A01, n01));
  PetscCall(nonzeros(_op1, n1));
  PetscCall(nonzeros(_op2, n2));
  bool equal;
  if (_settings.form == Form::ConsistentMass)
  {
    PetscCall(nonzeros(_fc, nf));
    equal = (n10 == n11 && n01 == n11 && n1 == n11 && n2 == n11 && nf == n11);
  }
  else
    equal = (n1 == n10 && n2 == n01);
  int flag = equal ? 1 : 0;
  PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, &flag, 1, MPI_INT, MPI_MIN, _comm));
  _refillable = (flag == 1);
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::createFactorSolver(PC pc, Mat op, const char * suffix, KSP * ksp) const
{
  PetscFunctionBeginUser;
  const char * prefix = nullptr;
  PetscCall(PCGetOptionsPrefix(pc, &prefix));
  PetscCall(KSPCreate(PetscObjectComm((PetscObject)op), ksp));
  PetscCall(KSPSetOptionsPrefix(*ksp, prefix));
  PetscCall(KSPAppendOptionsPrefix(*ksp, suffix));
  PetscCall(KSPSetType(*ksp, KSPPREONLY));
  PetscCall(KSPSetOperators(*ksp, op, op));
  PetscCall(KSPSetFromOptions(*ksp));
  // One application of the factor's PC per inverse, so that the outer Krylov method sees the same
  // linear operator at every apply and the apply reads a converged reason that is the same on every
  // process: an option that replaces preonly is refused. That the PC itself is a fixed linear
  // operator (one V-cycle or an LU) is documented, not checked
  PetscBool is_preonly = PETSC_FALSE;
  PetscCall(PetscObjectTypeCompare((PetscObject)*ksp, KSPPREONLY, &is_preonly));
  if (!is_preonly)
  {
    const char * full_prefix = nullptr;
    KSPType type = nullptr;
    PetscCall(KSPGetOptionsPrefix(*ksp, &full_prefix));
    PetscCall(KSPGetType(*ksp, &type));
    const std::string option_prefix = "-" + std::string(full_prefix ? full_prefix : "");
    mooseError(_settings.name,
               ": requires the factor solver under the options prefix '",
               option_prefix,
               "' to be preonly ('",
               option_prefix,
               "ksp_type preonly'), so that the outer Krylov method sees one fixed linear "
               "operator; found '",
               type ? type : "unset",
               "'");
  }
  PetscCall(KSPSetUp(*ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::destroyBuild()
{
  PetscFunctionBeginUser;
  PetscCall(KSPDestroy(&_k1));
  PetscCall(KSPDestroy(&_k2));
  PetscCall(MatDestroy(&_op1));
  PetscCall(MatDestroy(&_op2));
  PetscCall(MatDestroy(&_fc));
  PetscCall(VecDestroy(&_d));
  PetscCall(VecDestroy(&_t1));
  PetscCall(VecDestroy(&_t2));
  // Release the references taken on the blocks
  PetscCall(MatDestroy(&_A00));
  PetscCall(MatDestroy(&_A01));
  PetscCall(MatDestroy(&_A10));
  PetscCall(MatDestroy(&_A11));
  _built = false;
  _refillable = false;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::blockFingerprint(const std::array<Mat, 4> & blocks,
                                       std::array<Real, 4> & norms,
                                       std::array<Real, 4> & states) const
{
  PetscFunctionBeginUser;
  for (std::size_t b = 0; b < blocks.size(); ++b)
  {
    PetscReal norm;
    PetscObjectState state;
    PetscCall(MatNorm(blocks[b], NORM_FROBENIUS, &norm));
    PetscCall(MatGetNonzeroState(blocks[b], &state));
    norms[b] = norm;
    states[b] = static_cast<Real>(state);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::verifyApply(PC pc)
{
  PetscFunctionBeginUser;
  // The problem's dt, read for the print and the reuse dump's trigger only
  const Real dt_now = _problem.dt();
  {
    std::ostringstream os;
    os << std::setprecision(17) << _settings.name << " apply " << _n_apply << ": generation "
       << _generation << ", captured dt " << _dt_build << ", current dt " << dt_now;
    _console << os.str() << std::endl;
  }
  // The per-rebuild dump: at the first apply after post-setup k
  if (_dumped_post_setup != _n_post_setup)
  {
    PetscCall(writeDump(pc, /*reuse=*/false, _n_post_setup));
    _dumped_post_setup = _n_post_setup;
  }
  // The reuse dump: at the first apply of this generation whose dt differs from the captured one
  if (dt_now != _dt_build && _reuse_dumped_generation != _generation)
  {
    PetscCall(writeDump(pc, /*reuse=*/true, _generation));
    _reuse_dumped_generation = _generation;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::writeDump(PC pc, bool reuse, unsigned int index)
{
  PetscFunctionBeginUser;
  // The blocks the FieldSplit holds now, which may differ from those of the last build
  Mat S, A00, A01, A10, A11;
  PetscCall(PCGetOperators(pc, &S, nullptr));
  PetscCall(MatSchurComplementGetSubMatrices(S, &A00, nullptr, &A01, &A10, &A11));
  const std::array<Mat, 4> blocks{{A00, A01, A10, A11}};

  PetscInt n_c;
  PetscCall(MatGetSize(A11, &n_c, nullptr));
  if (n_c > split_ch_schur_dump_max_rows)
    mooseError(_settings.name,
               ": 'verify_apply' refuses a problem of ",
               n_c,
               " c rows (at most ",
               split_ch_schur_dump_max_rows,
               "): it writes the shell's operator applied to every unit vector");

  const std::string path = _problem.getMooseApp().getOutputFileBase() + "_verify_apply_" +
                           std::to_string(index) + (reuse ? "_reuse" : "") + ".bin";
  PetscMPIInt rank, n_ranks;
  PetscCallMPI(MPI_Comm_rank(_comm, &rank));
  PetscCallMPI(MPI_Comm_size(_comm, &n_ranks));
  int exists = 0;
  if (rank == 0)
    exists = MooseUtils::pathExists(path) ? 1 : 0;
  PetscCallMPI(MPI_Bcast(&exists, 1, MPI_INT, 0, _comm));
  if (exists)
    mooseError(_settings.name,
               ": 'verify_apply' refuses to overwrite the existing file '",
               path,
               "'; remove the earlier dumps of this output file base first");

  std::array<Real, 4> norm_now, state_now;
  PetscCall(blockFingerprint(blocks, norm_now, state_now));

  // The header; the layout is documented in verify_apply.py
  std::array<Real, split_ch_schur_header_length> h{};
  const bool answered = (_generation == index);
  h[0] = split_ch_schur_dump_version;
  h[1] = index;
  h[2] = _generation;
  h[3] = (reuse || answered) ? static_cast<int>(_branch) : 0;
  h[4] = (reuse || answered) ? static_cast<int>(_path) : 0;
  h[5] = (_settings.form == Form::ConsistentMass) ? 0 : 1;
  h[6] = _problem.dt();
  h[7] = _dt_build;
  h[8] = _settings.mobility;
  h[9] = _settings.kappa;
  h[10] = _settings.split_parameter;
  h[11] = _tau;
  h[12] = _eps;
  h[13] = _a;
  h[14] = _g;
  h[15] = n_c;
  h[16] = n_c;
  h[17] = _n_post_setup;
  h[18] = _n_first;
  h[19] = _n_rebuild;
  h[20] = _n_refill;
  h[21] = _n_guard;
  h[22] = _n_apply - 1;
  for (std::size_t b = 0; b < 4; ++b)
  {
    h[23 + b] = _norm_before[b];
    h[27 + b] = _norm_after[b];
    h[31 + b] = _state_before[b];
    h[35 + b] = _state_after[b];
    h[52 + b] = norm_now[b];
    h[56 + b] = state_now[b];
  }
  h[39] = _r;
  h[40] = _dt_deviation;
  h[41] = _refillable ? 1 : 0;
  h[42] = split_ch_schur_seed_x;
  h[43] = split_ch_schur_seed_y;
  h[44] = n_ranks;
  h[45] = reuse ? 1 : 0;
  h[46] = _beta_fit;
  h[47] = _gamma_fit;
  h[48] = _fpp_min;
  h[49] = _fpp_max;
  h[50] = _fpp_negative_fraction;
  h[51] = _fit_ran ? 1 : 0;
  h[60] = _n_factor_fail;
  h[61] = _n_builds;
  h[62] = _settings.block_check_tolerance;
  h[63] = static_cast<int>(_settings.block_check);

  // Written collectively in PETSc's binary format, so one rank and four ranks read the same way
  PetscViewer viewer;
  PetscCall(PetscViewerCreate(_comm, &viewer));
  PetscCall(PetscViewerSetType(viewer, PETSCVIEWERBINARY));
  PetscCall(PetscViewerBinarySetSkipInfo(viewer, PETSC_TRUE));
  PetscCall(PetscViewerFileSetMode(viewer, FILE_MODE_WRITE));
  PetscCall(PetscViewerFileSetName(viewer, path.c_str()));

  Vec header;
  PetscCall(VecCreate(_comm, &header));
  PetscCall(VecSetSizes(header,
                        rank == 0 ? static_cast<PetscInt>(split_ch_schur_header_length) : 0,
                        static_cast<PetscInt>(split_ch_schur_header_length)));
  PetscCall(VecSetType(header, VECSTANDARD));
  PetscScalar * ha;
  PetscCall(VecGetArray(header, &ha));
  if (rank == 0)
    for (std::size_t i = 0; i < h.size(); ++i)
      ha[i] = h[i];
  PetscCall(VecRestoreArray(header, &ha));
  PetscCall(VecView(header, viewer));
  PetscCall(VecDestroy(&header));

  if (!reuse)
  {
    // The four blocks, then the factor matrices and Fc (consistent_mass) or D (diagonal_mass)
    for (const auto & block : blocks)
      PetscCall(MatView(block, viewer));
    PetscCall(MatView(_op1, viewer));
    PetscCall(MatView(_op2, viewer));
    if (_settings.form == Form::ConsistentMass)
      PetscCall(MatView(_fc, viewer));
    else
      PetscCall(VecView(_d, viewer));
  }

  Vec e, out;
  PetscInt lo, hi;
  PetscCall(VecDuplicate(_t1, &e));
  PetscCall(VecDuplicate(_t1, &out));
  PetscCall(VecGetOwnershipRange(e, &lo, &hi));
  // Passes: 0 and 1 the factor applications (build dumps only), 2 the shell's operator
  for (int pass = reuse ? 2 : 0; pass < 3; ++pass)
    for (PetscInt j = 0; j < n_c; ++j)
    {
      PetscCall(VecSet(e, 0.0));
      if (j >= lo && j < hi)
        PetscCall(VecSetValue(e, j, 1.0, INSERT_VALUES));
      PetscCall(VecAssemblyBegin(e));
      PetscCall(VecAssemblyEnd(e));
      if (pass == 0)
        PetscCall(KSPSolve(_k1, e, out));
      else if (pass == 1)
        PetscCall(KSPSolve(_k2, e, out));
      else
        PetscCall(applyOperator(nullptr, e, out));
      PetscCall(VecView(out, viewer));
    }

  // Two seeded random vectors x and y and z = x + 2 y, then P x, P y and P z
  Vec x, y, z, px, py, pz;
  PetscCall(VecDuplicate(_t1, &x));
  PetscCall(VecDuplicate(_t1, &y));
  PetscCall(VecDuplicate(_t1, &z));
  PetscCall(VecDuplicate(_t1, &px));
  PetscCall(VecDuplicate(_t1, &py));
  PetscCall(VecDuplicate(_t1, &pz));
  PetscRandom random;
  PetscCall(PetscRandomCreate(_comm, &random));
  PetscCall(PetscRandomSetType(random, PETSCRANDER48));
  PetscCall(PetscRandomSetSeed(random, split_ch_schur_seed_x + 7919 * rank));
  PetscCall(PetscRandomSeed(random));
  PetscCall(VecSetRandom(x, random));
  PetscCall(PetscRandomSetSeed(random, split_ch_schur_seed_y + 7919 * rank));
  PetscCall(PetscRandomSeed(random));
  PetscCall(VecSetRandom(y, random));
  PetscCall(PetscRandomDestroy(&random));
  PetscCall(VecWAXPY(z, 2.0, y, x));
  PetscCall(applyOperator(nullptr, x, px));
  PetscCall(applyOperator(nullptr, y, py));
  PetscCall(applyOperator(nullptr, z, pz));
  for (Vec v : {x, y, px, py, pz})
    PetscCall(VecView(v, viewer));

  for (Vec * v : {&e, &out, &x, &y, &z, &px, &py, &pz})
    PetscCall(VecDestroy(v));
  PetscCall(PetscViewerDestroy(&viewer));

  _console << _settings.name << " verify_apply: wrote " << path << std::endl;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchur::view(PetscViewer viewer) const
{
  PetscFunctionBeginUser;
  PetscBool ascii = PETSC_FALSE;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERASCII, &ascii));
  if (!ascii)
    PetscFunctionReturn(PETSC_SUCCESS);
  if (!_built)
  {
    PetscCall(PetscViewerASCIIPrintf(viewer, "factors not built\n"));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,
                                   "form %s, split parameter %.17g, captured dt %.17g, build "
                                   "generation %u\n",
                                   formName().c_str(),
                                   (double)_settings.split_parameter,
                                   (double)_dt_build,
                                   _generation));
  PetscCall(PetscViewerASCIIPrintf(
      viewer, "stiffness factor (mass plus a K), a = %.17g:\n", (double)_a));
  PetscCall(PetscViewerASCIIPushTab(viewer));
  PetscCall(KSPView(_k1, viewer));
  PetscCall(PetscViewerASCIIPopTab(viewer));
  PetscCall(PetscViewerASCIIPrintf(
      viewer, "Hessian factor (mass plus g H), g = %.17g:\n", (double)_g));
  PetscCall(PetscViewerASCIIPushTab(viewer));
  PetscCall(KSPView(_k2, viewer));
  PetscCall(PetscViewerASCIIPopTab(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

std::string
SplitCHFactoredSchur::counters(const std::string & label) const
{
  const unsigned int sum = _n_first + _n_rebuild + _n_refill + _n_guard;
  std::ostringstream os;
  os << _settings.name << " " << label << ": post-setup calls " << _n_post_setup
     << ", first builds " << _n_first << ", rebuilds " << _n_rebuild << ", refills " << _n_refill
     << ", guard builds " << _n_guard << " (post-setup calls = first builds + rebuilds + "
     << "refills + guard builds: " << (sum == _n_post_setup ? "balanced" : "NOT BALANCED")
     << "), applies " << _n_apply << ", failed factor solves " << _n_factor_fail;
  return os.str();
}

std::string
SplitCHFactoredSchur::formName() const
{
  return _settings.form == Form::ConsistentMass ? "consistent_mass" : "diagonal_mass";
}
#endif
