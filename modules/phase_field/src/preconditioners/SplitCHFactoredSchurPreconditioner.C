//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SplitCHFactoredSchurPreconditioner.h"

#include "Conversion.h"
#include "FEProblemBase.h"
#include "ImplicitEuler.h"
#include "MooseUtils.h"
#include "MooseVariableFieldBase.h"
#include "NonlinearSystemBase.h"
#include "Split.h"
#include "TimeIntegrator.h"

#include <petscksp.h>

registerMooseObject("PhaseFieldApp", SplitCHFactoredSchurPreconditioner);

InputParameters
SplitCHFactoredSchurPreconditioner::validParams()
{
  InputParameters params = FieldSplitPreconditioner::validParams();
  params.addClassDescription(
      "Schur complement FieldSplit of the split Cahn-Hilliard system (first split the chemical "
      "potential chem_pot_c, second the concentration c, schur_pre = A11) whose Schur complement "
      "sub-solver is a PCSHELL applying a product of two shifted operators built from the "
      "FieldSplit's blocks, the mobility, kappa and the dt of the build. Each of its two factor "
      "solvers (the options prefixes stiffness_factor_ and hessian_factor_ under the concentration "
      "split's) must be preonly, and its PC a fixed linear operator: one V-cycle or an LU, no "
      "tolerance-controlled PC, so that the outer Krylov method sees one linear operator.");
  params.addRequiredRangeCheckedParam<Real>(
      "mobility",
      "mobility > 0",
      "The constant mobility M of the concentration equation, checked against the blocks "
      "(block_check)");
  params.addRequiredRangeCheckedParam<Real>(
      "kappa",
      "kappa > 0",
      "The constant gradient energy coefficient kappa of the chemical potential equation, checked "
      "against the blocks (block_check)");
  MooseEnum form("consistent_mass diagonal_mass", "consistent_mass");
  params.addParam<MooseEnum>(
      "schur_form",
      form,
      "consistent_mass: the factors' mass is dt A11 (dt^2 Op2^-1 A11 Op1^-1); diagonal_mass: the "
      "factors' mass is D = -diag(A00) (dt Op2^-1 D Op1^-1)");
  params.addRangeCheckedParam<Real>("split_parameter",
                                    1.0,
                                    "split_parameter > 0",
                                    "The split parameter beta, a = beta sqrt(dt M kappa)");
  MooseEnum block_check("first_build every_build none", "first_build");
  params.addParam<MooseEnum>(
      "block_check",
      block_check,
      "When the fit A01 = beta_fit A10 + gamma A00 against kappa / mobility and the nodal f'' "
      "statistics run: at the run's first build, at every build, or never. The diagonal checks "
      "(the signs, r and the dt ratio) run at every build regardless");
  params.addRangeCheckedParam<Real>("block_check_tolerance",
                                    1e-2,
                                    "block_check_tolerance > 0",
                                    "Relative tolerance on beta_fit against kappa / mobility");
  params.addParam<bool>(
      "verify_apply",
      false,
      "A diagnostic for tests, off by default: after every FieldSplit rebuild k, at the first "
      "apply that follows it, write the blocks, the factor matrices, the factor applications and "
      "the shell's operator to <file_base>_verify_apply_<k>.bin; print one line per apply; and at "
      "the first apply whose dt differs from the build's, write the operator to "
      "<file_base>_verify_apply_<k>_reuse.bin. The files are written under the output file base "
      "(Outputs/file_base). At most 2000 concentration rows; refuses an existing file");
  params.addParam<bool>(
      "print_counters",
      false,
      "Also print the counter line at the start of every solve after the first, labeled "
      "'counters after solve <n>' (the end-of-run line is printed regardless), so a run stopped "
      "before its end keeps the counts of its finished solves");
  return params;
}

SplitCHFactoredSchurPreconditioner::SplitCHFactoredSchurPreconditioner(
    const InputParameters & parameters)
  : FieldSplitPreconditioner(parameters)
#if PETSC_RELEASE_LESS_THAN(3, 20, 0)
{
  mooseError("requires PETSc 3.20 or later: it installs its shell from the FieldSplit's "
             "post-setup callback (PCSetPostSetUp)");
}
#else
    ,
    _print_counters(processor_id() == 0),
    _print_solve_counters(getParam<bool>("print_counters"))
{
  _settings.name = name();
  _settings.mobility = getParam<Real>("mobility");
  _settings.kappa = getParam<Real>("kappa");
  _settings.split_parameter = getParam<Real>("split_parameter");
  _settings.form = getParam<MooseEnum>("schur_form") == "diagonal_mass"
                       ? SplitCHFactoredSchur::Form::DiagonalMass
                       : SplitCHFactoredSchur::Form::ConsistentMass;
  const auto & block_check = getParam<MooseEnum>("block_check");
  _settings.block_check = block_check == "every_build"
                              ? SplitCHFactoredSchur::BlockCheck::EveryBuild
                              : (block_check == "none" ? SplitCHFactoredSchur::BlockCheck::None
                                                       : SplitCHFactoredSchur::BlockCheck::FirstBuild);
  _settings.block_check_tolerance = getParam<Real>("block_check_tolerance");
  _settings.verify_apply = getParam<bool>("verify_apply");

  checkSplitTree();
  checkTimeIntegrator();
}

SplitCHFactoredSchurPreconditioner::~SplitCHFactoredSchurPreconditioner()
{
  // MOOSE has no end-of-run hook for a preconditioner: the counters are printed here, on the
  // first rank, through Moose::out (the console objects may be gone at teardown)
  if (_context && _print_counters)
    Moose::out << _context->counters() << std::endl;
}

void
SplitCHFactoredSchurPreconditioner::checkSplitTree()
{
  const auto top = _nl.getSplit(_decomposition_split);
  const auto & splitting = top->getParam<std::vector<std::string>>("splitting");
  if (splitting.size() != 2)
    top->paramError("splitting",
                    "The preconditioner '",
                    name(),
                    "' requires exactly two splits under its top split '",
                    _decomposition_split,
                    "': the chemical potential split first, then the concentration split");
  if (top->getParam<MooseEnum>("splitting_type") != "schur")
    top->paramError("splitting_type",
                    "The preconditioner '",
                    name(),
                    "' requires splitting_type = schur under its top split");
  if (top->getParam<MooseEnum>("schur_type") != "full")
    top->paramError("schur_type",
                    "The preconditioner '",
                    name(),
                    "' requires schur_type = full under its top split");
  if (top->getParam<MooseEnum>("schur_pre") != "A11")
    top->paramError("schur_pre",
                    "The preconditioner '",
                    name(),
                    "' requires schur_pre = A11 under its top split (its shell never reads the "
                    "Schur complement's preconditioning matrix)");

  const auto w = _nl.getSplit(splitting[0]);
  const auto c = _nl.getSplit(splitting[1]);
  const auto & w_vars = w->getParam<std::vector<NonlinearVariableName>>("vars");
  const auto & c_vars = c->getParam<std::vector<NonlinearVariableName>>("vars");
  if (w_vars.size() != 1 || w_vars[0] != "chem_pot_c" || c_vars.size() != 1 || c_vars[0] != "c")
    top->paramError("splitting",
                    "The preconditioner '",
                    name(),
                    "' must list the chemical potential split first: the first split ('",
                    splitting[0],
                    "') must hold the variable 'chem_pot_c' only and the second ('",
                    splitting[1],
                    "') the variable 'c' only; found '",
                    Moose::stringify(w_vars),
                    "' and '",
                    Moose::stringify(c_vars),
                    "'");
  for (const auto & split : {w, c})
    if (!split->getParam<std::vector<std::string>>("splitting").empty())
      split->paramError("splitting",
                        "The preconditioner '",
                        name(),
                        "' requires the two splits under '",
                        _decomposition_split,
                        "' to have no splits of their own");

  // The factors add the two variables' blocks entry by entry (A11, A10 and A01 into one matrix, the
  // diagonal of A00 onto A10 and A01), which pairs the i-th degree of freedom of c with the i-th of
  // chem_pot_c. The pairing rests on two things checked here: splits that restrict nothing but the
  // variable (no block or side restriction on the top split or the two under it), and the two
  // variables of one finite element type on the same blocks, which puts each pair on one node. The
  // first build also requires the two blocks to own the same number of rows on every process, a
  // count that a broken pairing could still meet
  for (const auto & split : {top, w, c})
  {
    const auto reject = [this, &split](const std::string & param, const bool set)
    {
      if (set)
        split->paramError(param,
                          "The preconditioner '",
                          name(),
                          "' requires its top split '",
                          _decomposition_split,
                          "' and the two splits under it to restrict nothing but the variable (no "
                          "'blocks', 'sides', 'unsides' or 'unside_by_var' restriction), because its "
                          "factors pair the degrees of freedom of 'c' and 'chem_pot_c' node by node");
    };
    reject("blocks", !split->getParam<std::vector<SubdomainName>>("blocks").empty());
    reject("sides", !split->getParam<std::vector<BoundaryName>>("sides").empty());
    reject("unsides", !split->getParam<std::vector<BoundaryName>>("unsides").empty());
    reject("unside_by_var_boundary_name",
           split->isParamValid("unside_by_var_boundary_name") &&
               !split->getParam<std::vector<BoundaryName>>("unside_by_var_boundary_name").empty());
    reject("unside_by_var_var_name",
           split->isParamValid("unside_by_var_var_name") &&
               !split->getParam<std::vector<NonlinearVariableName>>("unside_by_var_var_name")
                    .empty());
  }

  const auto & c_var = _nl.getVariable(0, "c");
  const auto & w_var = _nl.getVariable(0, "chem_pot_c");
  if (c_var.feType() != w_var.feType() || c_var.blockIDs() != w_var.blockIDs())
    mooseError("requires the variables 'c' and 'chem_pot_c' to have the same finite element type "
               "and the same block restriction, because its factors add the two variables' "
               "blocks entry by entry; found ",
               Moose::stringify(c_var.feType().family),
               " of order ",
               Moose::stringify(c_var.feType().order),
               " on the blocks ",
               Moose::stringify(c_var.blockIDs()),
               " for 'c' and ",
               Moose::stringify(w_var.feType().family),
               " of order ",
               Moose::stringify(w_var.feType().order),
               " on the blocks ",
               Moose::stringify(w_var.blockIDs()),
               " for 'chem_pot_c'");
  _c_split = splitting[1];
}

void
SplitCHFactoredSchurPreconditioner::checkTimeIntegrator() const
{
  // The integrator exists here: setup_time_integrator runs before add_preconditioning
  const auto & c = _nl.getVariable(0, "c");
  const TimeIntegrator * const ti = _nl.queryTimeIntegrator(c.number());
  if (!dynamic_cast<const ImplicitEuler *>(ti))
    mooseError("requires the implicit Euler time integrator on the concentration 'c' "
               "('Executioner/scheme = implicit-euler', the default), because its factors take the "
               "mass as dt A11; found ",
               ti ? ti->type() : std::string("none (a steady problem)"));
}

void
SplitCHFactoredSchurPreconditioner::initialSetup()
{
  // The base call stores this block's PETSc options and applies its solve_type, if any
  FieldSplitPreconditioner::initialSetup();
  if (_fe_problem.solverParams(_nl.number())._type != Moose::ST_NEWTON)
    mooseError("requires solve_type = NEWTON ('Executioner/solve_type' or this block's "
               "'solve_type'): its factors are built from the assembled Jacobian's blocks");
}

void
SplitCHFactoredSchurPreconditioner::checkShellOptions() const
{
  // Read back from the options database, where the command line wins over the deck
  const std::string prefix = "-" + _nl.prefix() + "fieldsplit_" + _c_split + "_";
  const std::vector<std::pair<std::string, std::string>> required = {{"ksp_type", "preonly"},
                                                                     {"pc_type", "shell"}};
  for (const auto & [key, value] : required)
  {
    const std::string option = prefix + key;
    char buffer[256];
    PetscBool set = PETSC_FALSE;
    LibmeshPetscCall(
        PetscOptionsGetString(nullptr, nullptr, option.c_str(), buffer, sizeof(buffer), &set));
    if (!set || value != buffer)
      mooseError("requires '",
                 option,
                 " ",
                 value,
                 "' under the concentration split '",
                 _c_split,
                 "' (its petsc_options_iname and petsc_options_value, or the command line); found ",
                 set ? "'" + std::string(buffer) + "'" : std::string("the option unset"));
  }

  // The Schur factorization of the top split, under the prefix Split::setup writes it with and
  // the FieldSplit reads it with: PETSc reads the deprecated '_factorization_type' first and then
  // '_fact_type', which overrides it, both without regard to case; with neither set, PETSc's
  // default is full
  const std::string top_prefix = "-" + _nl.prefix() + "pc_fieldsplit_schur_";
  std::string fact_option, fact_value;
  for (const char * const key : {"factorization_type", "fact_type"})
  {
    const std::string option = top_prefix + key;
    char buffer[256];
    PetscBool set = PETSC_FALSE;
    LibmeshPetscCall(
        PetscOptionsGetString(nullptr, nullptr, option.c_str(), buffer, sizeof(buffer), &set));
    if (set)
    {
      fact_option = option;
      fact_value = buffer;
    }
  }
  if (!fact_option.empty() && MooseUtils::toLower(fact_value) != "full")
    mooseError("requires the full Schur factorization under its top split '",
               _decomposition_split,
               "' (its schur_type = full, written as '",
               top_prefix,
               "fact_type full'), because its shell preconditions the Schur complement of the full "
               "factorization; found '",
               fact_option,
               " ",
               fact_value,
               "' (a petsc_options entry of the deck or the command line overrides schur_type)");
}

std::unique_ptr<SplitCHFactoredSchur>
SplitCHFactoredSchurPreconditioner::buildContext(const SplitCHFactoredSchur::Settings & settings)
{
  return std::make_unique<SplitCHFactoredSchur>(_fe_problem, _console, settings);
}

void
SplitCHFactoredSchurPreconditioner::setupDM()
{
  // With print_counters, the counters as they stand after the solves so far, at the start of
  // every solve after the first, through Moose::out with a flush: a run stopped before its end (a
  // time limit, an abort) never reaches the destructor's line
  if (_print_solve_counters && _print_counters && _context)
    Moose::out << _context->counters("counters after solve " + std::to_string(_n_solves))
               << std::endl;
  ++_n_solves;

  FieldSplitPreconditioner::setupDM();
  checkShellOptions();
  if (!_context)
    _context = buildContext(_settings);

  // The post-setup runs after every FieldSplit rebuild; it is set again at every solve, since
  // the SNES (and its PC) is new after a solve without reuse
  PC pc;
  LibmeshPetscCall(KSPGetPC(getKSP(), &pc));
  LibmeshPetscCall(PCSetApplicationContext(pc, this));
  LibmeshPetscCall(PCSetPostSetUp(pc, &SplitCHFactoredSchurPreconditioner::postSetUpCallback));
}

PetscErrorCode
SplitCHFactoredSchurPreconditioner::postSetUpCallback(PC pc)
{
  PetscFunctionBeginUser;
  void * ctx = nullptr;
  PetscCall(PCGetApplicationContext(pc, &ctx));
  PetscCall(static_cast<SplitCHFactoredSchurPreconditioner *>(ctx)->postSetUp(pc));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
SplitCHFactoredSchurPreconditioner::postSetUp(PC pc)
{
  PetscFunctionBeginUser;
  PetscBool is_fieldsplit = PETSC_FALSE;
  PetscCall(PetscObjectTypeCompare((PetscObject)pc, PCFIELDSPLIT, &is_fieldsplit));
  if (!is_fieldsplit)
    mooseError("found a nonlinear system PC that is not a FieldSplit; a '-pc_type' option "
               "replaced the split tree");
  PCCompositeType type;
  PetscCall(PCFieldSplitGetType(pc, &type));
  if (type != PC_COMPOSITE_SCHUR)
    mooseError("found a FieldSplit that is not of type Schur; a '-pc_fieldsplit_type' option "
               "replaced the top split's splitting_type");
  PCFieldSplitSchurPreType schur_pre;
  PetscCall(PCFieldSplitGetSchurPre(pc, &schur_pre, nullptr));
  if (schur_pre != PC_FIELDSPLIT_SCHUR_PRE_A11)
    mooseError("found a FieldSplit whose Schur complement preconditioning matrix is built from '",
               PCFieldSplitSchurPreTypes[schur_pre],
               "', not from A11; a '-pc_fieldsplit_schur_precondition' option replaced the top "
               "split's schur_pre = A11. The shell never reads that matrix: another choice is "
               "either unused or has PETSc form one at every FieldSplit rebuild (Sp for selfp, "
               "the explicit Schur complement for full)");

  // The Schur complement sub-solver is the second of the FieldSplit's KSPs
  PetscInt n = 0;
  KSP * subksp = nullptr;
  PetscCall(PCFieldSplitSchurGetSubKSP(pc, &n, &subksp));
  KSP schur = subksp[1];
  PetscCall(PetscFree(subksp));
  PC shell;
  PetscCall(KSPGetPC(schur, &shell));
  PetscBool is_shell = PETSC_FALSE;
  PetscCall(PetscObjectTypeCompare((PetscObject)shell, PCSHELL, &is_shell));
  if (!is_shell)
  {
    PCType shell_type = nullptr;
    PetscCall(PCGetType(shell, &shell_type));
    mooseError("found the Schur complement sub-solver's PC of type '",
               shell_type ? shell_type : "unset",
               "', not shell: the concentration split '",
               _c_split,
               "' needs '-pc_type shell'");
  }
  PetscCall(_context->install(shell));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif
