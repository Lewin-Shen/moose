#!/usr/bin/env python3
"""Reference check of the verify_apply dumps of SplitCHFactoredSchurPreconditioner.

Usage (from the folder that holds the dumps):
  python3 -I verify_apply.py FILE_BASE (--build K | --all-builds) --mobility M --kappa KAPPA --beta BETA
      --form {consistent_mass,diagonal_mass} --mode {exact,approximate}
      [--fit] [--expect-fpp F] [--require-branches LIST] [--require-reuse] [--require-builds N]

It forms the shell's operator independently from the dumped blocks A00 = -Mm, A01 = H = kappa K + F,
A10 = M K and A11 = Mm / dt, with a = beta sqrt(dt M kappa) and g = dt M / a, and checks each dump:
  - always: the dump's implied dt, -diag(A00)_i / diag(A11)_i on every row, against the dt captured at the
    build within 1e-10 (a dump from a step that cut back holds the matrices after the cut); the lifecycle (the
    rebuild k answered by a build, the counters balanced, the blocks' Frobenius norms and nonzero states equal
    bitwise before and after the build and at the dump, the captured dt equal to the dt at the dump); the
    factor matrices against those formed from the blocks and the captured dt within 1e-10; the operator
    against dt^2 P2 A11 P1 (diagonal_mass: dt P2 D P1), P1 and P2 the dumped factor applications, within
    1e-10; the operator on the random vectors against its unit-vector columns within 1e-10; linearity,
    |P(x + 2y) - Px - 2Py| / |Px| at most 1e-12;
  - --mode exact (factors preonly lu): the operator against the dense inverse of
    (1/dt)(Mm + aK) Mm^-1 (Mm + gH) (diagonal_mass: D in place of Mm, D = diag(-A00)) within 1e-10, and the
    identity S~ - S = (1/dt)(aK + gH) against S = A11 - A10 A00^-1 A01 (diagonal_mass:
    S^ - Sp = (1/dt)(aK + gH + D - Mm), Sp = A11 - A10 diag(A00)^-1 A01, Mm = dt A11) within 1e-10 of |S|;
  - --fit: the least-squares fit A01 = beta_fit A10 + gamma A00 over A01's stored entries and the nodal
    f'' = rowsum(A01) / (-rowsum(A00)); beta_fit against kappa/M within the dump's block_check_tolerance, and
    against the object's own fit within 1e-10 when the fit ran at that build; with --expect-fpp F, beta_fit
    equal to kappa/M and -gamma equal to F within 1e-10 (the constant-f'' deck);
  - --all-builds: every build dump k = 1, 2, ... and every reuse dump against its build's operator within
    1e-10; --require-branches (comma list of refill, rebuild, guard), --require-reuse and --require-builds
    require what a case must have run (--require-reuse also fails on a rebuild from scratch, which voids a
    retained-build case; --require-builds N fails unless the run made exactly N builds: N build dumps and N
    builds counted at the last dump).
No tolerance is loosened for V-cycle factors. Exit status 0 when every check passes, 1 otherwise.

The dump layout (PETSc binary: big-endian; integers of 32 or 64 bits, the width read from the first class id
as PETSc's own reader does; real 8-byte double scalars, checked, because the format does not record the
scalar type; written collectively, so one and four ranks read the same way; no .info file):
  build dump FILE_BASE_verify_apply_<k>.bin, k the FieldSplit rebuild (post-setup) count:
    the header Vec (64 reals); Mats A00, A01, A10, A11 (the blocks the FieldSplit holds at the dump); Mats
    Op1, Op2; Mat Fc (consistent_mass) or Vec D (diagonal_mass); n_c Vecs P1 e_j; n_c Vecs P2 e_j; n_c Vecs
    P e_j (the shell's operator); Vecs x, y (seeded random), P x, P y, P (x + 2 y).
  reuse dump FILE_BASE_verify_apply_<g>_reuse.bin, g the generation of the build applied:
    the header; n_c Vecs P e_j; Vecs x, y, P x, P y, P (x + 2 y).
  header: 0 version (1), 1 the dump's index (k, or g), 2 generation (the post-setup count the last build
    answered), 3 the branch counted (0 none, 1 first build, 2 refill, 3 rebuild, 4 guard build), 4 the path
    taken (0 none, 1 first build, 2 refill, 3 rebuild), 5 form (0 consistent_mass, 1 diagonal_mass), 6 the
    problem's dt at the dump, 7 the dt captured at the build, 8 M, 9 kappa, 10 beta, 11 tau, 12 eps, 13 a,
    14 g, 15 n_c, 16 the number of unit vectors, 17 post-setup calls, 18 first builds, 19 rebuilds,
    20 refills, 21 guard builds, 22 applies before this one, 23-26 the Frobenius norms of A00, A01, A10, A11
    before the last build, 27-30 after it, 31-34 their nonzero states before it, 35-38 after it, 39 r, 40 the
    dt ratio deviation, 41 refillable, 42-43 the seeds of x and y, 44 ranks, 45 kind (0 build, 1 reuse),
    46 beta_fit, 47 gamma_fit, 48 f'' min, 49 f'' max, 50 f'' negative fraction, 51 whether the fit ran at
    the last build, 52-55 the norms at the dump, 56-59 the nonzero states at the dump, 60 failed factor
    solves, 61 builds, 62 block_check_tolerance, 63 block_check (0 first_build, 1 every_build, 2 none).

The PETSc binary reader below is inline, because python3 -I does not put the script's folder on the import
path.
Needs numpy and scipy only.
"""
import argparse
import glob
import math
import os
import re
import sys

import numpy as np
import scipy.sparse as sp

MAT_CLASSID = 1211216
VEC_CLASSID = 1211214
CLASS_IDS = (MAT_CLASSID, VEC_CLASSID)
SCALAR_BYTES = 8  # real double scalars, the only kind this script reads
HEADER_LENGTH = 64
BRANCHES = {0: "none", 1: "first build", 2: "refill", 3: "rebuild", 4: "guard build"}
TOL = 1e-10  # every operator, factor and dt check; each holds to rounding, V-cycle factors too
TOL_LINEAR = 1e-12  # linearity on the random vectors


def dump_format(path):
    """The width in bytes (4 or 8) of a dump's integers, after a check that its scalars are real doubles.

    The width comes from the first class id, the header Vec's, the way PETSc's own reader tells a file
    written with 32-bit integers from one written with 64-bit integers: read as a 32-bit integer, the first
    four bytes are the class id in the first case and 0 in the second. The binary format does not record the
    scalar type, so the scalars are checked through the header: its HEADER_LENGTH entries must end, at 8
    bytes each (not 4 or 16), where the next object's class id begins, and its first entry must read the
    version 1 exactly as a double (which a complex single-precision pair does not).
    """
    with open(path, "rb") as fh:
        data = fh.read(2 * 8 + 16 * HEADER_LENGTH + 8)

    def integer_at(offset, width):
        chunk = data[offset:offset + width]
        return int(np.frombuffer(chunk, f">i{width}")[0]) if len(chunk) == width else None

    if integer_at(0, 4) in CLASS_IDS:
        width = 4
    elif integer_at(0, 8) in CLASS_IDS:
        width = 8
    else:
        raise ValueError(f"{path}: not a PETSc binary file (its first class id reads as neither a 32-bit nor a "
                         f"64-bit Vec or Mat class id)")
    if integer_at(0, width) != VEC_CLASSID or integer_at(width, width) != HEADER_LENGTH:
        raise ValueError(f"{path}: the first object is not the {HEADER_LENGTH}-entry header Vec")
    start = 2 * width
    spans = [s for s in (4, SCALAR_BYTES, 16) if integer_at(start + HEADER_LENGTH * s, width) in CLASS_IDS]
    chunk = data[start:start + SCALAR_BYTES]
    version = float(np.frombuffer(chunk, ">f8")[0]) if len(chunk) == SCALAR_BYTES else math.nan
    if spans != [SCALAR_BYTES] or version != 1.0:
        raise ValueError(f"{path}: the scalars are not real 8-byte doubles (the header's entries end at the next "
                         f"object for {' or '.join(map(str, spans)) or 'no'} bytes per entry, and its first "
                         f"entry reads {version!r} as a double); this script reads the dumps of a PETSc build "
                         f"with real double scalars only")
    return width


def iter_objects(path):
    """Yields ('mat', csr_matrix) or ('vec', ndarray) for each object of a PETSc binary file, in order, with
    every integer field of the width dump_format() reads."""
    width = dump_format(path)
    itype = np.dtype(f">i{width}")
    rtype = np.dtype(">f8")
    with open(path, "rb") as fh:

        def read(count, dtype):
            raw = fh.read(count * dtype.itemsize)
            if len(raw) != count * dtype.itemsize:
                raise ValueError(f"{path}: the file ends inside an object")
            return np.frombuffer(raw, dtype)

        while True:
            head = fh.read(width)
            if not head:
                return
            if len(head) != width:
                raise ValueError(f"{path}: the file ends inside a class id")
            cid = int(np.frombuffer(head, itype)[0])
            if cid == MAT_CLASSID:
                m, n, nz = (int(x) for x in read(3, itype))
                rl = read(m, itype).astype(np.int64)
                ci = read(nz, itype).astype(np.int64)
                va = read(nz, rtype).astype(np.float64)
                indptr = np.concatenate(([0], np.cumsum(rl)))
                yield "mat", sp.csr_matrix((va, ci, indptr), shape=(m, n))
            elif cid == VEC_CLASSID:
                n = int(read(1, itype)[0])
                yield "vec", read(n, rtype).astype(np.float64)
            else:
                raise ValueError(f"{path}: unknown class id {cid} at offset {fh.tell() - width} "
                                 f"({8 * width}-bit integers)")


class Dump:
    """One dump, its header decoded and its objects split by role."""

    def __init__(self, path):
        self.path = path
        objs = list(iter_objects(path))
        kinds = [k for k, _ in objs]
        if not objs or kinds[0] != "vec" or objs[0][1].size != HEADER_LENGTH:
            raise ValueError(f"{path}: no {HEADER_LENGTH}-real header")
        h = objs[0][1]
        self.h = h
        if int(h[0]) != 1:
            raise ValueError(f"{path}: header version {h[0]}, expected 1")
        self.reuse = int(h[45]) == 1
        self.n = int(h[15])
        n = self.n
        pos = 1
        if not self.reuse:
            want = ["mat"] * 6 + ["mat" if int(h[5]) == 0 else "vec"] + ["vec"] * (3 * n + 5)
        else:
            want = ["vec"] * (n + 5)
        if kinds[1:] != want:
            raise ValueError(f"{path}: {len(objs) - 1} objects after the header, not the layout of a "
                             f"{'reuse' if self.reuse else 'build'} dump with n_c = {n}")
        o = [x for _, x in objs[1:]]
        if not self.reuse:
            self.A00, self.A01, self.A10, self.A11, self.Op1, self.Op2 = o[0:6]
            self.Fc = o[6] if int(h[5]) == 0 else None
            self.D = o[6] if int(h[5]) == 1 else None
            pos = 7
            self.P1 = np.column_stack(o[pos:pos + n])
            self.P2 = np.column_stack(o[pos + n:pos + 2 * n])
            pos += 2 * n
        else:
            pos = 0
        self.P = np.column_stack(o[pos:pos + n])
        self.x, self.y, self.Px, self.Py, self.Pz = o[pos + n:pos + n + 5]


def rel(a, b):
    """|a - b| / |b| in the Frobenius (2-) norm, for dense or sparse arrays."""
    def nrm(z):
        return sp.linalg.norm(z) if sp.issparse(z) else np.linalg.norm(z)
    d = a - b
    nb = nrm(b)
    return nrm(d) / nb if nb > 0 else nrm(d)


class Checks:
    def __init__(self):
        self.failed = []

    def check(self, label, ok, detail=""):
        print(f"  {'PASS' if ok else 'FAIL'} {label}{(': ' + detail) if detail else ''}")
        if not ok:
            self.failed.append(label)

    def tol(self, label, value, limit):
        # written so that a NaN fails
        self.check(label, bool(value <= limit), f"{value:.3e} (limit {limit:.0e})")

    @staticmethod
    def info(label, detail):
        print(f"  info {label}: {detail}")


def constants(dt, M, kappa, beta):
    tau = dt * M
    eps = math.sqrt(tau * kappa)
    a = beta * eps
    return tau, eps, a, tau / a


def fit(A00, A01, A10):
    """Least squares A01 = beta A10 + gamma A00 over A01's stored entries; the nodal f''."""
    P = A01.tocoo()
    c10 = np.asarray(A10.tocsr()[P.row, P.col]).ravel()
    c00 = np.asarray(A00.tocsr()[P.row, P.col]).ravel()
    coef, *_ = np.linalg.lstsq(np.vstack([c10, c00]).T, P.data, rcond=None)
    rs = lambda X: np.asarray(X.sum(axis=1)).ravel()
    fnode = rs(A01) / (-rs(A00))
    return float(coef[0]), float(coef[1]), fnode


def check_build(d, k, a, C):
    h = d.h
    form = int(h[5])
    print(f"build dump {d.path}: index {int(h[1])}, generation {int(h[2])}, branch {BRANCHES[int(h[3])]}, "
          f"path {BRANCHES[int(h[4])]}, ranks {int(h[44])}, n_c {d.n}")
    C.check("dump kind and form", (not d.reuse) and form == (0 if a.form == "consistent_mass" else 1),
           f"kind {'reuse' if d.reuse else 'build'}, form {form}")
    C.check("dump parameters equal the arguments",
           all(abs(h[i] - v) <= 1e-14 * abs(v) for i, v in ((8, a.mobility), (9, a.kappa), (10, a.beta))),
           f"M {h[8]:.17g}, kappa {h[9]:.17g}, beta {h[10]:.17g}")
    dt = float(h[7])
    tau, eps, aa, g = constants(dt, a.mobility, a.kappa, a.beta)
    C.tol("the object's a and g against the arguments and the captured dt",
          max(abs(h[13] - aa) / aa, abs(h[14] - g) / g), TOL)

    # the implied dt, first
    q = -d.A00.diagonal() / d.A11.diagonal()
    C.tol(f"implied dt -diag(A00)/diag(A11) against the captured dt {dt:.17g}", float(np.max(np.abs(q - dt)) / dt), TOL)

    # the lifecycle
    C.check("rebuild k answered by a build", int(h[2]) == k and int(h[3]) != 0,
           f"generation {int(h[2])}, branch {BRANCHES[int(h[3])]}")
    if k == 1:
        C.check("the first rebuild takes the first-build path", int(h[4]) == 1, BRANCHES[int(h[4])])
    C.check("counters balance at the dump", int(h[17]) == int(h[18] + h[19] + h[20] + h[21]),
           f"post-setup {int(h[17])} = first {int(h[18])} + rebuild {int(h[19])} + refill {int(h[20])} "
           f"+ guard {int(h[21])}")
    C.check("the blocks' norms and nonzero states equal before and after the build",
           np.array_equal(h[23:27], h[27:31]) and np.array_equal(h[31:35], h[35:39]))
    C.check("the blocks' norms and nonzero states unchanged from the build to the dump",
           np.array_equal(h[27:31], h[52:56]) and np.array_equal(h[35:39], h[56:60]))
    C.check("the captured dt equals the dt at the dump", h[6] == h[7], f"{h[6]:.17g} and {h[7]:.17g}")
    C.info("r and dt ratio deviation at the build", f"{h[39]:.17g}, {h[40]:.3e}")

    # the factor matrices, formed from the blocks and the captured dt
    M, kap = a.mobility, a.kappa
    A00, A01, A10, A11 = d.A00, d.A01, d.A10, d.A11
    if form == 0:
        op1 = dt * A11 + (aa / M) * A10
        fc = A01 - (kap / M) * A10
        op2 = dt * A11 + ((eps * eps / aa) / M) * A10 + (tau / aa) * fc
        C.tol("Op1 = dt A11 + (a/M) A10", rel(d.Op1, op1), TOL)
        C.tol("Fc = A01 - (kappa/M) A10", rel(d.Fc, fc), TOL)
        C.tol("Op2 = dt A11 + ((eps^2/a)/M) A10 + (tau/a) Fc", rel(d.Op2, op2), TOL)
    else:
        D = -A00.diagonal()
        op1 = (aa / M) * A10 + sp.diags(D)
        op2 = g * A01 + sp.diags(D)
        C.tol("D = -diag(A00)", rel(d.D, D), TOL)
        C.tol("Op1 = (a/M) A10 + D", rel(d.Op1, op1), TOL)
        C.tol("Op2 = g A01 + D", rel(d.Op2, op2), TOL)

    # the operator against its factor applications
    P = d.P
    if form == 0:
        C.tol("P against dt^2 P2 A11 P1", rel(P, dt * dt * (d.P2 @ (A11 @ d.P1))), TOL)
    else:
        C.tol("P against dt P2 D P1", rel(P, dt * (d.P2 @ (d.D[:, None] * d.P1))), TOL)
    check_vectors(d, C)

    if a.mode == "exact":
        Mm = -A00.toarray()
        K = A10.toarray() / M
        H = A01.toarray()
        I = np.eye(d.n)
        mass = Mm if form == 0 else np.diag(-A00.diagonal())
        Pref = dt * np.linalg.solve(mass + g * H, mass @ np.linalg.solve(mass + aa * K, I))
        C.tol(f"P against the dense inverse of the {'(1/dt)(Mm + aK) Mm^-1 (Mm + gH)' if form == 0 else '(1/dt)(D + aK) D^-1 (D + gH)'}",
              rel(P, Pref), TOL)
        Sshell = np.linalg.inv(P)
        if form == 0:
            S = A11.toarray() - A10.toarray() @ np.linalg.solve(A00.toarray(), A01.toarray())
            C.tol("identity S~ - S = (1/dt)(aK + gH), relative to |S|",
                  np.linalg.norm((Sshell - S) - (aa * K + g * H) / dt) / np.linalg.norm(S), TOL)
        else:
            Sp = A11.toarray() - A10.toarray() @ (A01.toarray() / A00.diagonal()[:, None])
            Dm = np.diag(-A00.diagonal())
            C.tol("identity S^ - Sp = (1/dt)(aK + gH + D - Mm), Mm = dt A11, relative to |Sp|",
                  np.linalg.norm((Sshell - Sp) - ((aa * K + g * H + Dm) / dt - A11.toarray())) / np.linalg.norm(Sp),
                  TOL)

    if a.fit:
        beta, gamma, fnode = fit(A00, A01, A10)
        expected = kap / M
        dev = abs(beta - expected) / expected
        C.info("fit", f"beta_fit {beta:.17g} against kappa/M {expected:.17g} (relative {dev:.3e}), gamma {gamma:.17g}; "
               f"nodal f'' in [{fnode.min():.6g}, {fnode.max():.6g}], negative fraction {np.mean(fnode < 0):.4f}")
        C.tol("beta_fit against kappa/M within the dump's block_check_tolerance", dev, float(h[62]))
        if int(h[51]) == 1:
            C.tol("the object's beta_fit against the script's", abs(h[46] - beta) / abs(beta), TOL)
            C.tol("the object's gamma_fit against the script's", abs(h[47] - gamma) / abs(gamma), TOL)
            C.info("the object's nodal f''", f"[{h[48]:.6g}, {h[49]:.6g}], negative fraction {h[50]:.4f}")
        else:
            C.info("the object's fit", "did not run at this build (block_check)")
        if a.expect_fpp is not None:
            C.tol("beta_fit equal to kappa/M (constant f'')", dev, TOL)
            C.tol(f"-gamma equal to f'' = {a.expect_fpp}", abs(-gamma - a.expect_fpp) / abs(a.expect_fpp), TOL)
    return int(h[4])


def check_vectors(d, C):
    """The random vectors against the unit-vector columns, and linearity."""
    C.tol("P x against the unit-vector columns", rel(d.Px, d.P @ d.x), TOL)
    C.tol("P y against the unit-vector columns", rel(d.Py, d.P @ d.y), TOL)
    C.tol("P (x + 2y) against the unit-vector columns", rel(d.Pz, d.P @ (d.x + 2 * d.y)), TOL)
    C.tol("linearity |P(x + 2y) - Px - 2Py| / |Px|",
          np.linalg.norm(d.Pz - d.Px - 2 * d.Py) / np.linalg.norm(d.Px), TOL_LINEAR)


def check_reuse(r, builds, C):
    h = r.h
    gen = int(h[1])
    print(f"reuse dump {r.path}: generation {gen}, captured dt {h[7]:.17g}, dt at the dump {h[6]:.17g} "
          f"(ratio {h[6] / h[7]:.6g})")
    C.check("dump kind", r.reuse)
    C.check("generation unchanged since the build", int(h[2]) == gen)
    C.check("the dt at the dump differs from the captured dt", h[6] != h[7])
    if gen not in builds:
        C.check(f"the build dump of generation {gen} exists", False)
        return
    b = builds[gen]
    C.check("the captured dt equals its build's", h[7] == b.h[7], f"{h[7]:.17g} and {b.h[7]:.17g}")
    C.tol("the reused operator against its build's", rel(r.P, b.P), TOL)
    check_vectors(r, C)


def load(path, C):
    """The dump at path, or None after a failed check that says why it does not read as one."""
    try:
        return Dump(path)
    except ValueError as e:
        C.check(f"{os.path.basename(path)} reads as a dump", False, str(e))
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("file_base")
    sel = ap.add_mutually_exclusive_group(required=True)
    sel.add_argument("--build", type=int)
    sel.add_argument("--all-builds", action="store_true")
    ap.add_argument("--mobility", type=float, required=True)
    ap.add_argument("--kappa", type=float, required=True)
    ap.add_argument("--beta", type=float, required=True)
    ap.add_argument("--form", choices=("consistent_mass", "diagonal_mass"), required=True)
    ap.add_argument("--mode", choices=("exact", "approximate"), required=True)
    ap.add_argument("--fit", action="store_true")
    ap.add_argument("--expect-fpp", type=float)
    ap.add_argument("--require-branches", default="")
    ap.add_argument("--require-reuse", action="store_true")
    ap.add_argument("--require-builds", type=int, metavar="N")
    a = ap.parse_args()
    if a.expect_fpp is not None and not a.fit:
        ap.error("--expect-fpp needs --fit")
    if a.require_builds is not None and not a.all_builds:
        ap.error("--require-builds needs --all-builds")

    C = Checks()
    pattern = re.compile(re.escape(os.path.basename(a.file_base)) + r"_verify_apply_(\d+)\.bin$")
    found = sorted(int(m.group(1)) for f in glob.glob(a.file_base + "_verify_apply_*.bin")
                   for m in [pattern.search(f)] if m)
    ks = [a.build] if a.build is not None else found
    if not ks:
        C.check(f"build dumps {a.file_base}_verify_apply_<k>.bin exist", False)
    elif a.build is None and ks != list(range(1, len(ks) + 1)):
        C.check("the build dumps are numbered 1, 2, ... without a gap", False, str(ks))
    builds, paths = {}, []
    for k in ks:
        path = f"{a.file_base}_verify_apply_{k}.bin"
        if not os.path.exists(path):
            C.check(f"{path} exists", False)
            continue
        d = load(path, C)
        if d is None:
            continue
        builds[k] = d
        paths.append(check_build(d, k, a, C))
    if a.all_builds:
        reuse = sorted(glob.glob(a.file_base + "_verify_apply_*_reuse.bin"))
        reuse_dumps = []
        for f in reuse:
            r = load(f, C)
            if r is not None:
                reuse_dumps.append(r)
                check_reuse(r, builds, C)
        seen = {BRANCHES[p] for p in paths}
        counts = {name: sum(BRANCHES[p] == name for p in paths) for name in ("first build", "refill", "rebuild")}
        C.info("paths taken by the builds", ", ".join(f"{n} {c}" for n, c in counts.items()))
        last = builds[max(builds)].h if builds else None
        if last is not None:
            C.info("counters at the last dump",
                   f"post-setup {int(last[17])}, first {int(last[18])}, rebuild {int(last[19])}, "
                   f"refill {int(last[20])}, guard {int(last[21])}, failed factor solves {int(last[60])}")
        for want in filter(None, a.require_branches.split(",")):
            name = {"refill": "refill", "rebuild": "rebuild", "guard": "guard build"}[want]
            if want == "guard":
                ok = any(int(b.h[3]) == 4 for b in builds.values())
            else:
                ok = name in seen
            C.check(f"at least one {name} (else the case is void)", ok)
        if a.require_reuse:
            C.check("at least one reuse dump (else the case is void)", len(reuse) > 0, f"{len(reuse)} found")
            C.check("no rebuild from scratch (a rebuild voids the retained-build case)", "rebuild" not in seen)
        if a.require_builds is not None:
            # header 61: the builds counted so far, at the dump; the counters only grow
            n = a.require_builds
            counted = max((int(d.h[61]) for d in list(builds.values()) + reuse_dumps), default=0)
            C.check(f"exactly {n} build{'' if n == 1 else 's'} (else the case is void)",
                    len(builds) == n and counted == n,
                    f"{len(builds)} build dumps, {counted} builds counted at the last dump")
    print(f"verify_apply.py: {'PASS' if not C.failed else 'FAIL'} ({len(C.failed)} failed)")
    for f in C.failed:
        print(f"  failed: {f}")
    return 0 if not C.failed else 1


if __name__ == "__main__":
    sys.exit(main())
