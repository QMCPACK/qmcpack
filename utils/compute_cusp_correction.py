#!/usr/bin/env python3

# Compute AO-based cusp correction of Manten and Luchow JCP 115,5362 (2001)
# The parameters for the correction are computed in the script and written
# to the HDF file.

# If the script is run without a -o option, nothing is written to an HDF file

# Short description of algorithm
# 1. Fit a*exp(-alpha*x) + c to 1s and 2s GTO's for r < 0.2  (initial_rc)
# 2. Find points where second derivative of above function meets the 2nd derivative the GTO
#    These are the 'roots'.  From the paper, it looks like the largest root would be the
#    best.  There is a command line option to choose a different root.
# 3. Compute interpolation function and modified cusp function: a*exp(-alpha*x) + b'*x + c'
# 4. Write parameters for these functions to the HDF file


import h5py
import argparse
import shutil
import math

# import read_qmcpack
from sympy import symbols, Symbol, diff, integrate, solve, exp, Eq, Function

# from gaussian_orbitals import GTO
import numpy as np
from scipy import optimize
import dataclasses
from collections import namedtuple


# Extracted and simplified from gaussian_orbitals.py

CG_basis = namedtuple(
    "CG_basis", ["name", "orbtype", "nbasis", "zeta", "contraction_coeff"]
)


def gen_GTO_expr():
    from sympy import symbols, Symbol, diff, sqrt, factorial, pi, S
    from sympy.utilities.lambdify import lambdastr

    x, y, z = symbols("x y z")
    alpha = Symbol("alpha", positive=True, real=True)
    r = Symbol("r", real=True, nonnegative=True)
    i, j, k = symbols("i j k", integer=True)
    N = Symbol("N")

    # gto_sym_raw = N * x**i * y**j * z**k * exp(-alpha *r**2)
    # Only need s-type orbitals, where i=j=k=0
    gto_sym_raw = N * exp(-alpha * r ** 2)
    gto_expr = gto_sym_raw.subs(r ** 2, x * x + y * y + z * z)
    grad = diff(gto_expr, x, 1)
    lap = diff(gto_expr, x, 2)
    dlap = diff(gto_expr, x, 3)
    # Will only evaluate along x
    s = lambdastr((x, y, z, alpha), gto_expr.subs({y: 0, z: 0}))
    ds = lambdastr((x, y, z, alpha), grad.subs({y: 0, z: 0}))
    d2s = lambdastr((x, y, z, alpha), lap.subs({y: 0, z: 0}))
    d3s = lambdastr((x, y, z, alpha), dlap.subs({y: 0, z: 0}))
    print(s)
    print(ds)
    print(d2s)
    print(d3s)

    # Normalization taken from
    # https://pyscf.org/pyscf_api_docs/pyscf.gto.html#pyscf.gto.mole.gto_norm
    l = Symbol("l")
    norm_sym = sqrt(
        S(2) ** (2 * l + 3)
        * factorial(l + 1)
        * (S(2) * alpha) ** (l + 1.5)
        / factorial(2 * l + 2)
        / sqrt(pi)
    )
    print("Normalization: ", lambdastr(alpha, norm_sym.subs(l, 0)))


# The expressions for the following functions were taken from the output of the above code
# (which can be run with --symbolic)


def eval_single_v(x, alpha):
    N = 3.36358566101486 * alpha ** 0.75 / math.pi ** (1 / 4)
    v = N * math.exp(-alpha * x ** 2)
    return v


def eval_single_vgld(x, alpha):
    N = 3.36358566101486 * alpha ** 0.75 / math.pi ** (1 / 4)
    v = N * math.exp(-alpha * x ** 2)
    dv = -2 * N * alpha * x * math.exp(-alpha * x ** 2)
    d2v = 2 * N * alpha * (2 * alpha * x ** 2 - 1) * math.exp(-alpha * x ** 2)
    d3v = -4 * N * alpha ** 2 * x * (2 * alpha * x ** 2 - 3) * math.exp(-alpha * x ** 2)

    return v, dv, d2v, d3v


def eval_contraction_v(x, basis):
    tmp_vals = [eval_single_v(x, basis.zeta[idx]) for idx in range(basis.nbasis)]
    val = np.dot(basis.contraction_coeff, tmp_vals)
    return val


def eval_contraction_vgld(x, basis):
    val = 0.0
    grad = 0.0
    lap = 0.0
    dlap = 0.0
    for idx in range(basis.nbasis):
        c = basis.contraction_coeff[idx]
        v, g, l, d = eval_single_vgld(x, basis.zeta[idx])
        val += c * v
        grad += c * g
        lap += c * l
        dlap += c * d
    return val, grad, lap, dlap


# Extracted from read_qmcpack.py


def read_from_hdf(fname):
    basis_sets = dict()
    f = h5py.File(fname)
    bss = f["basisset"]
    nbelem = bss["NbElements"][0]

    for ib in range(nbelem):
        basis_set = []
        bs_name = "atomicBasisSet" + str(ib)
        # print("bs_name = ",bs_name)
        bs = bss[bs_name]
        ngroups = bs["NbBasisGroups"][0]
        element = bs["elementType"][0].decode("utf-8")
        angular = bs["angular"][0]
        for ig in range(ngroups):
            bg_name = "basisGroup" + str(ig)
            bg = bs[bg_name]
            nradfunc = bg["NbRadFunc"][0]
            # print("reading basis group",bg_name," nradfunc =",nradfunc)
            ang_mom_l = bg["l"][0]
            n_val = bg["n"][0]
            radfunc = bg["radfunctions"]
            rid = bg["rid"][0]
            rtype = bg["type"][0].decode("utf-8")
            if rtype != "Gaussian":
                print("Expecting Gaussian type basisGroup, but got: ", rtype)
            zeta_list = []
            coeff_list = []
            for ir in range(nradfunc):
                radgroup = "DataRad" + str(ir)
                rg = radfunc[radgroup]
                zeta_list.append(rg["exponent"][0])
                coeff_list.append(rg["contraction"][0])
            cg = CG_basis(rid, ang_mom_l, len(zeta_list), zeta_list, coeff_list)
            basis_set.append(cg)

        basis_sets[element] = basis_set

    f.close()
    return basis_sets


def find_s_orbitals(fname):
    f = h5py.File(fname, "r")
    bs = f["basisset"]
    nbasis = bs["NbElements"][0]
    print("nbasis = ", nbasis)
    # bg = f["basisset/atomicBasisSet0/basisGroup0"]

    s_orbs = list()

    for ib in range(nbasis):
        atbs = bs["atomicBasisSet" + str(ib)]
        elem = atbs["elementType"][0].decode("utf-8")
        nb = atbs["NbBasisGroups"][0]
        print(ib, elem, nb)
        for ig in range(nb):
            bas_group = atbs["basisGroup" + str(ig)]
            l = bas_group["l"][0]
            n = bas_group["n"][0]
            if l == 0 and n in [0, 1]:
                print("  ", ig, n, l)
                s_orbs.append((ib, elem, ig, n, l))

    f.close()
    return s_orbs


# Attempt at creating block of symbols
@dataclasses.dataclass
class SmoothSymbols:
    a: object = Symbol("a")
    c: object = Symbol("c")
    r: object = Symbol("r")
    alpha: object = Symbol("alpha", positive=True)
    Delta: object = Symbol("Delta")
    rc: object = Symbol("r_c", positive=True)
    d0: object = Symbol("d_0")
    d1: object = Symbol("d_1")
    d2: object = Symbol("d_2")
    d3: object = Symbol("d_3")
    d4: object = Symbol("d_4")
    d5: object = Symbol("d_5")
    x: object = Symbol("x")
    # GTO value and derivatives at rc + Delta
    g: object = Symbol("g")
    dg: object = Symbol("dg")
    dg2: object = Symbol("dg2")
    dg3: object = Symbol("dg3")

    C1: object = Symbol("C_1")
    C2: object = Symbol("C_2")
    bp: object = Symbol("b'")
    cp: object = Symbol("c'")


# Todo - write the solution as a Python expression so
# sympy is not needed for normal operation of the script


def solve_smoothing():
    ss = SmoothSymbols()
    (
        a,
        c,
        r,
        alpha,
        Delta,
        rc,
        d0,
        d1,
        d2,
        d3,
        d4,
        d5,
        x,
        g,
        dg,
        dg2,
        dg3,
        C1,
        C2,
        bp,
        cp,
    ) = dataclasses.astuple(ss)

    fc = a * exp(-alpha * r) + c

    # Define the second derivative of fc
    dfc = diff(fc, r, 2)

    # Interpolation function
    p = (
        d5 * (x - rc) ** 5
        + d4 * (x - rc) ** 4
        + d3 * (x - rc) ** 3
        + d2 * (x - rc) ** 2
        + d1 * (x - rc)
        + d0
    )
    d1p = diff(p, x)
    d2p = diff(p, x, 2)
    d3p = diff(p, x, 3)

    # Constraints
    # Continuity at rc+Delta (with the GTOs)
    eq1 = Eq(d2p.subs(x, rc + Delta), dg2)
    eq2 = Eq(d3p.subs(x, rc + Delta), dg3)
    # Continuity at rc (with the cusp function fc)
    eq3 = Eq(dfc.subs(r, rc), d2p.subs(x, rc))
    eq4 = Eq(diff(fc, r, 3).subs(r, rc), d3p.subs(x, rc))
    # print('eq1',eq1)
    # print('eq2',eq2)
    # print('eq3',eq3)
    # print('eq4',eq4)

    # Integrate to get the modified cusp function

    dfc1 = integrate(dfc, r) + bp
    eq5 = Eq(d1p.subs(x, rc), dfc1.subs(r, rc))
    eq6 = Eq(d1p.subs(x, rc + Delta), dg)

    dfc2 = integrate(dfc1, r) + cp

    eq7 = Eq(p.subs(x, rc), dfc2.subs(r, rc))
    eq8 = Eq(p.subs(x, rc + Delta), g)

    # sln2 = solve([eq7,eq8],[C2,cp])

    sln = solve(
        [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8], [d0, d1, d2, d3, d4, d5, bp, cp]
    )

    # return ss,sln,sln1,sln2,p2
    return ss, sln, p


def get_reference_vals(basis, rc=0.2):
    npts = 80
    uplim = np.log(rc) / np.log(10)
    # xpts = np.logspace(-2.0, uplim, npts)
    # Get better results if the fit starts from 0.1 rather than a smaller value.
    xpts = np.linspace(0.1, rc, npts)
    vals = np.zeros(npts)

    for i, x in enumerate(xpts):
        v = eval_contraction_v(x, basis)
        vals[i] = v

    return xpts, vals


# Initial cusp correction function
def target_f(x, a, alpha, c):
    return a * np.exp(-alpha * x) + c


# Derivatives for root finding
def func_diff(x, a, alpha, c, basis):
    v, g, l, d = eval_contraction_vgld(x, basis)
    df2 = alpha * alpha * a * np.exp(-alpha * x)
    return l - df2


def func_diff_abs(x, a, alpha, c, basis):
    return np.abs(func_diff(x, a, alpha, c, basis))


def dfunc_diff(x, a, alpha, c, basis):
    v, g, l, d = eval_contraction_vgld(x, basis)
    df3 = -alpha * alpha * alpha * a * np.exp(-alpha * x)
    return d - df3


def compute_smoothing(basis, rc, a, alpha, c, delta, ss, sln):
    x = rc + delta
    g0, g1, g2, g3 = eval_contraction_vgld(x, basis)

    subs_dict = {ss.a: a, ss.alpha: alpha, ss.c: c, ss.rc: rc, ss.Delta: delta}
    subs_dict[ss.g] = g0
    subs_dict[ss.dg] = g1
    subs_dict[ss.dg2] = g2
    subs_dict[ss.dg3] = g3
    d0 = sln[ss.d0].subs(subs_dict)
    d1 = sln[ss.d1].subs(subs_dict)
    d2 = sln[ss.d2].subs(subs_dict)
    d3 = sln[ss.d3].subs(subs_dict)
    d4 = sln[ss.d4].subs(subs_dict)
    d5 = sln[ss.d5].subs(subs_dict)
    print("  d:", d0, d1, d2, d3, d4, d5)

    bp = sln[ss.bp].subs(subs_dict)

    subs_dict[ss.bp] = bp

    cp = sln[ss.cp].subs(subs_dict)

    subs_dict[ss.cp] = cp
    print("  b',c': ", bp, cp)

    return [float(x) for x in [d0, d1, d2, d3, d4, d5, bp, cp]]


def find_roots(initial_rc, popt, basis):
    nrootgrid = 30
    roots = []
    tol = 1e-5

    root_xpts = np.linspace(0.01, initial_rc, nrootgrid)

    for x0 in root_xpts:
        ri = optimize.root_scalar(
            func_diff, x0=x0, fprime=dfunc_diff, args=(*popt, basis)
        )
        root = ri.root
        exists = False
        for er in roots:
            if abs(er - root) < tol:
                exists = True
                break
        if not exists and root > 0.0 and root < initial_rc:
            roots.append(root)

    roots.sort()

    return roots


def write_param(
    fout, atomic_idx, bs_idx, rc, delta, a, alpha, d0, d1, d2, d3, d4, d5, bp, cp
):
    bg = fout[f"basisset/atomicBasisSet{atomic_idx}/basisGroup{bs_idx}"]
    cusp = bg.create_group("shortrangecusp")

    cusp["rcut"] = rc
    cusp["delta"] = delta
    cusp["alpha"] = alpha
    cusp["a"] = a
    print("d0", d0, type(d0))
    cusp["d0"] = d0
    cusp["d1"] = d1
    cusp["d2"] = d2
    cusp["d3"] = d3
    cusp["d4"] = d4
    cusp["d5"] = d5
    cusp["bp"] = bp
    cusp["cp"] = cp


def compute_cusp_correction(fname_in, s_orbs, fname_out=None, root_idx=None):
    bs = read_from_hdf(fname_in)

    ss, soln, p = solve_smoothing()

    fout = None
    if fname_out:
        fout = h5py.File(fname_out, "a")

    for (ib, elem, bs_idx, n, l) in s_orbs:
        print(f"Processing {elem} orbital {bs_idx} n = {n} l = {l}")
        initial_rc = 0.2
        basis = bs[elem]

        # Find initial_rc adaptively, so it results in at least one root
        for itry in range(5):
            xpts, vals = get_reference_vals(basis[bs_idx], initial_rc)
            res = optimize.curve_fit(target_f, xpts, vals, maxfev=3000)
            popt = res[0]
            a, alpha, c = popt
            print(f" optimized a = {a} alpha = {alpha} c = {c}")

            # Find possible roots to choose for actual r_c
            # to keep the second derivative continuous

            roots = find_roots(initial_rc, popt, basis[bs_idx])

            if len(roots) > 0:
                break

            initial_rc *= 2
            print("No roots found, try rc = ", initial_rc)

        # print('roots',sorted(roots))
        # print('roots',roots)
        print("roots", len(roots))
        for ir, root in enumerate(roots):
            print(" ", ir, root)

        delta = 0.001

        # Should allow an interactive choice option
        # if do_interactive:
        #    root_str = input("Root number: ")
        #    root_idx = int(root_str)
        #    if root_idx < 0 or root_idx > len(roots):
        #        print("Chosen root outside range")
        #        root_idx = -1
        #    rc = roots[root_idx]
        # else:
        #    # Choose the last root.
        #    rc = roots[-1]
        if root_idx is not None:
            if root_idx < 0 or root_idx >= len(roots):
                print(f"Chosen root ({root_idx}) outside range: 0-{len(roots)-1}")
                rc = roots[-1]
            else:
                rc = roots[root_idx]
        else:
            # Choose the last root.
            rc = roots[-1]
        print(f"Using root: {rc}")

        new_param = compute_smoothing(basis[bs_idx], rc, a, alpha, c, delta, ss, soln)

        if fout:
            write_param(fout, ib, bs_idx, rc, delta, a, alpha, *new_param)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute cusp correction to GTO basis sets"
    )
    parser.add_argument("input_file", help="input HDF file with basis set")
    parser.add_argument("-o", "--output-file", help="output HDF file with correction")
    # parser.add_argument("-i","--interactive",help="choose roots interactively",action="store_true")
    parser.add_argument("-r", "--root", help="Root index")
    parser.add_argument(
        "--symbolic", help="Output symbolic derivatives for GTO's", action="store_true"
    )

    args = parser.parse_args()

    if args.symbolic:
        gen_GTO_expr()
        exit(0)

    fname_in = args.input_file
    fname_out = args.output_file
    if fname_out is not None:
        shutil.copyfile(fname_in, fname_out)
    else:
        print("No output file specified, no parameters will be written")

    root_idx = None
    if args.root is not None:
        root_idx = int(args.root)
        print("root idx", root_idx)

    s_orbs = find_s_orbitals(fname_in)
    # print(s_orbs)
    compute_cusp_correction(fname_in, s_orbs, fname_out, root_idx)
