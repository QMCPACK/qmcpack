
# Evaluate STO's starting from a symbolic representation

from sympy import *
from sympy.utilities.lambdify import lambdastr
from collections import namedtuple, defaultdict

# import numpy as np
from autograd import grad
import autograd.numpy as np
import math


# n, zeta, and contraction_coeff are lists of size nbasis
CG_basis = namedtuple(
    "CG_basis", ["orbtype", "nbasis", "n", "zeta", "contraction_coeff"]
)


class STO:
    def __init__(self):
        x, y, z = symbols("x y z", real=True)
        zeta = Symbol("zeta", positive=True, real=True)
        r = Symbol("r", real=True, nonnegative=True)
        N = Symbol("N")
        self.N = N
        n = Symbol("n", integer=True, positive=True)
        norm = (2 * zeta) ** n * sqrt(2 * zeta / factorial(2 * n))

        sto_sym_raw = N * r ** (n - 1) * exp(-zeta * r)

        stosym = sto_sym_raw.subs(N, norm)

        nmax = 3
        self.sto = dict()
        self.sto_dr = dict()
        for nval in range(1, nmax + 1):
            subs_list = {n: nval}
            csto = stosym.subs(subs_list).evalf()

            sto_str = lambdastr((r, zeta), csto).replace("math", "np")
            self.sto[nval] = eval(sto_str)
            # print('s',sto_str)

    def set_basis(self, basis):
        self.basis = basis

    def eval_v(self, x, y, z):
        r = np.sqrt(x * x + y * y + z * z)
        ang_norm = 1 / np.sqrt(4 * np.pi)
        v = 0.0
        for basis in self.basis:
            for i in range(basis.nbasis):
                v += (
                    ang_norm
                    * basis.contraction_coeff[i]
                    * self.sto[basis.n[i]](r, basis.zeta[i])
                )
        return v

    def eval_v2(self, basis_idx, r):
        r = np.sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2])
        ang_norm = 1 / np.sqrt(4 * np.pi)
        v = 0.0
        basis = self.basis[basis_idx]
        for i in range(basis.nbasis):
            v += (
                ang_norm
                * basis.contraction_coeff[i]
                * self.sto[basis.n[i]](r, basis.zeta[i])
            )
        return v
