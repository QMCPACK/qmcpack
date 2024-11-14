
# Compute wavefunction values and parameter derivatives
# for a wavefunction with STO Be orbitals, two determinants, and orbital rotation

import autograd.numpy as np
from autograd import hessian, grad

from run_qmc import run_qmc
import read_qmcpack
from slater_orbitals import STO

import scipy.linalg

# From construct_rot.py
def construct_antisym_ex(p):
    return np.array(
        [
            [0, -p[0], -p[1], -p[2], -p[3], -p[4], -p[5]],
            [p[0], 0, -p[6], -p[7], -p[8], -p[9], -p[10]],
            [p[1], p[6], 0, -p[11], -p[12], -p[13], -p[14]],
            [p[2], p[7], p[11], 0, 0, 0, 0],
            [p[3], p[8], p[12], 0, 0, 0, 0],
            [p[4], p[9], p[13], 0, 0, 0, 0],
            [p[5], p[10], p[14], 0, 0, 0, 0],
        ]
    )


# 2x2 determinant between two states
def det2_ex(phi, i, j):
    return phi[i, 0] * phi[j, 1] - phi[j, 0] * phi[i, 1]


def mat_exp(m):
    # Simple approximation good enough for derivatives at zero rotation
    # Might only need to go up to the linear term
    return np.eye(m.shape[0]) + m + np.dot(m, m) / 2


class Wavefunction_Be_STO:
    def __init__(self, basis, mo_coeff):
        self.sto = STO()
        self.sto.set_basis(basis)
        self.mo_coeff = mo_coeff
        self.hess = list()
        for i in range(4):
            self.hess.append(hessian(self.psi_internal, i))

        self.nmo = 7

        self.rot_param_size = 15  # size of p in construct_antisym_ex

        self.dpsi = grad(self.psi, 1)

        self.dlocal_energy = grad(self.local_energy, 1)

    def mag(self, r):
        return np.sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2])

    def psi_internal(self, r1, r2, r3, r4, VP):

        param = VP[1 : self.rot_param_size + 1]

        rot = construct_antisym_ex(param)
        # Can use this line if not doing autodiff
        # rot_mat = scipy.linalg.expm(-rot)
        rot_mat = mat_exp(-rot)
        rot_mo = np.dot(rot_mat, self.mo_coeff)

        mo_size = self.mo_coeff.shape[0]
        phi0_1 = [self.sto.eval_v2(i, r1) for i in range(mo_size)]
        phi0_2 = [self.sto.eval_v2(i, r2) for i in range(mo_size)]
        phi0_a = np.array([phi0_1, phi0_2])
        phi0 = np.dot(rot_mo, phi0_a.T)

        phi1_1 = [self.sto.eval_v2(i, r3) for i in range(mo_size)]
        phi1_2 = [self.sto.eval_v2(i, r4) for i in range(mo_size)]
        phi1_a = np.array([phi1_1, phi1_2])
        phi1 = np.dot(rot_mo, phi1_a.T)

        d1 = det2_ex(phi0, 0, 1)
        d2 = det2_ex(phi1, 0, 1)
        c1 = 1.0

        d3 = det2_ex(phi0, 0, 2)
        d4 = det2_ex(phi1, 0, 2)
        c2 = VP[0]

        return c1 * d1 * d2 + c2 * d3 * d4

    def psi(self, r, VP):
        r1 = r[0, :]
        r2 = r[1, :]
        r3 = r[2, :]
        r4 = r[3, :]
        return self.psi_internal(r1, r2, r3, r4, VP)

    def dpsi(self, r, VP):
        r1 = r[0, :]
        r2 = r[1, :]
        r3 = r[2, :]
        r4 = r[3, :]
        return self.dpsi_internal(r1, r2, r3, r4, VP)

    def en_pot(self, r):
        Z = 4.0
        total = 0.0
        for i in range(r.shape[0]):
            r_mag = self.mag(r[i, :])
            total += -Z / r_mag
        return total

    def ee_pot(self, r):
        total = 0.0
        for i in range(r.shape[0]):
            for j in range(i):
                rij = r[j, :] - r[i, :]
                rij_mag = self.mag(rij)
                total += 1.0 / rij_mag
        return total

    def lap(self, r1, r2, r3, r4, VP):
        h = 0.0
        for i in range(4):
            h += np.sum(np.diag(self.hess[i](r1, r2, r3, r4, VP)))
        return h

    def local_energy(self, r, VP):
        r1 = r[0, :]
        r2 = r[1, :]
        r3 = r[2, :]
        r4 = r[3, :]
        pot = self.en_pot(r) + self.ee_pot(r)
        psi_val = self.psi_internal(r1, r2, r3, r4, VP)
        lapl = self.lap(r1, r2, r3, r4, VP)
        h = -0.5 * lapl / psi_val + pot
        return h

    def dlocal_energy(self, r, VP):
        r1 = r[0, :]
        r2 = r[1, :]
        r3 = r[2, :]
        r4 = r[3, :]
        pot = self.en_pot(r) + self.ee_pot(r)
        psi_val = self.psi_internal(r1, r2, r3, r4, VP)
        lapl = self.lap(r1, r2, r3, r4, VP)
        h = -0.5 * lapl / psi_val + pot
        return h


# Create reference values for
# "Rotated LCAO Be two determinant" in test_RotatedSPOs_LCAO.cpp
def gen_point_derivatives():
    # only uses the basis set
    fname = "rot_multi_2det_Be_STO.wfnoj.xml"
    basis, mo = read_qmcpack.parse_qmc_wf(fname, ["Be"])
    wf = Wavefunction_Be_STO(basis["Be"], mo)

    r = np.array([[0.7, 2.0, 3.0], [1.2, 1.5, 0.5], [1.5, 1.6, 1.5], [0.7, 1.0, 1.2]])

    VP = np.zeros(wf.rot_param_size + 1)
    VP[0] = 0.1
    p = wf.psi(r, VP)
    print("psi = ", p)
    print("log psi = ", np.log(abs(p)))

    dp = wf.dpsi(r, VP)
    print("dpsi = ", dp)

    print("dlogpsi = ", dp / p)
    dlogpsi = dp / p
    for i in range(dlogpsi.shape[0]):
        print("    CHECK(dlogpsi[{}] == ValueApprox({}));".format(i, dlogpsi[i]))

    en = wf.local_energy(r, VP)
    print("en = ", en)

    den = wf.dlocal_energy(r, VP)
    print("den = ", den)

    for i in range(den.shape[0]):
        print("   CHECK(dhpsioverpsi[{}] == ValueApprox({}));".format(i, den[i]))


def run():
    fname = "rot_multi_2det_Be_STO.wfnoj.xml"
    basis, mo = read_qmcpack.parse_qmc_wf(fname, ["Be"])
    wf = Wavefunction_Be_STO(basis["Be"], mo)

    r = np.array([[0.7, 2.0, 3.0], [1.2, 1.5, 0.5], [1.5, 1.6, 1.5], [0.7, 1.0, 1.2]])
    VP = np.zeros(wf.rot_param_size + 1)
    VP[0] = 0.1
    run_qmc(r, wf, VP, nstep=10, nblock=10)


if __name__ == "__main__":
    gen_point_derivatives()
    #run()
