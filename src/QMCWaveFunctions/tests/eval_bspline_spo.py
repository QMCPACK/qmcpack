# 3D B-splines
# Python implementation of SPO spline evaluation
# Reads spline coefficients from HDF file of saved coefficients

import autograd.numpy as np
from autograd import hessian,grad
from run_qmc import run_qmc
import read_qmcpack
import h5py
from scipy.interpolate import interp1d


# Integration on a sphere

class AngularIntegration:
    # Use 4 points, corresponds to nrule=2
    def __init__(self):
        self.n = 4
        q = 1/np.sqrt(3)
        self.ang_wt = 0.25
        self.ang_pts = np.array([[ q, q, q],
                                 [ q,-q,-q],
                                 [-q, q,-q],
                                 [-q,-q, q]])
        self.rot_pts = self.ang_pts

    def apply_rot(self, rot):
        self.rot_pts = np.empty_like(self.ang_pts)
        for i in range(self.ang_pts.shape[0]):
            self.rot_pts[i,:] = np.dot(rot, self.ang_pts[i,:])


class Grid1D:
  def __init__(self, num, start, stop):
    self.num = num
    self.start = start
    self.stop = stop

    self.delta = (stop-start)/num
    self.delta_inv = 1.0/self.delta

# Parameters for spline basis functions.
Ad = np.array(
  [ -1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
     3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
    -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
     1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0
  ])

def apply_rotation_to_coeffs(theta, spline1, spline2):
    tmp1 = spline1.coeffs[:,:,:]
    tmp2 = spline2.coeffs[:,:,:]

    spline1.coeffs =  tmp1 * np.cos(theta) + tmp2 * np.sin(theta)
    spline2.coeffs =  -tmp1 * np.sin(theta) + tmp2 * np.cos(theta)


class Spline3D:
    def __init__(self, grids, coeffs, G):
        self.grids = grids
        self.coeffs = coeffs
        self.G = G
        self.evaluate_grad = grad(self.convert_and_evaluate_v, 0)
        self.evaluate_hess = hessian(self.convert_and_evaluate_v, 0)

    def convert_and_evaluate_v(self, x):
        xu = convertPos(x, self.G)
        return self.evaluate_v(xu)

    def evaluate_v(self, x):

        # Unrolled over dimensions to make autograd work
        xx0 = x[0] - self.grids[0].start
        u = xx0 * self.grids[0].delta_inv
        t = u % 1
        i = np.int64(u)
        tp1 = np.array([t*t*t, t*t, t, 1.0])

        xx1 = x[1] - self.grids[1].start
        u1 = xx1 * self.grids[1].delta_inv
        t = u1 % 1
        j = np.int64(u1)
        tp2 = np.array([t*t*t, t*t, t, 1.0])

        xx2 = x[2] - self.grids[2].start
        u2 = xx2 * self.grids[2].delta_inv
        t = u2 % 1
        k = np.int64(u2)
        tp3 = np.array([t*t*t, t*t, t, 1.0])

        a = np.array([np.dot(Ad[0:4], tp1),
                      np.dot(Ad[4:8], tp1),
                      np.dot(Ad[8:12], tp1),
                      np.dot(Ad[12:16], tp1)])
        b = np.array([np.dot(Ad[0:4], tp2),
                      np.dot(Ad[4:8], tp2),
                      np.dot(Ad[8:12], tp2),
                      np.dot(Ad[12:16], tp2)])
        c = np.array([np.dot(Ad[0:4], tp3),
                      np.dot(Ad[4:8], tp3),
                      np.dot(Ad[8:12], tp3),
                      np.dot(Ad[12:16], tp3)])

        val = 0.0
        for idx in range(4):
          for jdx in range(4):
            pre00 = a[idx]*b[jdx]
            sum00 = c[0]*self.coeffs[i+idx, j+jdx, k  ] + \
                    c[1]*self.coeffs[i+idx, j+jdx, k+1] + \
                    c[2]*self.coeffs[i+idx, j+jdx, k+2] + \
                    c[3]*self.coeffs[i+idx, j+jdx, k+3]
            val += pre00 * sum00

        return val

# Taken from src/QMCWaveFunctions/BsplineFactory/SplineR2R.h convertPos
def convertPos(r, G):
    ru = np.dot(r, G)
    img = np.floor(ru)
    return ru - img

def apply_bc(r, lat, G):
    ru = np.dot(r, G)
    img = np.floor(ru + 0.5)
    x = ru - img
    return np.dot(x,lat)


def setup_splines(pwscf_file, coeff_file):
    # HDF file with system info
    f2 = h5py.File(pwscf_file)
    prim_vec = f2["supercell/primitive_vectors"]
    G = np.linalg.inv(prim_vec)
    print("primitive vectors = ",np.array(prim_vec))

    # HDF file with spline coefficients.
    # Generate by adding the attribute save_coefs="yes" to determinantset or
    # sposet_builder tag in test_einset.cpp or a QMCPACK run.
    f = h5py.File(coeff_file, "r")

    g = f["spline_0"]

    # Splines are defined from 0 to 1.
    # The conversion from cell coordinates to the unit cube happens in convertPos
    grid1 = Grid1D(g.shape[0] - 3, 0.0, 1.0)
    grid2 = Grid1D(g.shape[1] - 3, 0.0, 1.0)
    grid3 = Grid1D(g.shape[2] - 3, 0.0, 1.0)
    grids = [grid1, grid2, grid3]

    coeffs1 = g[:,:,:,0]
    coeffs2 = g[:,:,:,1]
    print("coefficients shape = ", g.shape)

    spline1 = Spline3D(grids, coeffs1, G)
    spline2 = Spline3D(grids, coeffs2, G)

    return spline1, spline2, prim_vec, G

class Wavefunction_hcpBe:
    def __init__(self, spline1, spline2, PrimVec, G):
        self.spline1 = spline1
        self.spline2 = spline2
        self.prim_vec = PrimVec
        self.G = G

        self.hess0 = hessian(self.psi_internal, 0)
        self.hess1 = hessian(self.psi_internal, 1)

        self.dpsi = grad(self.psi, 1)
        self.dkinetic_energy = grad(self.kinetic_energy, 1)
        self.dlocal_energy = grad(self.local_energy, 1)

        # Derivatives of a single orbital
        self.dorb = grad(self.orb, 1)
        self.orb_hess0 = hessian(self.orb, 0)
        self.grad0 = grad(self.orb, 0)
        self.dlap0 = grad(self.lap0, 1)

        self.ang = AngularIntegration()

    def mag(self, r):
        return np.sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2])

    def orb(self, r, theta):
        o1 = self.spline1.convert_and_evaluate_v(r)
        o2 = self.spline2.convert_and_evaluate_v(r)
        return o1 * np.cos(theta) + o2 * np.sin(theta)

    def psi_internal(self, r1, r2, VP):
        theta1 = VP[0]
        theta2 = VP[1]
        o1 = self.orb(r1, theta1)
        o2 = self.orb(r2, theta2)
        return o1*o2

    def psi(self, r, VP):
        r1 = r[0,:]
        r2 = r[1,:]
        return self.psi_internal(r1, r2, VP)

    def kinetic_energy(self, r, VP):
        r1 = r[0,:]
        r2 = r[1,:]
        psi_val = self.psi_internal(r1, r2, VP)
        h0 = np.sum(np.diag(self.hess0(r1, r2, VP)))
        h1 = np.sum(np.diag(self.hess1(r1, r2, VP)))

        # only care about the parameter derivative for now, so the potential term doesn't matter
        h = -0.5*(h0+h1)/psi_val
        return h

    def local_energy(self, r, VP):
        r1 = r[0,:]
        r2 = r[1,:]
        psi_val = self.psi_internal(r1, r2, VP)
        h0 = np.sum(np.diag(self.hess0(r1, r2, VP)))
        h1 = np.sum(np.diag(self.hess1(r1, r2, VP)))

        pp_en1 = self.pp_energy(r1)
        pp_en2 = self.pp_energy(r2)

        nlpp_en = self.nlpp_en(r, VP)

        h = -0.5*(h0+h1)/psi_val + pp_en1 + pp_en2 + nlpp_en
        return h

    # Laplacian for a single orbital
    def lap0(self, r1, VP):
        psi_val = self.orb(r1, VP)
        h0 = np.sum(np.diag(self.orb_hess0(r1, VP)))
        g0 = self.grad0(r1,VP)/psi_val
        return -0.5*(h0/psi_val - np.dot(g0, g0))


    def set_pp_func(self, pp_data):
        # Local/nonlocal channels specific to Be
        (grid,vals) = pp_data[('0','s')]
        (grid_loc,vals_loc) = pp_data[('0','p')]
        Zeff = 2

        # From ECPComponentBuilder.2.cpp doBreakUp
        # It appears the nonlocal part is the nonlocal input channels minus the local channel

        # Spline V*r/Zeff
        self.pp_func = interp1d(grid, vals_loc/Zeff, kind='cubic',bounds_error=False,fill_value=(0.0,0.0))

        # Spline V
        self.nlpp_func = interp1d(grid, (vals-vals_loc)/(grid+1e-20), kind='cubic',bounds_error=False,fill_value=(0.0,1.0))

    # Potential from local part of the pseudopotential
    def pp_energy(self, r1):
        Zeff = 2
        r1p = apply_bc(r1, self.prim_vec, self.G)
        r = self.mag(r1p)

        # Spline V*r/Zeff
        en = Zeff*self.pp_func(r)/r
        return en

    # Non-local potential
    def nlpp_pot(self, r1):
        r1p = apply_bc(r1, self.prim_vec, self.G)
        r = self.mag(r1p)

        # Spline V
        en = self.nlpp_func(r)
        return en

    # NLPP energy integrated over the sphere
    def nlpp_en(self, r, VP):
        nl = 0.0
        rt = np.zeros(r.shape)
        p0 = self.psi(r, VP)
        for i in range(2):
            ri = r[i,:]
            rt[:,:] = r[:,:]
            rm = self.mag(ri)
            proj_00 = 0.0
            pot = self.nlpp_pot(ri)
            for j in range(self.ang.n):
                rt[i,:] = rm * self.ang.rot_pts[j,:]
                proj_00 += self.ang.ang_wt * self.psi(rt, VP)/p0 * pot
            nl += proj_00

        return nl




def generate_point_values_diamondC():
    pwscf_file = "diamondC_1x1x1.pwscf.h5"
    coeff_file = "einspline.tile_100010001.spin_0.tw_0.l0u8.g40x40x40.h5"
    spline1, spline2, _, _ = setup_splines(pwscf_file, coeff_file)

    r1 = np.array([0.0, 0.0, 0.0])
    r2 = np.array([0.0, 1.0, 0.0])

    o1 = spline1.convert_and_evaluate_v(r1)
    print("orbital 1 for particle 1 = ", o1)

    o2 = spline1.convert_and_evaluate_v(r2)
    print("orbital 1 for particle 2 = ", o2)

    o11 = spline2.convert_and_evaluate_v(r1)
    print("orbital 2 for particle 1 = ", o11)

    o21 = spline2.convert_and_evaluate_v(r2)
    print("orbital 2 for particle 2 = ", o21)

    go1 = spline1.evaluate_grad(r2)
    print("gradient for robital 1, particle 2 = ", go1)

    ho1 = spline1.evaluate_hess(r1)
    lap = ho1[0,0] + ho1[1,1] + ho1[2,2]
    print("laplacian for particle 1 = ", lap)

    ho2 = spline1.evaluate_hess(r2)
    lap = ho2[0,0] + ho2[1,1] + ho2[2,2]
    print("lap for particle 2 = ", lap)
    print("hessian for particle 2", ho2)

    ho2 = spline2.evaluate_hess(r2)
    lap = ho2[0,0] + ho2[1,1] + ho2[2,2]
    print("lap for orbital 2, particle 2 = ", lap)
    print("hessian for orbital 2, particle 2", ho2)


def generate_point_values_hcpBe():
    pwscf_file = "hcpBe.pwscf.h5"
    coeff_file = "einspline.tile_100010001.spin_0.tw_0.l0u2.g40x40x68.h5"
    spline1, spline2, prim_vec, G = setup_splines(pwscf_file, coeff_file)

    r1 = np.array([0.0, 0.0, 0.0])
    o1 = spline1.convert_and_evaluate_v(r1)
    print("orbital 1 for particle 1 = ", o1)
    o2 = spline2.convert_and_evaluate_v(r1)
    print("orbital 2 for particle 1 = ", o2)

    ho1 = spline1.evaluate_hess(r1)
    lap = ho1[0,0] + ho1[1,1] + ho1[2,2]
    print("laplacian for particle 1 = ", lap)

    apply_rotation_to_coeffs(0.1, spline1, spline2)
    print("After rotation")
    o1 = spline1.convert_and_evaluate_v(r1)
    print("orbital 1 for particle 1 = ", o1)
    o2 = spline2.convert_and_evaluate_v(r1)
    print("orbital 2 for particle 1 = ", o2)

    ho1 = spline1.evaluate_hess(r1)
    lap = ho1[0,0] + ho1[1,1] + ho1[2,2]
    print("laplacian for particle 1 = ", lap)

    VP = np.array([0.0])
    wf = Wavefunction_hcpBe(spline1, spline2, prim_vec, G)
    psi = wf.orb(r1, VP[0])
    dp = wf.dorb(r1, VP[0])
    print("dpsi for r1 = ", dp)
    print("d log(psi) for r1 = ", dp/psi)

    lap0 = wf.lap0(r1, VP[0])
    print("lap 0 = ",lap0)
    dlap0 = wf.dlap0(r1, VP[0])
    print("dlap 0 = ",dlap0)

def generate_point_values_hcpBe_NLPP():
    pwscf_file = "hcpBe.pwscf.h5"
    coeff_file = "einspline.tile_100010001.spin_0.tw_0.l0u2.g40x40x68.h5"
    spline1, spline2, prim_vec, G = setup_splines(pwscf_file, coeff_file)

    pp_data = read_qmcpack.read_pp_file("Be.BFD.xml")

    r = np.array([[0.1, 0.2, 0.3],
                  [1.0, 1.0, 1.0]])
    VP = np.array([0.0, 0.0])

    wf = Wavefunction_hcpBe(spline1, spline2, prim_vec, G)
    wf.set_pp_func(pp_data)
    psi = wf.psi(r, VP)
    print("log(psi) = ",np.log(psi))

    ke = wf.kinetic_energy(r, VP)
    print("ke = ",ke)

    eloc = wf.local_energy(r, VP)
    print("loc_e = ",eloc)

    dp = wf.dpsi(r, VP)
    print("d log(psi) = ", dp/psi)
    dke = wf.dkinetic_energy(r, VP)
    print("d kinetic energy (no NLPP) = ", dke)
    de = wf.dlocal_energy(r, VP)
    print("d energy (with NLPP) = ", de)


def run_qmc_parameter_derivative_hcpBe():
    pwscf_file = "hcpBe.pwscf.h5"
    coeff_file = "einspline.tile_100010001.spin_0.tw_0.l0u2.g40x40x68.h5"
    spline1, spline2, prim_vec, G = setup_splines(pwscf_file, coeff_file)
    apply_rotation_to_coeffs(0.1, spline1, spline2)
    VP = np.array([0.0, 0.0])
    wf = Wavefunction_hcpBe(spline1, spline2, prim_vec, G)
    #r = np.array([[0.0, 0.0, 0.0],
    #              [1.0, 2.1, 0.1]])
    r = np.array([[0.0, 0.0, 0.0],
                  [0.0, 0.0, 0.0]])

    #de = wf.dlocal_energy(r, VP)
    #print('de = ',de)
    run_qmc(r, wf, VP, nstep=40, nblock=10)


if __name__ == "__main__":
    #generate_point_values_diamondC()
    #generate_point_values_hcpBe()
    #run_qmc_parameter_derivative_hcpBe()
    generate_point_values_hcpBe_NLPP()
