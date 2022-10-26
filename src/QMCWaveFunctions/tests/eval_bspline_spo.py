# 3D B-splines
# Python implementation of SPO spline evaluation
# Reads spline coefficients from HDF file of saved coefficients

import autograd.numpy as np
from autograd import hessian,grad
import h5py


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

def convertPos(r, G):
    # Compact form, but autograd doesn't like item assignment to arrays
    # ru = np.dot(r, G)
    # for i in range(3):
    #     if ru[i] < 0.0:
    #         img = np.floor(ru[i])
    #         ru[i] -= img

    ru = np.dot(r, G)
    # Unroll loop over dimensions and avoid array assignment
    if ru[0] < 0.0:
        img = np.floor(ru[0])
        x = ru[0] - img
    else:
        x = ru[0]

    if ru[1] < 0.0:
        img = np.floor(ru[1])
        y = ru[1] - img
    else:
        y = ru[1]

    if ru[2] < 0.0:
        img = np.floor(ru[2])
        z = ru[2] - img
    else:
        z = ru[2]

    return np.array([x,y,z])



def setup_splines():
    # HDF file with system info
    f2 = h5py.File("diamondC_1x1x1.pwscf.h5")
    prim_vec = f2["supercell/primitive_vectors"]
    G = np.linalg.inv(prim_vec)
    print("primitive vectors = ",np.array(prim_vec))

    # HDF file with spline coefficients.
    # Generate by adding the attribute save_coefs="yes" to determinantset or
    # sposet_builder tag in test_einset.cpp or a QMCPACK run.
    f = h5py.File("einspline.tile_100010001.spin_0.tw_0.l0u8.g40x40x40.h5", "r")

    a = 1.0
    grid1 = Grid1D(40, 0.0, a)
    grids = [grid1, grid1, grid1]

    g = f["spline_0"]
    coeffs1 = g[:,:,:,0]
    coeffs2 = g[:,:,:,1]
    print("coefficients shape = ", g.shape)

    spline1 = Spline3D(grids, coeffs1, G)
    spline2 = Spline3D(grids, coeffs2, G)

    return spline1, spline2


def generate_point_values():
    spline1, spline2 = setup_splines()

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


if __name__ == "__main__":
    generate_point_values()
