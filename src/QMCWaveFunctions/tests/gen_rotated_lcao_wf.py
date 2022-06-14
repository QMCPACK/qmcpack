
# Helium atom with a combination of two orbitals and simple jastrow factor

# Uses automatic differentiation via the autograd package to
#  compute spatial and parameter derivatives
import autograd.numpy as np
from autograd import hessian,grad

# Values used in test_RotatedSPOs_LCAO.cpp

def mag(r):
    return np.sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2])

# normalized STO's correspond to the 'normalized="no"' part of the input
#     <atomicBasisSet type="STO" elementType="He" normalized="no">


def sto_norm1(zeta):
  return 2*np.sqrt(zeta**3)

def sto_norm2(zeta):
  return 2*np.sqrt(3)*np.sqrt(zeta**5)/3


def orb1(R):
    r = mag(R)
    Z = 2.0
    y00 = 1/np.sqrt(4 * np.pi)
    snorm1 = sto_norm1(Z)
    return y00 * snorm1 * np.exp(-Z*r)

def orb2(R):
    r = mag(R)
    zeta = 1.0
    y00 = 1/np.sqrt(4*np.pi)
    snorm2 = sto_norm2(zeta)
    return snorm2* y00 * r* np.exp(-zeta*r)

def rot_orb(R, theta):
    return orb1(R) * np.cos(theta) + orb2(R) * np.sin(theta)


def jastrow(r12, B):
    A = 0.5
    return np.exp(A*r12/(1.0 + B*r12) - A/B)

def psi_no_jastrow(r1, r2, VP):
    theta1 = VP[0]
    theta2 = VP[1]
    o1 = rot_orb(r1,theta1)
    o2 = rot_orb(r2,theta2)
    return o1*o2

def psi_with_jastrow(r1, r2, VP):
    theta1 = VP[0]
    theta2 = VP[1]
    B = VP[2]
    o1 = rot_orb(r1,theta1)
    o2 = rot_orb(r2,theta2)
    r12 = r2 - r1
    j = jastrow(r12, B)
    return o1*o2*j

def psi(r1, r2, VP):
    global use_jastrow
    theta1 = VP[0]
    theta2 = VP[1]
    j = 1.0
    if use_jastrow:
        B = VP[2]
        r12 = mag(r2 - r1)
        j = jastrow(r12, B)

    o1 = rot_orb(r1,theta1)
    o2 = rot_orb(r2,theta2)
    return o1*o2*j

def log_psi(r1, r2, B):
    return np.log(psi(r1, r2, B))


# Spatial derivatives
hess0 = hessian(psi, 0)
hess1 = hessian(psi, 1)


hess_log_0 = hessian(log_psi, 0)
hess_log_1 = hessian(log_psi, 1)

grad0 = grad(psi, 0)
grad1 = grad(psi, 1)

# Derivative wrt parameters
dpsi = grad(psi, 2)

def lap0(r1, r2, VP):
    h0 = np.sum(np.diag(hess_log_0(r1, r2, VP)))
    return h0

def lap1(r1, r2, VP):
    h1 = np.sum(np.diag(hess_log_1(r1, r2, VP)))
    return h1

def lap(r1, r2, VP):
    h0 = np.sum(np.diag(hess0(r1, r2, VP)))
    h1 = np.sum(np.diag(hess1(r1, r2, VP)))
    return h0 + h1

def en_pot(r1, r2):
    r1_mag = mag(r1)
    r2_mag = mag(r2)
    Z = 2.0
    return -Z/r1_mag - Z/r2_mag

def ee_pot(r1, r2):
    r12 = r2 - r1
    r12_mag = mag(r12)
    return 1.0/r12_mag


def local_energy(r1, r2, VP):
    pot = en_pot(r1, r2) + ee_pot(r1, r2)
    psi_val = psi(r1, r2, VP)
    lapl = lap(r1, r2, VP)

    h = -0.5*lapl/psi_val + pot
    return h

dlocal_energy = grad(local_energy, 2)



def print_wf_values(theta1=0.0, theta2=0.0,  use_j=False, B=0.0):
    global use_jastrow

    # Adjust numpy output so arrays are printed with higher precision
    float_formatter = "{:.15g}".format
    np.set_printoptions(formatter={'float_kind':float_formatter})

    if use_j:
        use_jastrow = True
        VP = np.array([theta1, theta2, B])
        print("Values for theta = ",theta1,theta2," and jastrow B = ",B)
    else:
        use_jastrow = False
        VP = np.array([theta1, theta2])
        print("Values for theta = ",theta1,theta2," and no jastrow")



    r1 = np.array([1.0, 2.0, 3.0])
    r2 = np.array([0.0, 1.1, 2.2])

    psi_val = psi(r1, r2, VP)
    print("  wf = ",psi_val," log wf = ",np.log(np.abs(psi_val)))

    g0 = grad0(r1, r2, VP)/psi_val
    print("  grad/psi for particle 0 = ",g0[0],g0[1],g0[2])

    # Using the laplacian of log psi to match internal QMCPACK values
    lap_0 = lap0(r1, r2, VP)
    print(" laplacian of log psi for particle 0 = ",lap_0)

    lap_1 = lap1(r1, r2, VP)
    print(" laplacian for log psi particle 1 = ",lap_1)

    eloc = local_energy(r1, r2, VP)
    print("  local energy = ",eloc)

    dp = dpsi(r1, r2, VP)
    print("  parameter derivative of log psi = ",dp / psi_val)

    deloc = dlocal_energy(r1, r2, VP)
    print("  parameter derivative of local energy = ",deloc)

    print("")

def print_point_values():
    r1 = np.array([1.0, 2.0, 3.0])
    r2 = np.array([0.0, 1.1, 2.2])

    print_wf_values(theta1=0.1, theta2=0.2)

    print_wf_values(theta1=0.0, theta2=0.0)

    print_wf_values(theta1=0.0, theta2=0.0, use_j=True, B=0.1)


if __name__=='__main__':
    print_point_values()

