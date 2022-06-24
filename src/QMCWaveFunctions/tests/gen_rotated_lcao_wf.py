
# Helium atom with a combination of two orbitals and simple jastrow factor

# Uses automatic differentiation via the autograd package to
#  compute spatial and parameter derivatives
import autograd.numpy as np
from autograd import hessian,grad
from stats import averager

# Point values used in test_RotatedSPOs_LCAO.cpp
# QMC values used to validate tests/molecules/He_param/He_orb_rot_param_grad_legacy


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


# Generate the wavefunction values for a single set of electron positions
# used in test_RotatedSPOs_LCAO.cpp

def print_point_values():
    r1 = np.array([1.0, 2.0, 3.0])
    r2 = np.array([0.0, 1.1, 2.2])

    print_wf_values(theta1=0.1, theta2=0.2)

    print_wf_values(theta1=0.0, theta2=0.0)

    print_wf_values(theta1=0.0, theta2=0.0, use_j=True, B=0.1)

def run_qmc(VP):
    r1 = np.array([1.0, 2.0, 3.0])
    r2 = np.array([0.0, 1.1, 2.2])

    # Outer loop for statistics
    nblock = 20

    # Loop to gather data
    nstep = 100

    # Loop for decorrelation
    nsubstep = 10

    nelectron = 2

    # Step size for trial move
    delta = 1.2

    r = np.zeros((2,3))
    r[0,:] = r1
    r[1,:] = r2

    r_old = np.zeros(3)
    r_trial = np.zeros(3)

    block_dp_ave = averager()
    total_en_ave = averager()

    # Outer loop for statistics
    for nb in range(nblock):

        naccept = 0
        en_ave = averager()
        dpsi_ave = averager()
        deloc_ave = averager()
        eloc_dpsi_ave = averager()

        # Loop to gather data
        for ns in range(nstep):
            # Loop for decorrelation
            for nss in range(nsubstep):
                # Loop over electrons
                for ne in range(nelectron):

                    wf_old = psi(r[0,:], r[1,:], VP)
                    r_old[:] = r[ne,:]

                    change = delta*(np.random.rand(3)-0.5)

                    r_trial[:] = r_old[:] + change
                    r[ne,:] = r_trial[:]

                    wf_new = psi(r[0,:], r[1,:], VP)

                    wf_ratio = (wf_new/wf_old)**2

                    if wf_ratio > np.random.random():
                        naccept += 1
                    else:
                        r[ne,:] = r_old[:]

            eloc = local_energy(r[0,:], r[1,:], VP)
            en_ave.add_value(eloc)

            psi_val = psi(r[0,:], r[1,:], VP)
            dpsi_val = dpsi(r[0,:], r[1,:], VP)
            dpsi_ave.add_value(dpsi_val/psi_val)

            deloc_val = dlocal_energy(r[0,:], r[1,:], VP)
            deloc_ave.add_value(deloc_val)

            eloc_dpsi_ave.add_value(eloc * dpsi_val/psi_val)


        en = en_ave.average()
        en_err = en_ave.error()
        ar = naccept/(nstep*nsubstep*nelectron)
        print('block = ',nb, ' energy = ',en,en_err,' acceptance ratio = ',ar)

        dp = dpsi_ave.average()
        dp_err = dpsi_ave.error()

        deloc = deloc_ave.average()
        deloc_err = deloc_ave.error()

        eloc_dpsi = eloc_dpsi_ave.average()
        eloc_dpsi_err = eloc_dpsi_ave.error()

        # For the parameter derivative formula, see
        # https://github.com/QMCPACK/qmc_algorithms/blob/master/Variational/Parameter_Optimization.ipynb

        dg = deloc + 2*eloc_dpsi - 2*en*dp
        block_dp_ave.add_value(dg)

    dg = block_dp_ave.average()
    dg_err = block_dp_ave.error()
    print('parameter values = ',VP)
    print('parameter derivatives = ',dg)
    print('parameter derivative errors = ',dg_err)


def run_qmc_parameter_derivatives():
    global use_jastrow
    use_jastrow=False
    theta1 = 0.1
    theta2 = 0.1
    VP = np.array([theta1, theta2])

    run_qmc(VP)

# Some results from run_qmc_parameter_derivatives

# Run took about 10 minutes on laptop
# nblock=20, nstep=1000, nsubstep=10
# parameter values =  [0.1 0.1]
# parameter derivatives =  [-0.20164722 -0.18347461]
# parameter derivative errors =  [0.01201481 0.01314164]

# Run took about 40 minutes on laptop
# nblock=40, nstep=2000, nsubstep=10
# parameter values =  [0.1 0.1]
# parameter derivatives =  [-0.2204924  -0.21471184]
# parameter derivative errors =  [0.00493837 0.00571082]


if __name__=='__main__':
    #print_point_values()
    run_qmc_parameter_derivatives()

