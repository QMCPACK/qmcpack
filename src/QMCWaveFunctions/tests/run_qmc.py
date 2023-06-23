
# Minimal QMC loop for validation

# Uses automatic differentiation via the autograd package to
#  compute spatial and parameter derivatives.
# The autodiff is performed in the wavefunction class
from stats import averager
import autograd.numpy as np

# Parameters for run_qmc
# r - numpy array of shape (number of electrons, 3)
# wf - class implementing a wavefunction.
#   It should implement the functions:
#        psi(self, r, VP)
#        dpsi(self, r, VP) - derivative of psi wrt to VP
#        local_energy(self, r, VP)
#        dlocal_energy(self, r, VP) - derivative of local_energy wrt VP
# VP - numpy array of the variational parameters
# nblock - Loop for statistics
# nstep - Loop to gather data
# nsubstep - Loop for decorrelation

def run_qmc(r, wf, VP, nblock=20, nstep=10, nsubstep=10):
    nelectron = r.shape[0]

    # Step size for trial move
    delta = 1.2

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

                    wf_old = wf.psi(r, VP)
                    r_old[:] = r[ne,:]

                    change = delta*(np.random.rand(3)-0.5)

                    r_trial[:] = r_old[:] + change
                    r[ne,:] = r_trial[:]

                    wf_new = wf.psi(r, VP)

                    wf_ratio = (wf_new/wf_old)**2

                    if wf_ratio > np.random.random():
                        naccept += 1
                    else:
                        r[ne,:] = r_old[:]

            eloc = wf.local_energy(r, VP)
            en_ave.add_value(eloc)

            psi_val = wf.psi(r, VP)
            dpsi_val = wf.dpsi(r, VP)
            dpsi_ave.add_value(dpsi_val/psi_val)

            deloc_val = wf.dlocal_energy(r, VP)
            deloc_ave.add_value(deloc_val)

            eloc_dpsi_ave.add_value(eloc * dpsi_val/psi_val)


        en = en_ave.average()
        en_err = en_ave.error()
        ar = naccept/(nstep*nsubstep*nelectron)
        print('block = ',nb, ' energy = ',en,en_err,' acceptance ratio = ',ar, flush=True)


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

