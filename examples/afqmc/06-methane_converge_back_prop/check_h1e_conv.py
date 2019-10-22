# Analyse the AFQMC back propagated RDM.
import glob
import h5py
import numpy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy.stats
from afqmctools.analysis import rdm

def compute_one_body(filename, hcore, skip, group):
    energies = []
    # weights can be ignored if not using free projection.
    dm, weights = rdm.extract_rdm(filename, dm_name=group)
    for d in dm[skip:]:
        e1b = numpy.einsum('ij,ij->', hcore, d[0]+d[1])
        energies.append(e1b.real)
    av = numpy.mean(energies)
    err = scipy.stats.sem(energies)
    return av, err

def plot_convergence(filename):
    with h5py.File('afqmc.h5', 'r') as fh5:
        hcore = fh5['Hamiltonian/hcore'][:].view(numpy.complex128)[:,:,0]
    energies = []
    errs = []
    tau_bps = []
    filename = filename
    metadata = rdm.get_metadata(filename)
    taus = metadata['num_bp']
    for i, t in enumerate(taus):
        # skip the first block for equilibration
        skip = 1
        nelec = metadata['nalpha'] + metadata['nbeta']
        # We can extract the averaged 1RDM.
        group = 'Observables/BackPropagated/FullOneRDM/Average_{}'.format(str(i))
        rdm_av, rdm_errs = rdm.get_one_rdm_av(filename, skip, dm_name=group)
        nelec_rdm = (rdm_av[0].trace() + rdm_av[1].trace()).real
        assert(nelec_rdm-nelec < 1e-12)
        # Often it is simpler to compute error bars if we first contract the 1RDM.
        # For example, we can compute the one-body energy from the averaged RDM as.
        e1b = numpy.einsum('ij,sij->', hcore, rdm_av).real
        # Or we can instead compute the average of the one-body energies.
        e1b_series, err = compute_one_body(filename, hcore, 1, group)
        energies.append(e1b_series)
        errs.append(err)
        # get back propagation time
        tau_bps.append(t*metadata['dt'])
        assert(e1b-e1b_series < 1e-12)

    # Finally plot the one-body energy and check the estimator is converged with
    # respect to back propagation time.
    plt.errorbar(tau_bps, energies, yerr=errs, fmt='o')
    plt.xlabel(r'$\tau_{BP}$')
    plt.ylabel(r'$E_{1B}$ (Ha)')
    plt.savefig('h1e_conv.pdf', fmt='pdf', bbox_inches='tight')

if __name__ == '__main__':
    plot_convergence('qmc.s000.stat.h5')
