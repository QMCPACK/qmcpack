"""Write mean-field trial wavefunctions to file."""
import numpy
import h5py
import time
import sys
from afqmctools.wavefunction.mol import write_qmcpack_wfn

def write_wfn_pbc(scf_data, ortho_ao, filename, rediag=True, verbose=False,
                  energy_sort=False):
    """Generate QMCPACK trial wavefunction for PBC simulation.

    Parameters
    ----------
    scf_data : dict
        Dictionary containing scf data extracted from pyscf checkpoint file.
    ortho_ao : bool
        Whether we are working in orthogonalised AO basis or not.
    filename : string
        HDF5 file path to store wavefunction to.
    rediag : bool
        Whether to rediagonalise Fock matrix to compute MO coeffs in
        orthoagonlised AO basis. Default: True.
    verbose : bool
        Print additional information. Default: False.
    energy_sort : bool
        Whether to sort MOs globally by energy. Default: Wavefunction will be
        block diagonal in k-space.

    Returns
    -------
    wfn : :class:`numpy.ndarray'
        Wavefunction as numpy array. Format depends on wavefunction.
    """
    # Single determinant for the moment.
    ghf = False
    cell = scf_data['cell']
    C = scf_data['mo_coeff']
    X = scf_data['X']
    mo_occ = numpy.array(scf_data['Xocc'])
    fock = scf_data['fock']
    nmo_pk = scf_data['nmo_pk']
    nmo_tot = numpy.sum(nmo_pk)
    kpts = scf_data['kpts']
    nkpts = len(kpts)
    nalpha = cell.nelec[0] * nkpts
    nbeta = cell.nelec[1] * nkpts
    nelec = nalpha + nbeta
    uhf = scf_data['isUHF']
    if verbose:
        print(" # Generating QMCPACK trial wavefunction.")
        print(" # Shape of supercell wavefunction "
              "({:d},{:d})".format(nmo_tot, nelec if uhf else nalpha))
        if uhf:
            print(" # UHF-like trial wavefunction.")
        else:
            print(" # RHF-like trial wavefunction.")
    # For RHF only nalpha entries will be filled.
    wfn = numpy.zeros((1,nmo_tot,nalpha+nbeta), dtype=numpy.complex128)
    eigs_a = []
    eigs_b = []
    occ_b = numpy.zeros(1, dtype=numpy.int64)
    ks = []
    if ortho_ao:
        row = 0
        col = 0
        col_b = nalpha
        for k in range(nkpts):
            if verbose:
                print(" # Generating trial wavefunction for momentum: "
                      "{:d}".format(k))
            if rediag:
                start = time.time()
                if uhf:
                    occ_a = mo_occ[0,k] > 0
                    ea, wfn_a = rediag_fock(fock[0,k], X[k][:,:nmo_pk[k]], occ_a)
                    eigs_a.append(ea)
                    occ_b = mo_occ[1,k] > 0
                    eb, wfn_b = rediag_fock(fock[1,k], X[k][:,:nmo_pk[k]], occ_b)
                    eigs_b.append(eb)
                    ks.append([k]*len(eb))
                else:
                    occ_a = mo_occ[k] > 0
                    ea, wfn_a = rediag_fock(fock[k], X[k][:,:nmo_pk[k]], occ_a)
                    eigs_a.append(ea)
                    ks.append([k]*len(ea))
                if verbose:
                    print(" # Time to rediagonalise fock: "
                          "  {:13.8e} s".format(time.time()-start))
            else:
                # Potentially dangerous, I'm not sure if this is valid.
                print("rediag = False not implemented yet.")
                sys.exit()
                # Xinv = scipy.linalg.pinv(X[ik,:,:nmo_pk[ik]])
                # if uhf:
                    # # We are assuming C matrix is energy ordered.
                    # occ_a = mo_occ[0] > 0
                    # wfn_a = numpy.dot(Xinv, C[k,0])[:,occ_a]
                    # occ_b = mo_occ[1] > 0
                    # wfn_b = numpy.dot(Xinv, C[k,1])[:,occ_b]
                # else:
                    # occ_a = mo_occ[0] > 0
                    # wfn_a = numpy.dot(Xinv, C[k])[:,occ_a]
            row_end = row + nmo_pk[k]
            col_end = col + sum(occ_a)
            wfn[0,row:row_end,col:col_end] = wfn_a
            if uhf:
                col_end = col_b + sum(occ_b)
                wfn[0,row:row_end,col:col_end] = wfn_b
                col_b += sum(occ_b)
            row += nmo_pk[k]
            col += sum(occ_a)
        # Sort columns by energy.
        eigs_a = numpy.array(eigs_a).ravel()
        ks = numpy.array(ks).ravel()
        col_sort = list(numpy.argsort(eigs_a))
        if verbose:
            print(" # Recomputed MO energies.")
            print(" # {:5s}  {:>17s}".format("kpoint", "e_mo"))
            for e, ix in zip(eigs_a[col_sort], numpy.array(ks)[col_sort]):
                print(" #    {:3d}    {: 13.8e}".format(ix, e))
        if uhf:
            eigs_b = numpy.array(eigs_b).ravel()
            col_sort_b = list(nalpha + numpy.argsort(eigs_b))
        if energy_sort:
            wfn = wfn[0,:,col_sort]
    else:
        # Assuming we are working in MO basis, only works for RHF, ROHF trials.
        print("Not correct.")
        sys.exit()
        I = numpy.identity(C.shape[-1], dtype=numpy.float64)
        wfn[0,:,:nalpha] = I[:,:nalpha]
        if uhf:
            print(" # Warning: UHF trial wavefunction can only be used of "
                  "working in ortho AO basis.")
            sys.exit()
    coeff = numpy.array([1.0+0j])
    write_qmcpack_wfn(filename, (coeff,wfn), uhf, (nalpha, nbeta), nmo_tot)
    return nelec


def rediag_fock(fock, X, occ):
    """Rediagonalise Fock matrix.

    Parameters
    ----------
    fock : :class:`numpy.ndarray'
        Fock matrix for given kpoint.
    X : :class:`numpy.ndarray'
        Transformation matrix.
    occ : :class:`numpy.ndarray`
        Boolean array specifying which MOs to occupy.

    Returns
    -------
    eigs : :class:`numpy.array'
        MO eigenvalues.
    eigv : :class:`numpy.ndarray'
        Eigenvectors.
    """
    fock_oao = numpy.dot(X.conj().T, numpy.dot(fock, X))
    e, c = numpy.linalg.eigh(fock_oao)
    return e[occ], c[:,occ]
