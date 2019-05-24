"""Write mean-field trial wavefunctions to file."""
import numpy
import h5py
import itertools
import time
import sys
import scipy.linalg
from afqmctools.wavefunction.mol import write_qmcpack_wfn, write_nomsd_wfn

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
    mo_energy = scf_data['mo_energy']
    nmo_tot = numpy.sum(nmo_pk)
    kpts = scf_data['kpts']
    nkpts = len(kpts)
    nalpha = cell.nelec[0] * nkpts
    nbeta = cell.nelec[1] * nkpts
    nelec = nalpha + nbeta
    uhf = scf_data['isUHF']
    # For metallic / system with partial occupancies.
    re_occ, trial, ndeg, srt, isrt = reoccupy(mo_occ, mo_energy, verbose, uhf)
    if uhf:
        nalpha = int(numpy.sum(re_occ[0]))
        nbeta = int(numpy.sum(re_occ[1]))
    else:
        nalpha = nbeta = numpy.sum(re_occ)
    if verbose:
        print(" # Generating QMCPACK trial wavefunction.")
        print(" # Shape of supercell wavefunction "
              "({:d},{:d})".format(nmo_tot, nelec if uhf else nalpha))
        print(" # Number of electrons (nalpha, nbeta) = ({}, {})".format(nalpha, nbeta))
        if uhf:
            print(" # UHF-like trial wavefunction.")
        else:
            print(" # RHF-like trial wavefunction.")
    # For RHF only nalpha entries will be filled.
    wfn = numpy.zeros((1,nmo_tot,nalpha+nbeta), dtype=numpy.complex128)
    eigs_a = []
    eigs_a_occ = []
    eigs_b = []
    eigs_b_occ = []
    occ_b = numpy.zeros(1, dtype=numpy.int64)
    ks = []
    bands = []
    row = 0
    col = 0
    col_b = nalpha
    full_mo_a = []
    full_mo_b = []
    for k in range(nkpts):
        if verbose:
            print(" # Generating trial wavefunction for momentum: "
                  "{:d}".format(k))
        if rediag and ortho_ao:
            start = time.time()
            if uhf:
                occ_a = re_occ[0][k] > 0
                ea, orb_a = rediag_fock(fock[0,k], X[k][:,:nmo_pk[k]], occ_a)
                full_mo_a.append(orb_a)
                eigs_a += list(ea)
                eigs_a_occ +=  list(ea[occ_a])
                occ_b = re_occ[1][k] > 0
                eb, orb_b = rediag_fock(fock[1,k], X[k][:,:nmo_pk[k]], occ_b)
                full_mo_b.append(orb_b)
                eigs_b += list(eb)
                eigs_b_occ += list(eb[occ_b])
                ks.append([k]*len(ea))
                bands.append([i for i in range(len(ea))])
            else:
                occ_a = re_occ[k] > 0
                ea, orb_a = rediag_fock(fock[k], X[k][:,:nmo_pk[k]], occ_a)
                full_mo_a.append(orb_a)
                eigs_a.append(ea)
                eigs_a_occ.append(ea[occ_a])
                ks.append([k]*len(ea))
                bands.append([i for i in range(len(ea))])
                if verbose:
                    print(" # Time to rediagonalise fock: "
                          "  {:13.8e} s".format(time.time()-start))
        elif not rediag and not ortho_ao:
            if uhf:
                print("WARNING: UHF wavefunction with ortho_ao = False not "
                      "allowed.")
                sys.exit()
            occ_a = re_occ[k] > 0
            ea = mo_energy[k]
            eigs_a.append(ea)
            eigs_a_occ.append(ea[occ_a])
            orb_a = numpy.eye(len(ea), dtype=numpy.complex128)
            full_mo_a.append(orb_a)
            ks.append([k]*len(ea))
            bands.append([i for i in range(len(ea))])
        elif ortho_ao and not rediag:
            print("WARNING: ortho_ao = True and rediag = False not implemented.")
            sys.exit()
        row_end = row + nmo_pk[k]
        col_end = col + sum(occ_a)
        wfn[0,row:row_end,col:col_end] = orb_a[:,occ_a].copy()
        if uhf:
            col_end = col_b + sum(occ_b)
            wfn[0,row:row_end,col_b:col_end] = orb_b[:,occ_b].copy()
            col_b += sum(occ_b)
        row += nmo_pk[k]
        col += sum(occ_a)
    # print(wfn[0,:26,:4])
    # Sort columns by energy.
    eigs_a = numpy.array(eigs_a).ravel()
    eigs_a_occ = numpy.array(eigs_a_occ).ravel()
    ks = numpy.array(ks).ravel()
    bands = numpy.array(bands).ravel()
    # TODO: Dangerous should subselect from already sorted array.
    col_sort_occ = numpy.argsort(eigs_a_occ)
    ihomo_a = len(eigs_a_occ[col_sort_occ]) - 1
    if uhf:
        eigs_b = numpy.array(eigs_b).ravel()
        eigs_b_occ = numpy.array(eigs_b_occ).ravel()
        col_sort_b = numpy.argsort(eigs_b)
        col_sort_occ_b = numpy.argsort(eigs_b_occ)
        ihomo_b = len(eigs_b_occ[col_sort_occ_b]) - 1
    if verbose > 1:
        print(" # Recomputed MO energies.")
        if uhf:
            header = ("index", "kpoint", "band", "e_mo(alpha)", "e_mo(beta)")
            fmt = " # {:5s}  {:5s}    {:5s}  {:>17s}  {:>17s}"
            print(fmt.format(*header))
            zipped = zip(ks[srt[0]], bands[srt[0]],
                         eigs_a[srt[0]], eigs_b[srt[1]])
            fmt = " {:5d}   {:5d}     {: 13.8e}    {: 13.8e}"
            for i, t in enumerate(zipped):
                tag = ''
                if i == ihomo_a:
                    tag = ' <--- HOMO(alpha) '
                elif i == ihomo_b:
                    tag += ' <--- HOMO(beta) '
                print(" # {:5d}  ".format(i)+fmt.format(*t)+tag)
        else:
            header = ("index", "kpoint", "band", "e_mo")
            fmt = " # {:5s}  {:5s}    {:5s}  {:>17s}"
            print(fmt.format(*header))
            zipped = zip(ks[srt], bands[srt], eigs_a[isrt])
            fmt = " {:5d}   {:5d}     {: 13.8e}"
            for i, t in enumerate(zipped):
                tag = ''
                if i == ihomo_a:
                    tag = ' <--- HOMO(alpha) '
                print(" # {:5d}  ".format(i)+fmt.format(*t)+tag)
    if energy_sort:
        print("Sorting orbitals by energy.")
        col_sort = list(col_sort_occ)
        if uhf:
            print(" # Warning: UHF trial wavefunction can only be used of "
                  "working in ortho AO basis.")
            sys.exit()
    orb_mat_a = scipy.linalg.block_diag(*full_mo_a)
    orb_mat_b = scipy.linalg.block_diag(*full_mo_b)
    if ndeg > 1:
        trial = trial
    else:
        coeff = numpy.array([1.0+0j])
        trial = (coeff, wfn)
    write_qmcpack_wfn(filename, trial, uhf, (nalpha, nbeta), nmo_tot,
                      orbmat=(orb_mat_a, orb_mat_b))
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
    return e, c

def reoccupy(mo_occ, mo_energy, verbose, uhf):
    if uhf:
        if verbose:
            print(" # Determining occupancies for alpha electrons.")
        re_occ_a, ndeg_a, msd_a, srt_a, isrt_a = (
                determine_occupancies(mo_occ[0],
                                      mo_energy[0],
                                      verbose=verbose>=1)
                )
        if verbose:
            print(" # Determining occupancies for beta electrons.")
        re_occ_b, ndeg_b, msd_b, srt_b, isrt_b = (
                determine_occupancies(mo_occ[1],
                                      mo_energy[1],
                                      verbose=verbose>=1)
                )
        ndeg = max(ndeg_a, ndeg_b)
        nalpha = int(numpy.sum(re_occ_a))
        nbeta = int(numpy.sum(re_occ_b))
        re_occ = [re_occ_a, re_occ_b]
        if ndeg_a > 1:
            occs_a = msd_a
            occs_b = [numpy.where(numpy.ravel(re_occ_b) > 0)[0]]*len(occs_a)
        else:
            occs_a = msd_b
            occs_b = [numpy.where(numpy.ravel(re_occ_a) > 0)[0]]*len(occs_b)
        ndet = len(occs_a)
        msd_coeff = numpy.array([1.0/ndet**0.5]*ndet,
                            dtype=numpy.complex128)
        srt = (srt_a,srt_b)
        isrt = (isrt_a,isrt_b)
        trial = (msd_coeff, occs_a, occs_b)
    else:
        re_occ = mo_occ
        ndeg = 1
        srt = numpy.ravel(mo_energy).argsort()
        isrt = srt.argsort()
        trial = None

    return re_occ, trial, ndeg, srt, isrt

def determine_occupancies(mo_occ, mo_energy, low=0.25,
                          high=0.95, verbose=False, refdet=0,
                          offset=0):
    nocc = 0
    nelec = numpy.sum(mo_occ)
    col_sort = numpy.ravel(mo_energy).argsort()
    col_sort_inv = col_sort.argsort()
    nmo_pk = []
    for occ in mo_occ:
        occ = occ > high
        nmo_pk.append(len(occ))
        nocc += sum(occ)
    if abs(nelec-nocc) < 1e-12:
        if verbose:
            print(" # Found closed shell system.")
        ndeg = 0
        msd = None
        return (mo_occ, ndeg, msd, col_sort,col_sort_inv)
    else:
        if verbose:
            print(" # Found orbital occupancies < 0.25.")
            print(" # Constructing multi determinant trial wavefunction from "
                  "degenerate orbitals.")
        nleft = int(round(nelec-nocc))
        mo_order = numpy.array(mo_occ).ravel()[col_sort]
        deg = (mo_order < high) & (mo_order > low)
        ndeg = sum(deg)
        # Supercell indexed.
        deg_orb = numpy.where(deg)[0]
        combs = [c for c in itertools.combinations(deg_orb, int(nleft))]
        if verbose:
            print(" # Distributing {} electrons in {} orbitals.".format(nleft,ndeg))
            print(" # Number of determinants: {}.".format(len(combs)))
        core = list(numpy.where(mo_order > high)[0])
        core = [c for c in core]
        msd = [core + list(d) for d in combs]
        # Remap to primitive cell (kpt,band) indexing.
        reordered = []
        for i, d in enumerate(msd):
            mo_occ_new = numpy.zeros(len(mo_order))
            mo_occ_new[d] = 1.0
            if i == refdet:
                ref = mo_occ_new[col_sort_inv]
            reordered.append(numpy.where(mo_occ_new[col_sort_inv])[0])
        locc = []
        s = 0
        e = nmo_pk[0]
        for ik in range(len(mo_occ)):
            locc.append(ref[s:e])
            s += nmo_pk[ik]
            try:
                e += nmo_pk[ik+1]
            except IndexError:
                e = -1
        return (locc, ndeg, reordered, col_sort, col_sort_inv)
