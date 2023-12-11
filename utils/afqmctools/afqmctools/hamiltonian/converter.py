import h5py
import numpy
import scipy.sparse
import scipy.linalg
import sys
from afqmctools.utils.io import (
        to_qmcpack_complex,
        from_qmcpack_complex
        )
from afqmctools.hamiltonian.io import (
        write_sparse_basic,
        write_sparse_chol_chunk
        )


def get_dset_simple(fh5, name):
    try:
        dset = fh5[name][:]
        return dset
    except KeyError:
        print("Error reading {:} dataset".format(name))
        return None

def read_fcidump(filename, symmetry=8, verbose=True):
    """Read in integrals from file.

    Parameters
    ----------
    filename : string
        File containing integrals in FCIDUMP format.
    symmetry : int
        Permutational symmetry of two electron integrals.
    verbose : bool
        Controls printing verbosity. Optional. Default: False.

    Returns
    -------
    h1e : :class:`numpy.ndarray`
        One-body part of the Hamiltonian.
    h2e : :class:`numpy.ndarray`
        Two-electron integrals.
    ecore : float
        Core contribution to the total energy.
    nelec : tuple
        Number of electrons.
    """
    assert(symmetry==1 or symmetry==4 or symmetry==8)
    if verbose:
        print ("# Reading integrals in plain text FCIDUMP format.")
    with open(filename) as f:
        while True:
            line = f.readline()
            if 'END' in line or '/' in line:
                break
            for i in line.split(','):
                if 'NORB' in i:
                    nbasis = int(i.split('=')[1])
                elif 'NELEC' in i:
                    nelec = int(i.split('=')[1])
                elif 'MS2' in i:
                    ms2 = int(i.split('=')[1])
        if verbose:
            print("# Number of orbitals: {}".format(nbasis))
            print("# Number of electrons: {}".format(nelec))
        h1e = numpy.zeros((nbasis, nbasis), dtype=numpy.complex128)
        h2e = numpy.zeros((nbasis, nbasis, nbasis, nbasis), dtype=numpy.complex128)
        lines = f.readlines()
        for l in lines:
            s = l.split()
            # ascii fcidump uses Chemist's notation for integrals.
            # each line contains v_{ijkl} i k j l
            # Note (ik|jl) = <ij|kl>.
            if len(s) == 6:
                # FCIDUMP from quantum package.
                integral = float(s[0]) + 1j*float(s[1])
                s = s[1:]
            else:
                try:
                    integral = float(s[0])
                except ValueError:
                    ig = ast.literal_eval(s[0].strip())
                    integral = ig[0] + 1j*ig[1]
            i, k, j, l = [int(x) for x in s[1:]]
            if i == j == k == l == 0:
                ecore = integral
            elif j == 0 and l == 0:
                # <i|k> = <k|i>
                h1e[i-1,k-1] = integral
                h1e[k-1,i-1] = integral.conjugate()
            elif i > 0  and j > 0 and k > 0 and l > 0:
                # Assuming 8 fold symmetry in integrals.
                # <ij|kl> = <ji|lk> = <kl|ij> = <lk|ji> =
                # <kj|il> = <li|jk> = <il|kj> = <jk|li>
                # (ik|jl)
                h2e[i-1,k-1,j-1,l-1] = integral
                if symmetry == 1:
                    continue
                # (jl|ik)
                h2e[j-1,l-1,i-1,k-1] = integral
                # (ki|lj)
                h2e[k-1,i-1,l-1,j-1] = integral.conjugate()
                # (lj|ki)
                h2e[l-1,j-1,k-1,i-1] = integral.conjugate()
                if symmetry == 4:
                    continue
                # (ki|jl)
                h2e[k-1,i-1,j-1,l-1] = integral
                # (lj|ik)
                h2e[l-1,j-1,i-1,k-1] = integral
                # (ik|lj)
                h2e[i-1,k-1,l-1,j-1] = integral
                # (jl|ki)
                h2e[j-1,l-1,k-1,i-1] = integral
    if symmetry == 8:
        if numpy.any(numpy.abs(h1e.imag)) > 1e-18:
            print("# Found complex numbers in one-body Hamiltonian but 8-fold"
                  " symmetry specified.")
        if numpy.any(numpy.abs(h2e.imag)) > 1e-18:
            print("# Found complex numbers in two-body Hamiltonian but 8-fold"
                  " symmetry specified.")
    nalpha = (nelec + ms2) // 2
    nbeta = nalpha - ms2
    return h1e, h2e, ecore, (nalpha, nbeta)

def read_hamil_type(filename):
    ham_type = None
    with h5py.File(filename, 'r') as fh5:
        try:
            ds = fh5['Hamiltonian/DenseFactorized/L']
            return 'dense'
        except KeyError:
            ham_type = None
        try:
            ds = fh5['Hamiltonian/Factorized/vals_0']
            return 'sparse'
        except KeyError:
            ham_type = None
        try:
            ds = fh5['Hamiltonian/KPFactorized/L0']
            return 'kpoint'
        except KeyError:
            ham_type = None
        try:
            ds = fh5['Hamiltonian/THC/Luv']
            return 'thc'
        except KeyError:
            ham_type = None
    return ham_type

def read_qmcpack_hamiltonian(filename, get_chol=True):
    """Read Hamiltonian from QMCPACK format.

    Parameters
    ----------
    filename : string
        QMPACK Hamiltonian file.

    Returns
    -------
    hamil : dict
        Data read from file.
    """
    ham_type = read_hamil_type(filename)
    if ham_type == 'kpoint':
        hc, chol, enuc, nmo, nelec, nmok, qkk2, nchol_pk, minus_k = (
                read_qmcpack_cholesky_kpoint(filename, get_chol=get_chol)
                )
        hamil = {
            'hcore': hc,
            'chol': chol,
            'enuc': enuc,
            'nelec': nelec,
            'nmo': nmo,
            'nmo_pk': nmok,
            'nchol_pk': nchol_pk,
            'minus_k': minus_k,
            'qk_k2': qkk2
            }
    elif ham_type == 'sparse':
        hc, chol, enuc, nmo, nelec = read_qmcpack_sparse(filename,
                                                         get_chol=get_chol)
        hamil = {
            'hcore': hc,
            'chol': chol,
            'enuc': enuc,
            'nmo': nmo,
            'nelec': nelec
            }
        return hamil
    elif ham_type == 'dense':
        hcore, chol, enuc = read_qmcpack_dense(filename)
        hamil = {'hcore': hcore, 'chol': chol, 'enuc': enuc}
        return hamil
    elif ham_type == 'thc':
        print(" Please implement THC hamiltonian reader.")
        hamil = {}
    else:
        hamil = None
        raise TypeError("Unknown Hamiltonian type in read_qmcpack_hamiltonian.")
    return hamil

def read_qmcpack_sparse(filename, get_chol=True):
    """Read in integrals from hdf5.

    Parameters
    ----------
    filename : string
        File containing integrals in qmcpack format.

    Returns
    -------
    hcore : :class:`numpy.ndarray`
        One-body part of the Hamiltonian.
    chol_vecs : :class:`scipy.sparse.csr_matrix`
        Two-electron integrals. Shape: [nmo*nmo, nchol]
    ecore : float
        Core contribution to the total energy.
    nmo : int
        Number of orbitals.
    nelec : tuple
        Number of electrons.
    """
    enuc, dims, hcore, real_ints = read_common_input(filename)
    chunks = dims[2]
    nmo = dims[3]
    nalpha = dims[4]
    nbeta = dims[5]
    nchol = dims[7]
    if get_chol:
        chol_vecs = read_cholesky(filename, real_ints=real_ints)
    else:
        chol_vecs = None
    return (hcore, chol_vecs, enuc, int(nmo), (int(nalpha), int(nbeta)))

def read_cholesky(filename, full=True, ichunk=None, real_ints=False):
    with h5py.File(filename, 'r') as fh5:
        try:
            block_sizes = fh5['Hamiltonian/Factorized/block_sizes'][:]
        except KeyError:
            print(" Error reading Hamiltonian/Factorized/block_sizes")
        nval = sum(block_sizes)
        s = 0
        dims = fh5['Hamiltonian/dims'][:]
        nmo = dims[3]
        nchol = dims[-1]
        if full:
            if real_ints:
                vals = numpy.zeros(nval, dtype=numpy.float64)
            else:
                vals = numpy.zeros(nval, dtype=numpy.complex128)
            row_ix = numpy.zeros(nval, dtype=numpy.int32)
            col_ix = numpy.zeros(nval, dtype=numpy.int32)
            for ic, bs in enumerate(block_sizes):
                ixs, vchunk = get_chunk(fh5, ic, real_ints)
                row_ix[s:s+bs] = ixs[::2]
                col_ix[s:s+bs] = ixs[1::2]
                vals[s:s+bs] = vchunk
                s += bs
            chol_vecs = scipy.sparse.csr_matrix((vals, (row_ix, col_ix)),
                                                shape=(nmo*nmo,nchol))
            return chol_vecs
        else:
            return get_chunk(fh5, ichunk, real_ints)

def get_chunk(fh5, ichunk, real_ints):
    ixs = fh5['Hamiltonian/Factorized/index_%i'%ichunk][:]
    if real_ints:
        try:
            vals = fh5['Hamiltonian/Factorized/vals_%i'%ichunk][:].ravel()
        except KeyError:
            print(" Error reading Hamiltonian/Factorized/vals_{:}".format(ichunk))
    else:
        try:
            vals = fh5['Hamiltonian/Factorized/vals_%i'%ichunk][:].view(numpy.complex128).ravel()
        except KeyError:
            print(" Error reading Hamiltonian/Factorized/vals_{:}".format(ichunk))
    return ixs, vals

def check_sym(ikjl, nmo, sym):
    """Check permutational symmetry of integral

    Parameters
    ----------
    ikjl : tuple of ints
        Orbital indices of ERI.
    nmo : int
        Number of orbitals
    sym : int
        Desired permutational symmetry to check.

    Returns
    -------
    sym_allowed : bool
        True if integral is unique from set of equivalent.
    """
    if sym == 1:
        return True
    else:
        i, k, j, l = ikjl
        if sym == 4:
            kilj = (k,i,l,j)
            jlik = (j,l,i,k)
            ljki = (l,j,k,i)
            if (ikjl > jlik) or (ikjl > kilj) or (ikjl > ljki):
                return False
            else:
                return True
        else:
            ik = i + k*nmo
            jl = j + l*nmo
            return (i >= k and j >= l) and ik >= jl

def fmt_integral(intg, i, k, j, l, cplx, paren=False):
    if cplx:
        if paren:
            fmt = '  ({: 13.8e}, {: 13.8e}) {:4d}  {:4d}  {:4d}  {:4d}\n'
        else:
            fmt = '  {: 13.8e}    {: 13.8e}  {:4d}  {:4d}  {:4d}  {:4d}\n'
        out = fmt.format(intg.real, intg.imag, i+1, k+1, j+1, l+1)
    else:
        fmt = '  {: 13.8e}    {:4d}  {:4d}  {:4d}  {:4d}\n'
        out = fmt.format(intg.real, i+1, k+1, j+1, l+1)
    return out


def write_fcidump(filename, hcore, chol, enuc, nmo, nelec, tol=1e-8,
                  sym=1, cplx=True, paren=False):
    """Write FCIDUMP based from Cholesky factorised integrals.

    Parameters
    ----------
    filename : string
        Filename to write FCIDUMP to.
    hcore : :class:`numpy.ndarray`
        One-body hamiltonian.
    chol : :class:`numpy.ndarray`
        Cholesky matrix L[ik,n]
    enuc : float
        Nuclear repulsion energy.
    nmo : int
        Total number of MOs.
    nelec : tuple
        Number of alpha and beta electrons.
    tol : float
        Only print eris above tol. Optional. Default 1e-8.
    sym : int
        Controls whether to only print symmetry inequivalent ERIS.
        Optional. Default 1, i.e. print everything.
    cplx : bool
        Write in complex format. Optional. Default : True.
    paren : bool
        Write complex numbers in parenthesis.
    """
    header = fcidump_header(sum(nelec), nmo, nelec[0]-nelec[1])
    if cplx and sym > 4:
        print("Warning: Requested 8-fold permutational "
              "symmetry with complex integrals.")
        cplx = False

    with open(filename, 'w') as f:
        f.write(header)
        # Generate M_{(ik),(lj)} = (ik|jl)
        eris = chol.dot(chol.conj().T).toarray().reshape((nmo,nmo,nmo,nmo))
        for i in range(0,nmo):
            for k in range(0,nmo):
                for j in range(0,nmo):
                    for l in range(0,nmo):
                        sym_allowed = check_sym((i,k,j,l), nmo, sym)
                        if abs(eris[i,k,l,j]) > tol and sym_allowed:
                            if not cplx:
                                if abs(eris[i,k,l,j].imag > 1e-12):
                                    print("# Found complex integrals with cplx==False.")
                                    #sys.exit()
                            out = fmt_integral(eris[i,k,l,j], i, k, j, l,
                                               cplx, paren=paren)
                            f.write(out)
        for i in range(0,nmo):
            for j in range(0,i+1):
                if abs(hcore[i,j]) > tol:
                    out = fmt_integral(hcore[i,j], i, j, -1, -1,
                                       cplx, paren=paren)
                    f.write(out)

        f.write(fmt_integral(enuc+0j,-1,-1,-1,-1, cplx, paren=paren))


def read_qmcpack_cholesky_kpoint(filename, get_chol=True):
    """Read in integrals from qmcpack hdf5 format. kpoint dependent case.

    Parameters
    ----------
    filename : string
        File containing integrals in qmcpack format.

    Returns
    -------
    hcore : :class:`numpy.ndarray`
        One-body part of the Hamiltonian.
    chol_vecs : :class:`scipy.sparse.csr_matrix`
        Two-electron integrals. Shape: [nmo*nmo, nchol]
    ecore : float
        Core contribution to the total energy.
    nmo : int
        Number of orbitals.
    nelec : tuple
        Number of electrons.
    nmo_pk : :class:`numpy.ndarray`
        Number of orbitals per kpoint.
    qk_k2 : :class:`numpy.ndarray`
        Array mapping (q,k) pair to kpoint: Q = k_i - k_k + G.
        qk_k2[iQ,ik_i] = i_kk.
    """
    enuc, dims, hcore, real_ints = read_common_input(filename, get_hcore=False)
    with h5py.File(filename, 'r') as fh5:
        nmo_pk = get_dset_simple(fh5, 'Hamiltonian/NMOPerKP')
        nchol_pk = get_dset_simple(fh5, 'Hamiltonian/NCholPerKP')
        qk_k2 = get_dset_simple(fh5, 'Hamiltonian/QKTok2')
        minus_k = get_dset_simple(fh5, 'Hamiltonian/MinusK')
        hcore = []
        nkp = dims[2]
        nmo_tot = dims[3]
        nalpha = dims[4]
        nbeta = dims[5]
        if nmo_pk is None:
            raise KeyError("Could not read NMOPerKP dataset.")
        for i in range(0, nkp):
            hk = get_dset_simple(fh5, 'Hamiltonian/H1_kp{}'.format(i))
            if hk is None:
                raise KeyError("Could not read one-body hamiltonian.")
            nmo = nmo_pk[i]
            hcore.append(hk.view(numpy.complex128).reshape(nmo,nmo))
        chol_vecs = []
        if nmo_pk is None:
            raise KeyError("Error nmo_pk dataset does not exist.")
        nmo_max = max(nmo_pk)
    if minus_k is None:
        raise KeyError("Error MinusK dataset does not exist.")
    if nchol_pk is None:
        raise KeyError("Error NCholPerKP dataset does not exist.")
    if get_chol:
        for i in range(0, nkp):
            chol_vecs.append(get_kpoint_chol(filename, nchol_pk, minus_k, i))
    else:
        chol_vecs = None

    return (hcore, chol_vecs, enuc, int(nmo_tot), (int(nalpha), int(nbeta)),
            nmo_pk, qk_k2, nchol_pk, minus_k)

def get_kpoint_chol(filename, nchol_pk, minus_k, i):
    with h5py.File(filename, 'r') as fh5:
        try:
            Lk = get_dset_simple(fh5, 'Hamiltonian/KPFactorized/L{}'.format(i))
            if Lk is None:
                raise TypeError("Could not read Cholesky integrals from file.")
            nchol = nchol_pk[i]
            Lk = Lk.view(numpy.complex128)[:,:,0]
        except KeyError:
            Lk = get_dset_simple(fh5, 'Hamiltonian/KPFactorized/L{}'.format(i))
            if Lk is None:
                raise TypeError("Could not read Cholesky integrals from file.")
            nchol = nchol_pk[minus_k[i]]
            Lk = Lk.view(numpy.complex128).conj()[:,:,0]
    return Lk

def read_qmcpack_dense(filename):
    """Read in integrals from qmcpack hdf5 format.

    Parameters
    ----------
    filename : string
        File containing integrals in qmcpack format.

    Returns
    -------
    hcore : :class:`numpy.ndarray`
        One-body part of the Hamiltonian.
    chol_vecs : :class:`numpy.ndarray`
        Two-electron integrals. Shape: [nmo*nmo, nchol]
    ecore : float
        Core contribution to the total energy.
    """
    real_ints = False
    enuc, dims, hcore, real_ints = read_common_input(filename)
    with h5py.File(filename, 'r') as fh5:
        nmo = dims[3]
        nchol = dims[-1]
        if real_ints:
            try:
                chol = fh5['Hamiltonian/DenseFactorized/L'][:]
            except KeyError:
                print("Error Hamiltonian/DenseFactorized/L")
        else:
            try:
                chol = from_qmcpack_complex(fh5['Hamiltonian/DenseFactorized/L'][:],
                                            (nmo*nmo,nchol))
            except KeyError:
                print("Error Hamiltonian/DenseFactorized/L")
    assert chol.shape == (nmo*nmo, nchol)
    assert hcore.shape == (nmo, nmo)
    return hcore, chol, enuc

def read_common_input(filename, get_hcore=True):
    with h5py.File(filename, 'r') as fh5:
        try:
            enuc = fh5['Hamiltonian/Energies'][:][0]
        except:
            print(" Error reading Hamiltonian/Energies dataset.")
            enuc = None
        try:
            dims = fh5['Hamiltonian/dims'][:]
        except:
            dims = None
        assert dims is not None, "Error reading Hamiltonian/dims data set."
        assert len(dims) == 8, "Hamiltonian dims data set has incorrect length."
        nmo = dims[3]
        real_ints = False
        if get_hcore:
            try:
                hcore = from_qmcpack_complex(fh5['Hamiltonian/hcore'][:], (nmo,nmo))
            except ValueError:
                hcore = fh5['Hamiltonian/hcore'][:]
                real_ints = True
            except KeyError:
                hcore = None
                pass
            if hcore is None:
                try:
                    # old sparse format only for complex.
                    hcore = fh5['Hamiltonian/H1'][:].view(numpy.complex128).ravel()
                    idx = fh5['Hamiltonian/H1_indx'][:]
                    row_ix = idx[::2]
                    col_ix = idx[1::2]
                    hcore = scipy.sparse.csr_matrix((hcore, (row_ix, col_ix))).toarray()
                    hcore = numpy.tril(hcore, -1) + numpy.tril(hcore, 0).conj().T
                except:
                    hcore = None
                    print("Error reading Hamiltonian/hcore data set.")
            try:
                complex_ints = bool(fh5['Hamiltonian/ComplexIntegrals'][0])
            except KeyError:
                complex_ints = None

            if complex_ints is not None:
                if hcore is not None:
                    hc_type = hcore.dtype
                else:
                    hc_type = type(hcore)
                msg = ("ComplexIntegrals flag conflicts with integral data type. "
                       "dtype = {:} ComplexIntegrals = {:}.".format(hc_type, complex_ints))
                assert real_ints ^ complex_ints, msg
        else:
            hcore = None
    return enuc, dims, hcore, real_ints

def fcidump_header(nel, norb, spin):
    header = (
        "&FCI " +
        "NORB={:d}, ".format(norb) +
        "NELEC={:d}, ".format(nel) +
        "MS2={:d},\n".format(spin) +
        "ORBSYM=" +
        ",".join([str(1)]*norb) +
        ",\n" +
        "ISYM=1\n" +
        "&END\n"
        )
    return header


def write_fcidump_kpoint(filename, hcore, chol, enuc, nmo_tot, nelec,
                         nmo_pk, qk_k2, tol=1e-8, sym=1, paren=False,
                         cplx=True):
    """Write FCIDUMP based from Cholesky factorised integrals.

    Parameters
    ----------
    filename : string
        Filename to write FCIDUMP to.
    hcore : list
        One-body hamiltonian.
    chol : list
        Cholesky matrices L[Q][k_i][i,k]
    enuc : float
        Nuclear repulsion energy.
    nmo_tot : int
        Total number of MOs.
    nelec : tuple
        Number of alpha and beta electrons.
    nmo_pk : :class:`numpy.ndarray`
        Number of MOs per kpoint.
    qk_k2 : :class:`numpy.ndarray`
        Array mapping (q,k) pair to kpoint: Q = k_i - k_k + G.
        qk_k2[iQ,ik_i] = i_kk.
    tol : float
        Only print eris above tol. Optional. Default 1e-8.
    sym : int
        Controls whether to only print symmetry inequivalent ERIS.
        Optional. Default 1, i.e. print everything.
    paren : bool
        Write complex numbers in parenthesis.
    """
    header = fcidump_header(sum(nelec), nmo_tot, nelec[0]-nelec[1])
    nkp = len(nmo_pk)
    offsets = numpy.zeros(nkp, dtype=numpy.int32)
    for i in range(1,nkp):
        offsets[i] = offsets[i-1] + nmo_pk[i-1]

    with open(filename, 'w') as f:
        f.write(header)
        for iq, lq in enumerate(chol):
            for ki in range(nkp):
                for kl in range(nkp):
                    eri = numpy.dot(lq[ki], lq[kl].conj().T)
                    if not cplx:
                        if len(numpy.where(numpy.abs(eri.imag) > 1e-12)) > 0:
                            print("# Found complex integrals with cplx==False.")
                            #sys.exit()
                    ik = 0
                    for i in range(0, nmo_pk[ki]):
                        kk = qk_k2[iq,ki]
                        I = i + offsets[ki]
                        for k in range(0, nmo_pk[kk]):
                            kj = qk_k2[iq,kl]
                            K = k + offsets[kk]
                            lj = 0
                            for l in range(0, nmo_pk[kl]):
                                L = l + offsets[kl]
                                for j in range(0, nmo_pk[kj]):
                                    J = j + offsets[kj]
                                    sym_allowed = check_sym((I,K,J,L),
                                                            nmo_tot, sym)
                                    if abs(eri[ik,lj]) > tol and sym_allowed:
                                        out = fmt_integral(eri[ik,lj],
                                                           I, K, J, L,
                                                           cplx, paren=paren)
                                        f.write(out)
                                    lj += 1
                            ik += 1

        for ik, hk in enumerate(hcore):
            for i in range(nmo_pk[ik]):
                I = i + offsets[ik]
                for j in range(nmo_pk[ik]):
                    J = j + offsets[ik]
                    if I >= J and abs(hk[i,j]) > tol:
                        out = fmt_integral(hk[i,j], I, J, -1, -1,
                                           cplx, paren=paren)
                        f.write(out)

        out = fmt_integral(enuc+0j, -1, -1, -1, -1, cplx, paren=paren)
        f.write(out)

def sparse_to_dense(sparse_file, dense_file, real_chol=False):
    hcore, chol, enuc, nmo, nelec = read_qmcpack_sparse(sparse_file,
                                                        get_chol=False)
    with h5py.File(dense_file, 'w') as fh5:
        fh5['Hamiltonian/Energies'] = numpy.array([enuc,0])
        with h5py.File(sparse_file, 'r') as sph5:
            block_sizes = sph5['Hamiltonian/Factorized/block_sizes'][:]
            dims = sph5['Hamiltonian/dims'][:]
            nchol = dims[-1]
            if real_chol:
                fh5['Hamiltonian/hcore'] = hcore.real
                shape = (nmo*nmo,nchol)
            else:
                fh5['Hamiltonian/hcore'] = to_qmcpack_complex(hcore.astype(numpy.complex128))
                shape = (nmo*nmo,nchol,2)
            fh5['Hamiltonian/dims'] = dims
            chol_dset = fh5.create_dataset('Hamiltonian/DenseFactorized/L',
                                           shape,
                                           dtype=numpy.float64)
            s = 0
            for ic, bs in enumerate(block_sizes):
                ixs, vchunk = get_chunk(sph5, ic, real_chol)
                row_ix, col_ix = ixs[::2], ixs[1::2]
                sort = numpy.argsort(row_ix)
                row_ix = row_ix[sort]
                col_ix = col_ix[sort]
                vchunk = vchunk[sort]
                rows = numpy.unique(row_ix)
                for r in rows:
                    # H5PY array slicing is restrictive.
                    # index has to be in increasing order
                    rix = row_ix == r
                    cr = col_ix[rix]
                    vr = vchunk[rix]
                    sort = numpy.argsort(cr)
                    if real_chol:
                        chol_dset[r,cr[sort]] = numpy.real(vr[sort])
                    else:
                        chol_dset[r,cr[sort]] = to_qmcpack_complex(vr[sort])

def kpoint_to_sparse(kp_file, sp_file, real_chol=False,
                     verbose=False, thresh=1e-8):
    hamil = read_qmcpack_hamiltonian(kp_file, get_chol=False)
    hcore = hamil['hcore']
    enuc = hamil['enuc']
    nelec = hamil['nelec']
    nmo = hamil['nmo']
    nmo_pk = hamil['nmo_pk']
    minus_k = hamil['minus_k']
    qk_k2 = hamil['qk_k2']
    nkp = len(hcore)
    # 1. Extend hcore to block diagonal.
    hcore = scipy.linalg.block_diag(*hcore)
    # 2. Unpack chol
    dtype = numpy.float64 if real_chol else numpy.complex128
    nchol_pk = hamil['nchol_pk']
    for i, nc in enumerate(nchol_pk):
        if nc == 0:
            nchol_pk[i] = nchol_pk[minus_k[i]]
    nchol = sum(nchol_pk)
    orb_offset = numpy.zeros(nkp, dtype=numpy.int32)
    chol_offset = numpy.zeros(nkp, dtype=numpy.int32)
    write_sparse_basic(sp_file, hcore, enuc, nelec, real_chol=real_chol)
    for i in range(1, nkp):
        orb_offset[i] = orb_offset[i-1] + nmo_pk[i-1]
        chol_offset[i] = chol_offset[i-1] + nchol_pk[i-1]
    block_sizes = []
    for iq in range(nkp):
        s = chol_offset[iq]
        e = s + nchol_pk[iq]
        lq = get_kpoint_chol(kp_file, nchol_pk, minus_k, iq)
        vals = []
        ixs = []
        for ki in range(nkp):
            for i in range(0, nmo_pk[ki]):
                kk = qk_k2[iq,ki]
                I = i + orb_offset[ki]
                for k in range(0, nmo_pk[kk]):
                    K = k + orb_offset[kk]
                    for nc in range(0,nchol_pk[iq]):
                        ll = lq[ki,(i*nmo_pk[ki]+k)*nchol_pk[iq]+nc]
                        if abs(ll) > thresh:
                            vals.append(ll)
                            ixs.append([I*nmo+K,s+nc])
        block_sizes.append(len(vals))
        write_sparse_chol_chunk(numpy.array(ixs).ravel(), vals, iq, filename=sp_file)

    with h5py.File(sp_file, 'a') as fh5:
        fh5['Hamiltonian/Factorized/block_sizes'] = block_sizes
        fh5['Hamiltonian/dims'] = numpy.array([0, sum(block_sizes),
                                               nkp, nmo, nelec[0], nelec[1], 0,
                                               nchol])
