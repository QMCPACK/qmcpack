"""Generate AFQMC data from PYSCF (molecular) simulation."""
import h5py
import numpy
from pyscf import fci
import scipy.sparse
import sys
import time
from afqmctools.utils.io import (
        format_fixed_width_floats,
        format_fixed_width_strings,
        to_qmcpack_complex
        )
from afqmctools.utils.pyscf_utils import load_from_pyscf_chk_mol
from afqmctools.hamiltonian.io import (
        write_qmcpack_sparse,
        write_qmcpack_dense
        )


def write_hamil_mol(scf_data, hamil_file, chol_cut,
                    verbose=True, cas=None,
                    ortho_ao=False, nelec=None,
                    real_chol=False, df=False,
                    dense=False):
    """Write QMCPACK hamiltonian from pyscf scf calculation on mol object.
    """
    hcore, chol_vecs, nelec, enuc, X = generate_hamiltonian(scf_data,
                                                            verbose=verbose,
                                                            chol_cut=chol_cut,
                                                            cas=cas,
                                                            ortho_ao=ortho_ao,
                                                            nelec=nelec,
                                                            df=df)
    # Want L_{(ik),n}
    chol_vecs = chol_vecs.T
    nbasis = hcore.shape[-1]
    if dense:
        write_qmcpack_dense(hcore, chol_vecs, nelec,
                            nbasis, enuc, real_chol=real_chol,
                            filename=hamil_file,
                            ortho=X)
    else:
        write_qmcpack_sparse(hcore, chol_vecs, nelec, nbasis, enuc,
                             filename=hamil_file, real_chol=real_chol,
                             verbose=verbose, ortho=X)

def generate_hamiltonian(scf_data, chol_cut=1e-5, verbose=False, cas=None,
                         ortho_ao=False, nelec=None, df=False):
    # Unpack SCF data.
    # 1. core (1-body) Hamiltonian.
    hcore = scf_data['hcore']
    # 2. Rotation matrix to orthogonalised basis.
    if ortho_ao:
        X = scf_data['X']
    else:
        if scf_data['isUHF']:
            print(" # UHF integrals are not allowed. Use ortho AO option (-a/--ao).")
            sys.exit()
        X = scf_data['mo_coeff']
    df_ints = scf_data.get('df_ints', None)
    C = scf_data['mo_coeff']
    # 3. Pyscf mol object.
    mol = scf_data['mol']
    # Step 1. Rotate core Hamiltonian to orthogonal basis.
    if verbose:
        if ortho_ao:
            print(" # Transforming hcore and eri to ortho AO basis.")
        else:
            print(" # Transforming hcore and eri to MO basis.")
    h1e = numpy.dot(X.T, numpy.dot(hcore, X))
    nbasis = h1e.shape[-1]
    if verbose:
        print(" # Number of basis functions: {}.".format(nbasis))
    # Step 2. Genrate Cholesky decomposed ERIs in non-orthogonal AO basis.
    if df_ints is not None and df:
        chol_vecs = df_ints
        if verbose:
            print(" # Using DF integrals from checkpoint file.")
        assert chol_vecs.shape[1] == nbasis*nbasis
    else:
        if verbose:
            print (" # Performing modified Cholesky decomposition on ERI tensor.")
        chol_vecs = chunked_cholesky(mol, max_error=chol_cut, verbose=verbose)
    if verbose:
        print (" # Orthogonalising Cholesky vectors.")
    start = time.time()
    # Step 2.a Orthogonalise Cholesky vectors.
    chol_trans = ao2mo_chol(chol_vecs, X)
    if verbose:
        print(" # Time to orthogonalise: %f"%(time.time() - start))
    enuc = mol.energy_nuc()
    # Step 3. (Optionally) freeze core / virtuals.
    nelec = mol.nelec
    if cas is not None:
        nfzc = (sum(mol.nelec)-cas[0])//2
        ncas = cas[1]
        if ncas == -1:
            ncas = nbasis - nfzc
        nfzv = nbasis - ncas - nfzc
        h1e, chol_trans, enuc = freeze_core(h1e, chol_trans, enuc, nfzc, ncas,
                                           verbose)
        h1e = h1e[0]
        nelec = (mol.nelec[0]-nfzc, mol.nelec[1]-nfzc)
        # TODO: Don't do this. Pass (N,M) to wavefunction.
        mol.nelec = nelec
        orbs = numpy.identity(h1e.shape[-1])
        orbs = orbs[nfzc:nbasis-nfzv,nfzc:nbasis-nfzv]
        scf_data['mo_coeff'] = C[nfzc:nbasis-nfzv,nfzc:nbasis-nfzv]
    return h1e, chol_trans, nelec, enuc, X


def freeze_core(h1e, chol, ecore, nc, ncas, verbose=True):
    # 1. Construct one-body hamiltonian
    nbasis = h1e.shape[-1]
    chol = chol.reshape((-1,nbasis,nbasis))
    psi = numpy.identity(nbasis)[:,:nc]
    Gcore = gab(psi,psi)
    efzc = local_energy_generic_cholesky(h1e, chol, [Gcore,Gcore], ecore)
    efzc = efzc[0]
    (hc_a, hc_b) = core_contribution_cholesky(chol, [Gcore,Gcore])
    h1e = numpy.array([h1e,h1e])
    h1e[0] = h1e[0] + 2*hc_a
    h1e[1] = h1e[1] + 2*hc_b
    h1e = h1e[:,nc:nc+ncas,nc:nc+ncas]
    nchol = chol.shape[0]
    chol = chol[:,nc:nc+ncas,nc:nc+ncas]
    chol = chol.reshape((nchol,-1))
    # 4. Subtract one-body term from writing H2 as sum of squares.
    if verbose:
        print(" # Number of active orbitals: %d"%ncas)
        print(" # Freezing %d core electrons and %d virtuals."
              %(2*nc, nbasis-nc-ncas))
        print(" # Frozen core energy: %13.8e"%efzc)
    return h1e, chol, efzc

def ao2mo_chol(eri, C):
    nao = C.shape[0]
    nmo = C.shape[1]
    nik = nmo*nmo
    nchol = eri.shape[0]
    eri_ = eri.ravel()
    for i in range(nchol):
        cv = eri[i].reshape(nao,nao)
        half = numpy.dot(cv, C)
        # if nao < nmo we overwrite the data
        eri_[i*nik:(i+1)*nik] = numpy.dot(C.T, half).ravel()
    return eri_[:nchol*nik].reshape((nchol,nik))

def chunked_cholesky(mol, max_error=1e-6, verbose=False, cmax=10):
    """Modified cholesky decomposition from pyscf eris.

    See, e.g. [Motta17]_

    Only works for molecular systems.

    Parameters
    ----------
    mol : :class:`pyscf.mol`
        pyscf mol object.
    orthoAO: :class:`numpy.ndarray`
        Orthogonalising matrix for AOs. (e.g., mo_coeff).
    delta : float
        Accuracy desired.
    verbose : bool
        If true print out convergence progress.
    cmax : int
        nchol = cmax * M, where M is the number of basis functions.
        Controls buffer size for cholesky vectors.

    Returns
    -------
    chol_vecs : :class:`numpy.ndarray`
        Matrix of cholesky vectors in AO basis.
    """
    nao = mol.nao_nr()
    diag = numpy.zeros(nao*nao)
    nchol_max = cmax * nao
    # This shape is more convenient for pauxy.
    chol_vecs = numpy.zeros((nchol_max, nao*nao))
    # eri = numpy.zeros((nao,nao,nao,nao))
    ndiag = 0
    dims = [0]
    nao_per_i = 0
    for i in range(0,mol.nbas):
        l = mol.bas_angular(i)
        nc = mol.bas_nctr(i)
        nao_per_i += (2*l+1)*nc
        dims.append(nao_per_i)
    start = time.time()
    for i in range(0,mol.nbas):
        shls = (i,i+1,0,mol.nbas,i,i+1,0,mol.nbas)
        buf = mol.intor('int2e_sph', shls_slice=shls)
        di, dk, dj, dl = buf.shape
        diag[ndiag:ndiag+di*nao] = buf.reshape(di*nao,di*nao).diagonal()
        ndiag += di * nao
    nu = numpy.argmax(diag)
    delta_max = diag[nu]
    if verbose:
        print(" # Generating Cholesky decomposition of ERIs."%nchol_max)
        print(" # max number of cholesky vectors = %d"%nchol_max)
        header = ['iteration', 'max_residual', 'time']
        print(format_fixed_width_strings(header))
        init = [delta_max, time.time()-start]
        print('{:17d} '.format(0)+format_fixed_width_floats(init))
    j = nu // nao
    l = nu % nao
    sj = numpy.searchsorted(dims, j)
    sl = numpy.searchsorted(dims, l)
    if dims[sj] != j and j != 0:
        sj -= 1
    if dims[sl] != l and l != 0:
        sl -= 1
    Mapprox = numpy.zeros(nao*nao)
    # ERI[:,jl]
    eri_col = mol.intor('int2e_sph',
                         shls_slice=(0,mol.nbas,0,mol.nbas,sj,sj+1,sl,sl+1))
    cj, cl = max(j-dims[sj],0), max(l-dims[sl],0)
    chol_vecs[0] = numpy.copy(eri_col[:,:,cj,cl].reshape(nao*nao)) / delta_max**0.5

    nchol = 0
    while abs(delta_max) > max_error:
        # Update cholesky vector
        start = time.time()
        # M'_ii = \sum_x L_i^x L_i^x
        Mapprox += chol_vecs[nchol] * chol_vecs[nchol]
        # D_ii = M_ii - M'_ii
        delta = diag - Mapprox
        nu = numpy.argmax(numpy.abs(delta))
        delta_max = numpy.abs(delta[nu])
        # Compute ERI chunk.
        # shls_slice computes shells of integrals as determined by the angular
        # momentum of the basis function and the number of contraction
        # coefficients. Need to search for AO index within this shell indexing
        # scheme.
        # AO index.
        j = nu // nao
        l = nu % nao
        # Associated shell index.
        sj = numpy.searchsorted(dims, j)
        sl = numpy.searchsorted(dims, l)
        if dims[sj] != j and j != 0:
            sj -= 1
        if dims[sl] != l and l != 0:
            sl -= 1
        # Compute ERI chunk.
        eri_col = mol.intor('int2e_sph',
                            shls_slice=(0,mol.nbas,0,mol.nbas,sj,sj+1,sl,sl+1))
        # Select correct ERI chunk from shell.
        cj, cl = max(j-dims[sj],0), max(l-dims[sl],0)
        Munu0 = eri_col[:,:,cj,cl].reshape(nao*nao)
        # Updated residual = \sum_x L_i^x L_nu^x
        R = numpy.dot(chol_vecs[:nchol+1,nu], chol_vecs[:nchol+1,:])
        chol_vecs[nchol+1] = (Munu0 - R) / (delta_max)**0.5
        nchol += 1
        if verbose:
            step_time = time.time() - start

            out = [delta_max, step_time]
            print('{:17d} '.format(nchol)+format_fixed_width_floats(out))

    return chol_vecs[:nchol]

def write_qmcpack_trial_wfn(wfn, nelec, filename='wfn.dat'):
    UHF = len(wfn.shape) == 3
    namelist = "&FCI\n UHF = %d\n FullMO \n NCI = 1\n TYPE = matrix\n/"%UHF
    # Single determinant for the moment.
    with open(filename, 'w') as f:
        f.write(namelist+'\n')
        f.write('Coefficients: 1.0\n')
        f.write('Determinant: 1\n')
        nao = wfn.shape[-1]
        if UHF:
            nao = wfn[0].shape[-1]
            write_qmcpack_wfn_single(f, wfn[0], nao)
            nao = wfn[1].shape[-1]
            write_qmcpack_wfn_single(f, wfn[1], nao)
        else:
            write_qmcpack_wfn_single(f, wfn, nao)

def write_qmcpack_wfn_single(out, mos, nao):
    for i in range(0, nao):
        for j in range(0, nao):
            val = mos[i,j]
            out.write('(%.10e,%.10e) '%(val.real, val.imag))
        out.write('\n')

def local_energy_generic_cholesky(h1e, chol_vecs, G, ecore):
    r"""Calculate local for generic two-body hamiltonian.

    This uses the cholesky decomposed two-electron integrals.

    Parameters
    ----------
    system : :class:`hubbard`
        System information for the hubbard model.
    G : :class:`numpy.ndarray`
        Walker's "green's function"

    Returns
    -------
    (E, T, V): tuple
        Local, kinetic and potential energies.
    """
    # Element wise multiplication.
    e1b = numpy.sum(h1e*G[0]) + numpy.sum(h1e*G[1])
    cv = chol_vecs
    ecoul_uu = 0
    ecoul_dd = 0
    ecoul_ud = 0
    ecoul_du = 0
    exx_uu = 0
    exx_dd = 0
    # Below to compute exx_uu/dd we do
    # t1 = numpy.einsum('nik,il->nkl', cv, G[0])
    # t2 = numpy.einsum('nlj,jk->nlk', cv.conj(), G[0])
    # exx_uu = numpy.einsum('nkl,nlk->', t1, t2)
    exx_uu = 0
    for c in cv:
        ecoul_uu += numpy.sum(c*G[0]) * numpy.sum(c.conj().T*G[0])
        ecoul_dd += numpy.sum(c*G[1]) * numpy.sum(c.conj().T*G[1])
        ecoul_ud += numpy.sum(c*G[0]) * numpy.sum(c.conj().T*G[1])
        ecoul_du += numpy.sum(c*G[1]) * numpy.sum(c.conj().T*G[0])
        t1 = numpy.dot(c.T, G[0])
        # print(t1.sum())
        t2 = numpy.dot(c.conj(), G[0])
        # print(t2.sum())
        exx_uu += numpy.einsum('ij,ji->',t1,t2)
        # print("sum:", exx_uu)
        t1 = numpy.dot(c.T, G[1])
        t2 = numpy.dot(c.conj(), G[1])
        exx_dd += numpy.einsum('ij,ji->',t1,t2)
    euu = 0.5*(ecoul_uu-exx_uu)
    edd = 0.5*(ecoul_dd-exx_dd)
    eud = 0.5 * ecoul_ud
    edu = 0.5 * ecoul_du
    e2b = euu + edd + eud + edu
    return (e1b+e2b+ecore, e1b+ecore, e2b)

def core_contribution_cholesky(chol_vecs, G):
    cv = chol_vecs
    hca_j = numpy.einsum('l,lij->ij', numpy.sum(cv*G[0], axis=(1,2)), cv)
    ta_k = numpy.einsum('lpr,pq->lrq', cv, G[0])
    hca_k = 0.5*numpy.einsum('lrq,lsq->rs', ta_k, cv)
    hca = hca_j - hca_k
    hcb_j = numpy.einsum('l,lij->ij', numpy.sum(cv*G[1], axis=(1,2)), cv)
    tb_k = numpy.einsum('lpr,pq->lrq', cv, G[1])
    hcb_k = 0.5*numpy.einsum('lrq,lsq->rs', tb_k, cv)
    hcb = hcb_j - hcb_k
    return (hca, hcb)

def gab(A, B):
    r"""One-particle Green's function.

    This actually returns 1-G since it's more useful, i.e.,

    .. math::
        \langle \phi_A|c_i^{\dagger}c_j|\phi_B\rangle =
        [B(A^{\dagger}B)^{-1}A^{\dagger}]_{ji}

    where :math:`A,B` are the matrices representing the Slater determinants
    :math:`|\psi_{A,B}\rangle`.

    For example, usually A would represent (an element of) the trial wavefunction.

    .. warning::
        Assumes A and B are not orthogonal.

    Parameters
    ----------
    A : :class:`numpy.ndarray`
        Matrix representation of the bra used to construct G.
    B : :class:`numpy.ndarray`
        Matrix representation of the ket used to construct G.

    Returns
    -------
    GAB : :class:`numpy.ndarray`
        (One minus) the green's function.
    """
    # Todo: check energy evaluation at later point, i.e., if this needs to be
    # transposed. Shouldn't matter for Hubbard model.
    inv_O = scipy.linalg.inv((A.conj().T).dot(B))
    GAB = B.dot(inv_O.dot(A.conj().T))
    return GAB
