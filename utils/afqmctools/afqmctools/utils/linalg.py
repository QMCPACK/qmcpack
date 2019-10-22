import numpy
import time
from pyscf import lib
from afqmctools.utils.io import (
        format_fixed_width_floats,
        format_fixed_width_strings,
        )


def get_ortho_ao(cell, kpts, LINDEP_CUTOFF=0):
    """Generate canonical orthogonalization transformation matrix.

    Parameters
    ----------
    cell : :class:`pyscf.pbc.cell' object.
        PBC cell object.
    kpts : :class:`numpy.array'
        List of kpoints.
    LINDEP_CUTOFF : float
        Linear dependency cutoff. Basis functions whose eigenvalues lie below
        this value are removed from the basis set. Should be set in accordance
        with value in pyscf (pyscf.scf.addons.remove_linear_dep_).

    Returns
    -------
    X : :class:`numpy.array`
        Transformation matrix.
    nmo_per_kpt : :class:`numpy.array`
        Number of OAOs orbitals per kpoint.
    """
    kpts = numpy.reshape(kpts,(-1,3))
    nkpts = len(kpts)
    nao = cell.nao_nr()
    s1e = lib.asarray(cell.pbc_intor('cint1e_ovlp_sph',hermi=1,kpts=kpts))
    X = numpy.zeros((nkpts,nao,nao),dtype=numpy.complex128)
    nmo_per_kpt = numpy.zeros(nkpts,dtype=numpy.int32)
    for k in range(nkpts):
        sdiag, Us = numpy.linalg.eigh(s1e[k])
        nmo_per_kpt[k] = sdiag[sdiag>LINDEP_CUTOFF].size
        norm = numpy.sqrt(sdiag[sdiag>LINDEP_CUTOFF])
        X[k,:,0:nmo_per_kpt[k]] = Us[:,sdiag>LINDEP_CUTOFF] / norm
    return X, nmo_per_kpt

def get_ortho_ao_mol(S, LINDEP_CUTOFF=0):
    """Generate canonical orthogonalization transformation matrix.

    Parameters
    ----------
    S : :class:`numpy.ndarray`
        Overlap matrix.
    LINDEP_CUTOFF : float
        Linear dependency cutoff. Basis functions whose eigenvalues lie below
        this value are removed from the basis set. Should be set in accordance
        with value in pyscf (pyscf.scf.addons.remove_linear_dep_).

    Returns
    -------
    X : :class:`numpy.array`
        Transformation matrix.
    """
    sdiag, Us = numpy.linalg.eigh(S)
    X = Us[:,sdiag>LINDEP_CUTOFF] / numpy.sqrt(sdiag[sdiag>LINDEP_CUTOFF])
    return X

def modified_cholesky_direct(M, tol=1e-5, verbose=False, cmax=20):
    """Modified cholesky decomposition of matrix.

    See, e.g. [Motta17]_

    Parameters
    ----------
    M : :class:`numpy.ndarray`
        Positive semi-definite, symmetric matrix.
    tol : float
        Accuracy desired. Optional. Default : False.
    verbose : bool
        If true print out convergence progress. Optional. Default : False.
    cmax : int
        Number of cholesky vectors to store N_gamma = cmax M.

    Returns
    -------
    chol_vecs : :class:`numpy.ndarray`
        Matrix of cholesky vectors.
    """
    # matrix of residuals.
    delta = numpy.copy(M.diagonal())
    nchol_max = int(cmax*M.shape[0]**0.5)
    # index of largest diagonal element of residual matrix.
    nu = numpy.argmax(numpy.abs(delta))
    delta_max = delta[nu]
    if verbose:
        print ("# max number of cholesky vectors = %d"%nchol_max)
        header = ['iteration', 'max_residual', 'time']
        print(format_fixed_width_strings(header))
        init = [delta_max]
        print('{:17d} '.format(0)+format_fixed_width_floats(init))
        # print ("# iteration %d: delta_max = %f"%(0, delta_max.real))
    # Store for current approximation to input matrix.
    Mapprox = numpy.zeros(M.shape[0], dtype=M.dtype)
    chol_vecs = numpy.zeros((nchol_max, M.shape[0]), dtype=M.dtype)
    nchol = 0
    chol_vecs[0] = numpy.copy(M[:,nu])/delta_max**0.5
    while abs(delta_max) > tol:
        # Update cholesky vector
        start = time.time()
        Mapprox += chol_vecs[nchol]*chol_vecs[nchol].conj()
        delta = M.diagonal() - Mapprox
        nu = numpy.argmax(numpy.abs(delta))
        delta_max = numpy.abs(delta[nu])
        nchol += 1
        Munu0 = numpy.dot(chol_vecs[:nchol,nu].conj(), chol_vecs[:nchol,:])
        chol_vecs[nchol] = (M[:,nu] - Munu0) / (delta_max)**0.5
        if verbose:
            step_time = time.time() - start
            out = [delta_max, step_time]
            print('{:17d} '.format(nchol)+format_fixed_width_floats(out))

    return numpy.array(chol_vecs[:nchol])
