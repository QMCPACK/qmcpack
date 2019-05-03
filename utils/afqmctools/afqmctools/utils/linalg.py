import numpy
from pyscf import lib

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
