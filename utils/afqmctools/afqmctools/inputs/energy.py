import numpy
from afqmctools.hamiltonian.converter import read_qmcpack_hamiltonian
from afqmctools.wavefunction.converter import read_qmcpack_wavefunction


def calculate_hf_energy(hamil_file, wfn_file):
    hamil = read_qmcpack_hamiltonian(hamil_file)
    wfn, psi0, nelec = read_qmcpack_wavefunction(wfn_file)
    if isinstance(hamil['chol'], numpy.ndarray):
        chol = hamil['chol']
    else:
        chol = hamil['chol'].toarray()
    return local_energy_mol(hamil['hcore'], chol, hamil['enuc'], wfn[1][0], nelec)

def local_energy_mol(hcore, chol, enuc, psi, nelec):
    r"""Calculate local for generic two-body hamiltonian.

    This uses the cholesky decomposed two-electron integrals.

    Parameters
    ----------

    Returns
    -------
    """
    na, nb = nelec
    nbasis = hcore.shape[0]
    assert psi.shape[1] == na + nb
    assert psi.shape[0] == nbasis
    assert chol.shape[0] == nbasis*nbasis

    # (M^2, nchol)
    nchol = chol.shape[-1]
    # Half-rotated green's functions
    O = numpy.dot(psi[:,:na].T, psi[:,:na].conj())
    Ga = numpy.dot(numpy.linalg.inv(O), psi[:,:na].T)
    O = numpy.dot(psi[:,na:].T, psi[:,na:].conj())
    Gb = numpy.dot(numpy.linalg.inv(O), psi[:,na:].T)

    # Construct half rotated cholesky
    # L[ak,n] = \sum_i A*[i,a] L[i,k,n]
    rchol_a = numpy.tensordot(psi[:,:na].conj(),
                              chol.reshape((nbasis,nbasis,-1)),
                              axes=((0),(0))).reshape(nbasis*na,-1)
    rchol_b = numpy.tensordot(psi[:,na:].conj(),
                              chol.reshape((nbasis,nbasis,-1)),
                              axes=((0),(0))).reshape(nbasis*nb,-1)
    Xa = rchol_a.T.dot(Ga.ravel())
    Xb = rchol_b.T.dot(Gb.ravel())
    ecoul = numpy.dot(Xa,Xa)
    ecoul += numpy.dot(Xb,Xb)
    ecoul += 2*numpy.dot(Xa,Xb)
    # T_{abn} = \sum_k Theta_{ak} LL_{ak,n}
    # LL_{ak,n} = \sum_i L_{ik,n} A^*_{ia}
    Ta = numpy.tensordot(Ga, rchol_a.reshape((na,nbasis,-1)), axes=((1),(1)))
    exxa = numpy.tensordot(Ta, Ta, axes=((0,1,2),(1,0,2)))
    Tb = numpy.tensordot(Gb, rchol_b.reshape((nb,nbasis,-1)), axes=((1),(1)))
    exxb = numpy.tensordot(Tb, Tb, axes=((0,1,2),(1,0,2)))
    exx = exxa + exxb
    e2b = 0.5 * (ecoul - exx)
    Ga = numpy.dot(psi[:,:na].conj(), Ga)
    Gb = numpy.dot(psi[:,na:].conj(), Gb)
    e1b = numpy.sum(hcore*(Ga+Gb))
    return (e1b + e2b + enuc, e1b + enuc, e2b)
