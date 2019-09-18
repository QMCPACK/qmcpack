import numpy
from afqmctools.utils.linalg import modified_cholesky_direct

def generate_hamiltonian(nmo, nelec, cplx=False, sym=8):
    h1e = numpy.random.random((nmo,nmo))
    if cplx:
        h1e = h1e + 1j*numpy.random.random((nmo,nmo))
    eri = numpy.random.normal(scale=0.01, size=(nmo,nmo,nmo,nmo))
    if cplx:
        eri = eri + 1j*numpy.random.normal(scale=0.01, size=(nmo,nmo,nmo,nmo))
    # Restore symmetry to the integrals.
    if sym >= 4:
        # (ik|jl) = (jl|ik)
        # (ik|jl) = (ki|lj)*
        eri = eri + eri.transpose(2,3,0,1)
        eri = eri + eri.transpose(3,2,1,0).conj()
    if sym == 8:
        eri = eri + eri.transpose(1,0,2,3)
    # Construct hermitian matrix M_{ik,lj}.
    eri = eri.transpose((0,1,3,2))
    eri = eri.reshape((nmo*nmo,nmo*nmo))
    # Make positive semi-definite.
    eri = numpy.dot(eri,eri.conj().T)
    chol = modified_cholesky_direct(eri, tol=1e-4, verbose=False)
    chol = chol.reshape((-1,nmo,nmo))
    enuc = numpy.random.rand()
    return h1e, chol, enuc, eri
