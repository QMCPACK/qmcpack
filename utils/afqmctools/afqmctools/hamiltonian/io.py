import h5py
import numpy
from afqmctools.utils.io import to_qmcpack_complex

def write_qmcpack_sparse(hcore, chol, nelec, nmo, e0=0.0, filename='hamiltonian.h5',
                 real_chol=False, verbose=False, cutoff=1e-16, ortho=None):
    with h5py.File(filename, 'w') as fh5:
        fh5['Hamiltonian/Energies'] = numpy.array([e0,0])
        if real_chol:
            fh5['Hamiltonian/hcore'] = hcore
        else:
            shape = hcore.shape
            hcore = hcore.astype(numpy.complex128).view(numpy.float64)
            hcore = hcore.reshape(shape+(2,))
            fh5['Hamiltonian/hcore'] = hcore
        if ortho is not None:
            fh5['Hamiltonian/X'] = ortho
        # number of cholesky vectors
        nchol_vecs = chol.shape[-1]
        ix, vals = to_sparse(chol, cutoff=cutoff)
        nnz = len(vals)
        mem = (8 if real_chol else 16) * nnz / (1024.0**3)
        if verbose:
            print(" # Total number of non-zero elements in sparse cholesky ERI"
                   " tensor: %d"%nnz)
            nelem = chol.shape[0]*chol.shape[1]
            print(" # Sparsity of ERI Cholesky tensor: "
                   "%f"%(1-float(nnz)/nelem))
            print(" # Total memory required for ERI tensor: %13.8e GB"%(mem))
        fh5['Hamiltonian/Factorized/block_sizes'] = numpy.array([nnz])
        fh5['Hamiltonian/Factorized/index_0'] = numpy.array(ix)
        if real_chol:
            fh5['Hamiltonian/Factorized/vals_0'] = numpy.array(vals)
        else:
            fh5['Hamiltonian/Factorized/vals_0'] = (
                    to_qmcpack_complex(numpy.array(vals, dtype=numpy.complex128))
                    )
        # Number of integral blocks used for chunked HDF5 storage.
        # Currently hardcoded for simplicity.
        nint_block = 1
        (nalpha, nbeta) = nelec
        unused = 0
        fh5['Hamiltonian/dims'] = numpy.array([unused, nnz, nint_block, nmo,
                                               nalpha, nbeta, unused,
                                               nchol_vecs])
        fh5['Hamiltonian/ComplexIntegrals'] = numpy.array([not int(real_chol)],
                                                          dtype=numpy.int32)
        # TODO: FDM This is deprecated, remove eventually.
        occups = [i for i in range(0, nalpha)]
        occups += [i+nmo for i in range(0, nbeta)]
        fh5['Hamiltonian/occups'] = numpy.array(occups)

def write_sparse_chol_chunk(ixs, vals, chunk, filename='hamiltonian.h5', real_chol=False):
    with h5py.File(filename, 'a') as fh5:
        fh5['Hamiltonian/Factorized/index_{:d}'.format(chunk)] = numpy.array(ixs)
        if real_chol:
            fh5['Hamiltonian/Factorized/vals_{:d}'.format(chunk)] = numpy.array(chunk)
        else:
            fh5['Hamiltonian/Factorized/vals_{:d}'.format(chunk)] = (
                    to_qmcpack_complex(numpy.array(vals, dtype=numpy.complex128))
                    )

def write_sparse_basic(filename, hcore, e0, nelec, real_chol=False):
    nalpha, nbeta = nelec
    nmo = hcore.shape[-1]
    with h5py.File(filename, 'w') as fh5:
        fh5['Hamiltonian/Energies'] = numpy.array([e0,0])
        if real_chol:
            fh5['Hamiltonian/hcore'] = hcore
        else:
            shape = hcore.shape
            hcore = hcore.astype(numpy.complex128).view(numpy.float64)
            hcore = hcore.reshape(shape+(2,))
            fh5['Hamiltonian/hcore'] = hcore
        occups = [i for i in range(0, nalpha)]
        occups += [i+nmo for i in range(0, nbeta)]
        fh5['Hamiltonian/occups'] = numpy.array(occups)

def write_qmcpack_dense(hcore, chol, nelec, nmo, enuc=0.0,
                        filename='hamiltonian.h5', real_chol=False,
                        verbose=False, ortho=None):
    with h5py.File(filename, 'w') as fh5:
        fh5['Hamiltonian/Energies'] = numpy.array([enuc,0])
        if real_chol:
            fh5['Hamiltonian/hcore'] = numpy.real(hcore)
            fh5['Hamiltonian/DenseFactorized/L'] = numpy.real(chol)
        else:
            fh5['Hamiltonian/hcore'] = to_qmcpack_complex(hcore.astype(numpy.complex128))
            fh5['Hamiltonian/DenseFactorized/L'] = to_qmcpack_complex(chol.astype(numpy.complex128))
        fh5['Hamiltonian/dims'] = numpy.array([0, 0, 0, nmo,
                                               nelec[0], nelec[1], 0,
                                               chol.shape[-1]])
        fh5['Hamiltonian/ComplexIntegrals'] = numpy.array([not int(real_chol)],
                                                          dtype=numpy.int32)
        if ortho is not None:
            fh5['Hamiltonian/X'] = ortho


def to_sparse(vals, cutoff=1e-8):
    nz = numpy.where(numpy.abs(vals) > cutoff)
    ix = numpy.empty(nz[0].size+nz[1].size, dtype=numpy.int32)
    ix[0::2] = nz[0]
    ix[1::2] = nz[1]
    vals = vals[nz]
    return ix, vals
