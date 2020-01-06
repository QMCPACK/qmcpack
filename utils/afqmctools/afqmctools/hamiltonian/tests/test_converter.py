import numpy
import os
import scipy.sparse
import unittest
from afqmctools.hamiltonian.converter import (
        read_qmcpack_hamiltonian,
        read_fcidump,
        write_fcidump
        )
from afqmctools.utils.linalg import modified_cholesky_direct
from afqmctools.hamiltonian.io import write_qmcpack_sparse
from afqmctools.utils.testing import generate_hamiltonian

numpy.random.seed(7)

class TestConverter(unittest.TestCase):

    def test_convert_real(self):
        nmo = 17
        nelec = (3,3)
        h1e, chol, enuc, eri = generate_hamiltonian(nmo, nelec, cplx=False, sym=8)
        write_qmcpack_sparse(h1e, chol.reshape((-1,nmo*nmo)).T.copy(),
                             nelec, nmo, e0=enuc, real_chol=True)
        hamil = read_qmcpack_hamiltonian('hamiltonian.h5')
        write_fcidump('FCIDUMP', hamil['hcore'], hamil['chol'], hamil['enuc'],
                      hamil['nmo'], hamil['nelec'], sym=8, cplx=False)
        h1e_r, eri_r, enuc_r, nelec_r = read_fcidump('FCIDUMP', verbose=False)
        dm = numpy.zeros((nmo,nmo))
        dm[(0,1,2),(0,1,2)] = 1.0
        eri_r = eri_r.transpose((0,1,3,2)).reshape((nmo*nmo,nmo*nmo))
        chol_r = modified_cholesky_direct(eri_r, tol=1e-8, verbose=False)
        chol_r = chol_r.reshape((-1,nmo,nmo))
        self.assertAlmostEqual(numpy.einsum('ij,ij->', dm, h1e-h1e_r).real, 0.0)
        self.assertAlmostEqual(numpy.einsum('ij,nij->', dm, chol-chol_r).real, 0.0)
        # Test integral only appears once in file.
        h1e_r, eri_r, enuc_r, nelec_r = read_fcidump('FCIDUMP', symmetry=1,
                                                     verbose=False)
        i,j,k,l = (0,1,2,3)
        combs = [
            (i,j,k,l),
            (k,l,i,j),
            (j,i,l,k),
            (l,k,j,i),
            (j,i,k,l),
            (l,k,i,j),
            (i,j,l,k),
            (k,l,i,j),
            ]
        for c in combs:
            if abs(eri_r[c]) > 0:
                self.assertEqual(c,(l,k,j,i))

    def test_convert_cplx(self):
        nmo = 17
        nelec = (3,3)
        h1e, chol, enuc, eri = generate_hamiltonian(nmo, nelec, cplx=True, sym=4)
        write_qmcpack_sparse(h1e, chol.reshape((-1,nmo*nmo)).T.copy(),
                             nelec, nmo, e0=enuc, real_chol=False)
        hamil = read_qmcpack_hamiltonian('hamiltonian.h5')
        write_fcidump('FCIDUMP', hamil['hcore'], hamil['chol'], hamil['enuc'],
                      hamil['nmo'], hamil['nelec'], sym=4, cplx=True)
        h1e_r, eri_r, enuc_r, nelec_r = read_fcidump('FCIDUMP', symmetry=4,
                                                     verbose=False)
        dm = numpy.zeros((nmo,nmo))
        dm[(0,1,2),(0,1,2)] = 1.0
        eri_r = eri_r.transpose((0,1,3,2)).reshape((nmo*nmo,nmo*nmo))
        chol_r = modified_cholesky_direct(eri_r, tol=1e-8, verbose=False)
        chol_r = chol_r.reshape((-1,nmo,nmo))
        self.assertAlmostEqual(numpy.einsum('ij,ij->', dm, h1e-h1e_r).real, 0.0)
        self.assertAlmostEqual(numpy.einsum('ij,nij->', dm, chol-chol_r).real, 0.0)
        # Test integral only appears once in file.
        h1e_r, eri_r, enuc_r, nelec_r = read_fcidump('FCIDUMP', symmetry=1,
                                                     verbose=False)
        i,k,j,l = (1,0,0,0)
        ikjl = (i,k,j,l)
        jlik = (j,l,i,k)
        kilj = (k,i,l,j)
        ljki = (l,j,k,i)
        d1 = eri_r[ikjl] - eri_r[kilj].conj()
        d2 = eri_r[ikjl] - eri_r[jlik]
        d3 = eri_r[ikjl] - eri_r[ljki].conj()
        self.assertAlmostEqual(d1,0.0)
        self.assertAlmostEqual(d2,0.0)
        self.assertAlmostEqual(d3,-0.00254428836-0.00238852605j)

    def tearDown(self):
        cwd = os.getcwd()
        files = ['FCIDUMP', 'hamiltonian.h5']
        for f in files:
            try:
                os.remove(cwd+'/'+f)
            except OSError:
                pass


if __name__ == "__main__":
    unittest.main()
