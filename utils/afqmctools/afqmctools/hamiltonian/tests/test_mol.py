import numpy
import os
import unittest
from pyscf import gto, ao2mo, scf, mcscf
from afqmctools.hamiltonian.converter import read_qmcpack_sparse
import afqmctools.hamiltonian.mol as mol
from afqmctools.utils.linalg import modified_cholesky_direct
from afqmctools.utils.testing import generate_hamiltonian

class TestMol(unittest.TestCase):

    def test_cholesky_direct(self):
        atom = gto.M(atom='Ne 0 0 0', basis='sto-3g')
        eri = atom.intor('int2e', aosym='s1')
        self.assertEqual(eri.shape,(5,5,5,5))
        eri = eri.reshape(25,25)
        chol = modified_cholesky_direct(eri, 1e-5, cmax=20)
        Mrecon = numpy.dot(chol.T, chol)
        self.assertTrue(numpy.allclose(Mrecon, eri, atol=1e-8, rtol=1e-6))

    def test_chunked_cholesky(self):
        atom = gto.M(atom='Ne 0 0 0', basis='aug-ccpvdz')
        eri = atom.intor('int2e', aosym='s1')
        self.assertEqual(eri.shape,(23,23,23,23))
        eri = eri.reshape(529,529)
        chol = mol.chunked_cholesky(atom, max_error=1e-5)
        self.assertEqual(chol.shape,(122,529))
        eri_loc = numpy.dot(chol.T, chol)
        self.assertTrue(numpy.allclose(eri, eri_loc, atol=1e-5, rtol=1e-3))

    def test_ao2mo_chol(self):
        atom = gto.M(atom='Ne 0 0 0', basis='sto-3g', verbose=0)
        eri = atom.intor('int2e', aosym='s1')
        mf = scf.RHF(atom)
        energy = mf.kernel()
        self.assertAlmostEqual(energy, -126.60452499805)
        self.assertEqual(eri.shape,(5,5,5,5))
        eri = eri.reshape(25,25)
        chol = mol.chunked_cholesky(atom, max_error=1e-5)
        mol.ao2mo_chol(chol, mf.mo_coeff)
        self.assertAlmostEqual(numpy.linalg.norm(chol), 3.52947146946)

    def test_ao2mo_chol_rect(self):
        numpy.random.seed(7)
        h1e, chol, enuc, eri = generate_hamiltonian(17, 4)
        X = numpy.random.random((17,15))
        nchol = chol.shape[0]
        chol_ = mol.ao2mo_chol(chol, X)
        self.assertEqual(chol_.shape, (nchol, 15*15))

    def test_frozen_core(self):
        atom = gto.M(atom='Ne 0 0 0', basis='sto-3g', verbose=0)
        mf = scf.RHF(atom)
        energy = mf.kernel()
        C = mf.mo_coeff
        chol = mol.chunked_cholesky(atom, max_error=1e-5)
        mol.ao2mo_chol(chol, C)
        hcore = mf.get_hcore()
        h1e = numpy.dot(C.T, numpy.dot(hcore, C))
        ncore = 1
        ncas = 4
        h1e, chol, efzc = mol.freeze_core(h1e, chol, 0, 1, 4, verbose=False)
        self.assertEqual(h1e.shape, (2,4,4))
        self.assertEqual(chol.shape, (15,16))
        # Check from CASSCF object with same core.
        mc = mcscf.CASSCF(mf, 4, (4,4))
        h1eff, ecore = mc.get_h1eff()
        self.assertAlmostEqual(efzc, ecore)
        self.assertTrue(numpy.allclose(h1eff, h1e, atol=1e-8, rtol=1e-5))

    def test_generate_hamiltonian(self):
        atom = gto.M(atom='Ne 0 0 0', basis='sto-3g', verbose=0)
        mf = scf.RHF(atom)
        energy = mf.kernel()
        C = mf.mo_coeff
        scf_data = {'mo_coeff': C, 'mol': atom, 'hcore': mf.get_hcore(),
                    'isUHF': False}
        mol.generate_hamiltonian(scf_data)

    def test_write_hamil_mol(self):
        atom = gto.M(atom='Ne 0 0 0', basis='sto-3g', verbose=0)
        mf = scf.RHF(atom)
        energy = mf.kernel()
        C = mf.mo_coeff
        scf_data = {'mo_coeff': C, 'mol': atom, 'hcore': mf.get_hcore(),
                    'isUHF': False}
        h1e, chol, nelec, enuc, X = mol.generate_hamiltonian(scf_data)
        mol.write_hamil_mol(scf_data, 'ham.h5', 1e-5, verbose=False)
        h1e_f, chol_f, enuc_f, nmo, ne = read_qmcpack_sparse('ham.h5')
        self.assertTrue(numpy.allclose(h1e_f, h1e, atol=1e-12, rtol=1e-8))
        chol_f = chol_f.toarray().real.T
        self.assertTrue(numpy.allclose(chol_f, chol, atol=1e-12, rtol=1e-8))
        self.assertAlmostEqual(enuc, enuc_f)

    def tearDown(self):
        cwd = os.getcwd()
        files = ['ham.h5']
        for f in files:
            try:
                os.remove(cwd+'/'+f)
            except OSError:
                pass

if __name__ == "__main__":
    unittest.main()
