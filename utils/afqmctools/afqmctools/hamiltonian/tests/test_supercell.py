import h5py
import numpy
import os
import scipy.linalg
import scipy.sparse
import unittest
from pyscf.pbc import gto, dft, df, tools
import afqmctools.hamiltonian.supercell as sc
from afqmctools.utils.linalg import get_ortho_ao
try:
    from mpi4py import MPI
    no_mpi = False
except ImportError:
    no_mpi = True

class FakeComm:
    def __init__(self, size=1):
        self.size = size
        self.rank = 0

class TestSupercell(unittest.TestCase):

    def setUp(self):
        cell = gto.Cell()
        alat0 = 3.6
        cell.a = (numpy.ones((3,3))-numpy.eye(3))*alat0/2.0
        cell.atom = (('C',0,0,0),('C',numpy.array([0.25,0.25,0.25])*alat0))
        cell.basis = 'gth-szv'
        cell.pseudo = 'gth-pade'
        cell.mesh = [12]*3  # 10 grids on postive x direction, => 21^3 grids in total
        cell.verbose = 0
        cell.build()
        self.cell = cell

        nk = [2,1,1]
        kpts = cell.make_kpts(nk)

        mf = dft.KRKS(cell,kpts=kpts)
        mf.chkfile = 'scf.dump'
        ehf = mf.kernel()
        self.mf = mf
        self.nmo_pk = numpy.array([C.shape[-1] for C in mf.mo_coeff])
        self.nkpts = len(kpts)
        self.kpts = kpts
        self.nmo_max = numpy.max(self.nmo_pk)

    def test_basis_map(self):
        mf = self.mf
        Xocc = mf.mo_occ
        ik2n, nmo_tot = sc.setup_basis_map(Xocc, self.nmo_max, self.nkpts,
                                           self.nmo_pk, False)
        # Original mapping was broken for MO basis.
        # bench = [140, 172, 204, 236, 268, 300, 332, 364]
        # FDM: Update for faster unit tests. bench = [28, 92, 156, 220, 284, 348, 412, 476]
        bench = [28, 92]
        self.assertTrue(numpy.allclose(sum(ik2n), bench))
        self.assertEqual(nmo_tot, 16)

    def test_hcore(self):
        mf = self.mf
        Xocc = mf.mo_occ
        hcore = mf.get_hcore()
        h5f = h5py.File('ham.h5', 'w')
        h5grp = h5f.create_group("Hamiltonian")
        ik2n, nmo_tot = sc.setup_basis_map(Xocc, self.nmo_max,
                                           self.nkpts, self.nmo_pk, False)
        sc.write_one_body(hcore, mf.mo_coeff, self.nkpts,
                          self.nmo_pk, self.nmo_max, ik2n, h5grp, gtol=1e-16)
        # TODO test writing.

    def test_grid_shifts(self):
        gmap, Qi, ngs = sc.generate_grid_shifts(self.cell)
        self.assertEqual(gmap[0,-1], 1570)
        self.assertTrue(numpy.allclose(Qi[1], [0,0,-1.8471769]))
        self.assertEqual(ngs,1728)

    def test_gen_partition(self):
        comm = FakeComm(16)
        nmo_tot = numpy.sum(self.nmo_pk)
        maxvecs = 20 * nmo_tot
        part = sc.Partition(comm, maxvecs, nmo_tot,
                            self.nmo_max, self.nkpts)
        self.assertEqual(part.ij0, 0)
        self.assertEqual(part.ijN, 8)
        self.assertEqual(part.nij, 8)
        self.assertEqual(part.nkk, 2)
        self.assertEqual(part.n2k1[0], 0)
        self.assertEqual(part.n2k2[1], 1)
        self.assertEqual(part.kk0, 0)
        self.assertEqual(part.kkN, 2)

    def test_gen_orbital_products(self):
        cell = self.cell
        mydf = df.FFTDF(cell, self.kpts)
        comm = FakeComm(16)
        nmo_tot = numpy.sum(self.nmo_pk)
        maxvecs = 20 * nmo_tot
        part = sc.Partition(comm, maxvecs, nmo_tot,
                            self.nmo_max, self.nkpts)
        gmap, Qi, ngs = sc.generate_grid_shifts(self.cell)
        X = [numpy.identity(C.shape[0]) for C in self.mf.mo_coeff]
        xik, xlj = sc.gen_orbital_products(cell, mydf, X,
                                           self.nmo_pk, ngs, part, self.kpts,
                                           self.nmo_max)
        self.assertEqual(xik.shape,(2,1728,8))
        self.assertAlmostEqual(numpy.max(numpy.abs(xik)),
                               0.01239256218954048, places=8)
        self.assertAlmostEqual(numpy.max(numpy.abs(xlj)),
                               50.73269574297129, places=8)

    @unittest.skipIf(no_mpi, "MPI4PY not found")
    def test_modified_cholesky(self):
        cell = self.cell
        mydf = df.FFTDF(cell, self.kpts)
        comm = MPI.COMM_WORLD
        nmo_tot = numpy.sum(self.nmo_pk)
        maxvecs = 20 * nmo_tot
        part = sc.Partition(comm, maxvecs, nmo_tot,
                            self.nmo_max, self.nkpts)
        gmap, Qi, ngs = sc.generate_grid_shifts(self.cell)
        X = [numpy.identity(C.shape[0]) for C in self.mf.mo_coeff]
        xik, xlj = sc.gen_orbital_products(cell, mydf, X,
                                           self.nmo_pk, ngs, part, self.kpts,
                                           self.nmo_max)
        kconserv = tools.get_kconserv(cell, self.kpts)
        solver = sc.Cholesky(part, kconserv, 1e-3, verbose=False)
        chol = solver.run(comm, xik, xlj, part, self.kpts,
                          self.nmo_pk, self.nmo_max, Qi, gmap)
        self.assertEqual(chol.shape, (4,64,53))
        self.assertAlmostEqual(numpy.max(numpy.abs(chol)),
                               0.8077979869286759, places=8)
        self.assertEqual(len(chol[abs(chol)>1e-10]), 6047)

    def tearDown(self):
        cwd = os.getcwd()
        files = ['ham.h5', 'scf.dump', 'hamil.h5', 'test.h5']
        for f in files:
            try:
                os.remove(cwd+'/'+f)
            except OSError:
                pass


if __name__ == "__main__":
    unittest.main()
