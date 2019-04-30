import h5py
import numpy
import os
import unittest
from pyscf.pbc import gto, dft, df, tools
import afqmctools.hamiltonian.kpoint as kp
from afqmctools.utils.linalg import get_ortho_ao
try:
    from mpi4py import MPI
    no_mpi = False
except ImportError:
    no_mpi = True
try:
    test_phdf = h5py.h5p.PropFAID.set_fapl_mpio
    phdf = True
except AttributeError:
    phdf = False

cell = gto.Cell()
alat0 = 3.6
cell.a = (numpy.ones((3,3))-numpy.eye(3))*alat0/2.0
cell.atom = (('C',0,0,0),('C',numpy.array([0.25,0.25,0.25])*alat0))
cell.basis = 'gth-szv'
cell.pseudo = 'gth-pade'
cell.mesh = [25]*3  # 10 grids on postive x direction, => 21^3 grids in total
cell.verbose = 0
cell.build()
cell = cell

nk = [2,2,2]
kpts = cell.make_kpts(nk)

mf = dft.KRKS(cell,kpts=kpts)
mf.chkfile = 'scf.dump'
ehf = mf.kernel()
mf = mf
nmo_pk = numpy.array([C.shape[-1] for C in mf.mo_coeff])
nkpts = len(kpts)
kpts = kpts
nmo_max = numpy.max(nmo_pk)

class TestKpoint(unittest.TestCase):

    def test_momentum_maps(self):
        qk, km = kp.construct_qk_maps(cell, kpts)
        self.assertEqual(qk.shape, (8,8))
        self.assertEqual(km.shape, (8,))
        self.assertEqual(qk[1,5], 4)
        self.assertEqual(qk[6,3], 5)
        self.assertEqual(km[3], 3)

    @unittest.skipIf(not phdf, "Parallel hdf5 not found.")
    def test_file_handler(self):
        comm = MPI.COMM_WORLD
        handler = kp.FileHandler(comm, 'ham.h5')
        self.assertFalse(handler.error)
        handler.close()

    @unittest.skipIf(not phdf, "Parallel hdf5 not found.")
    def test_write_basic(self):
        comm = MPI.COMM_WORLD
        hcore = mf.get_hcore()
        X, nmo_pk = get_ortho_ao(cell, kpts)
        # Fake some data
        qk = numpy.zeros((8,8))
        km = numpy.ones((8,))
        h5file = kp.FileHandler(comm, 'ham.h5')
        kp.write_basic(comm, cell, kpts, hcore, h5file,
                       X, nmo_pk, qk, km)
        hcore_k3 = h5file.grp['H1_kp3'][:]
        self.assertEqual(hcore_k3.shape, (8,8,2))
        self.assertAlmostEqual(numpy.max(hcore_k3), 1.27436250577, places=8)
        h5file.close()


class TestKPCholesky(unittest.TestCase):

    @unittest.skipIf(no_mpi, "MPI4PY not found")
    def setUp(self):
        qk, km = kp.construct_qk_maps(cell, kpts)
        self.X = [numpy.identity(C.shape[0]) for C in mf.mo_coeff]
        self.comm = MPI.COMM_WORLD
        self.chol = kp.KPCholesky(self.comm, cell, kpts, 20, nmo_pk, qk, km,
                                  verbose=False, gtol_chol=1e-4)

    @unittest.skipIf(no_mpi, "MPI4PY not found")
    def test_kpchol_constructor(self):
        self.assertEqual(self.chol.maxvecs, 160)
        self.assertEqual(self.chol.part.kkbounds[1], 8)
        self.assertEqual(self.chol.part.kkN, 8)
        self.assertEqual(self.chol.ngs, 15625)
        self.assertEqual(self.chol.gmap.shape, (27,15625))
        self.assertEqual(sum(self.chol.gmap[24]), 122062500)
        qtest = numpy.array([0.92358847, -0.92358847, -0.92358847])
        self.assertTrue(numpy.allclose(self.chol.Qi[4], qtest))
        self.assertEqual(self.chol.maxres_buff.shape, (5,))

    @unittest.skipIf(no_mpi, "MPI4PY not found")
    def test_orbital_products(self):
        part = self.chol.part
        ngs = self.chol.ngs
        Xaoik = numpy.zeros((part.nkk,ngs,part.nij),
                             dtype=numpy.complex128)
        Xaolj = numpy.zeros((part.nkk,ngs,part.nij),
                            dtype=numpy.complex128)
        self.chol.generate_orbital_products(3, self.X, Xaoik, Xaolj)
        self.assertAlmostEqual(numpy.max(numpy.abs(Xaoik)),
                               0.001189761073068490, places=8)
        self.assertAlmostEqual(numpy.max(numpy.abs(Xaolj)),
                               250.497175057856, places=8)

    # TODO isolate writing. This is extremely slow through unittest for some
    # reason.
    # @unittest.skipIf(not phdf, "MPI4PY not found")
    # def test_kpchol_solution(self):
        # h5file = kp.FileHandler(self.comm, 'test.h5')
        # self.chol.run(self.comm, self.X, h5file)
        # h5file.close()
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
