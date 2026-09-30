import h5py
import numpy
import os
import unittest
from pyscf.pbc import gto, scf, dft
from afqmctools.wavefunction import pbc

class TestMolWavefunction(unittest.TestCase):

    def test_write_wfn_pbc(self):
        cell = gto.Cell()
        alat0 = 3.6
        cell.a = (numpy.ones((3,3))-numpy.eye(3))*alat0/2.0
        cell.atom = (('C',0,0,0),('C',numpy.array([0.25,0.25,0.25])*alat0))
        cell.basis = 'gth-szv'
        cell.pseudo = 'gth-pade'
        cell.mesh = [12]*3  # 10 grids on postive x direction, => 21^3 grids in total
        cell.verbose = 0
        cell.build()

        nk = [2,2,1]
        kpts = cell.make_kpts(nk)

        mf = dft.KRKS(cell,kpts=kpts)
        mf.chkfile = 'scf.dump'
        emf = mf.kernel()
        nmo_pk = numpy.array([C.shape[-1] for C in mf.mo_coeff])
        nkpts = len(kpts)
        kpts = kpts
        nmo_max = numpy.max(nmo_pk)
        hcore = mf.get_hcore()
        fock = hcore + mf.get_veff()
        scf_data = {'cell': cell,
                    'mo_coeff': mf.mo_coeff,
                    'Xocc': mf.mo_occ,
                    'X': mf.mo_coeff,
                    'fock': fock,
                    'isUHF': False,
                    'hcore': hcore,
                    'nmo_pk': nmo_pk,
                    'mo_energy': mf.mo_energy,
                    'kpts': kpts}
        with h5py.File('wfn.h5', 'w') as fh5:
            pass
        pbc.write_wfn_pbc(scf_data, True, 'wfn.h5', rediag=True)
        with h5py.File('wfn.h5', 'r') as fh5:
            dims = fh5['Wavefunction/NOMSD/dims'][:]
        print(dims)
        self.assertTrue(numpy.allclose(dims, [32,16,16,1,1]))

    def tearDown(self):
        cwd = os.getcwd()
        files = ['wfn.dat', 'wfn.h5', 'scf.dump']
        for f in files:
            try:
                os.remove(cwd+'/'+f)
            except OSError:
                pass

if __name__ == "__main__":
    unittest.main()
