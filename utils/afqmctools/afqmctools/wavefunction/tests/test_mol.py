import h5py
import numpy
import os
import unittest
from pyscf import gto, scf
from afqmctools.wavefunction import mol

class TestMolWavefunction(unittest.TestCase):

    def test_write_wfn_mol(self):
        atom = gto.M(atom='Ne 0 0 0', basis='sto-3g', verbose=0)
        mf = scf.RHF(atom)
        energy = mf.kernel()
        scf_data = {'mol': atom, 'mo_coeff': mf.mo_coeff, 'mo_occ': mf.mo_occ,
                'X': mf.mo_coeff, 'isUHF': False, 'mo_energy': mf.mo_energy}
        with h5py.File('wfn.h5', 'w') as fh5:
            pass
        mol.write_wfn_mol(scf_data, False, 'wfn.h5')
        with h5py.File('wfn.h5', 'r') as fh5:
            dims = fh5['Wavefunction/NOMSD/dims'][:]
        self.assertTrue(numpy.allclose(dims, [5,5,5,1,1]))

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
