import os
import unittest
from pyscf import gto, scf
from afqmctools.hamiltonian.mol import generate_hamiltonian
from afqmctools.inputs.from_pyscf import write_qmcpack
from afqmctools.inputs.energy import calculate_hf_energy
from afqmctools.utils.pyscf_utils import load_from_pyscf_chk_mol

class TestEnergy(unittest.TestCase):

    def test_from_pyscf(self):
        atom = gto.M(atom='Ne 0 0 0', basis='sto-3g', verbose=0)
        mf = scf.RHF(atom)
        mf.chkfile = 'scf.chk'
        energy = mf.kernel()
        self.assertAlmostEqual(energy, -126.60452499805)
        write_qmcpack('scf.chk', 'afqmc.h5', 1e-5, wfn_file='afqmc.h5',
                      dense=True, real_chol=True, verbose=2)
        assert False
        # etot, e1b, e2b = calculate_hf_energy('afqmc.h5', 'afqmc.h5')
        # self.assertAlmostEqual(etot, -126.60452499805)

    def tearDown(self):
        cwd = os.getcwd()
        files = ['afqmc.h5', 'scf.chk']
        for f in files:
            try:
                os.remove(cwd+'/'+f)
            except OSError:
                pass
