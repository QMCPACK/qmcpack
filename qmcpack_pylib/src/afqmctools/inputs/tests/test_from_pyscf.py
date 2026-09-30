import os
import unittest
import h5py
from pyscf import gto, scf
from afqmctools.hamiltonian.mol import generate_hamiltonian
from afqmctools.inputs.from_pyscf import write_qmcpack
from afqmctools.inputs.energy import calculate_hf_energy
from afqmctools.utils.pyscf_utils import load_from_pyscf_chk_mol
from afqmctools.utils.linalg import get_ortho_ao_mol

class TestEnergy(unittest.TestCase):

    def test_from_pyscf(self):
        atom = gto.M(atom='Ne 0 0 0', basis='sto-3g', verbose=0)
        mf = scf.RHF(atom)
        mf.chkfile = 'scf.chk'
        energy = mf.kernel()
        self.assertAlmostEqual(energy, -126.60452499805)
        write_qmcpack('scf.chk', 'afqmc.h5', 1e-8, wfn_file='afqmc.h5',
                      dense=True, real_chol=True)
        etot, e1b, e2b = calculate_hf_energy('afqmc.h5', 'afqmc.h5')
        self.assertAlmostEqual(etot, -126.60452499805)

    def test_from_pyscf_lindep(self):
        atom = gto.M(atom='Ne 0 0 0', basis='aug-ccpvtz', verbose=0)
        from pyscf.scf.addons import remove_linear_dep_
        mf = scf.RHF(atom)
        remove_linear_dep_(mf, 0.1, 0.1)
        mf.chkfile = 'scf.chk'
        energy = mf.kernel()
        X = get_ortho_ao_mol(atom.intor('int1e_ovlp_sph'), 0.1)
        with h5py.File('scf.chk', 'r+') as fh5:
            fh5['scf/orthoAORot'] = X
        self.assertAlmostEqual(energy,  -128.11089429105502)
        write_qmcpack('scf.chk', 'afqmc.h5', 1e-8, wfn_file='afqmc.h5',
                      dense=True, real_chol=True, ortho_ao=True)
        etot, e1b, e2b = calculate_hf_energy('afqmc.h5', 'afqmc.h5')
        self.assertAlmostEqual(etot, -128.11089429105502)

    def tearDown(self):
        cwd = os.getcwd()
        files = ['afqmc.h5', 'scf.chk']
        for f in files:
            try:
                os.remove(cwd+'/'+f)
            except OSError:
                pass
