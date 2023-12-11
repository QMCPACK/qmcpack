#! /usr/bin/env python3

import sys
import h5py
import numpy
from pyscf import gto, scf, mcscf, fci, lib
from afqmctools.utils.pyscf_utils import load_from_pyscf_chk_mol
from afqmctools.hamiltonian.mol import write_hamil_mol
from afqmctools.wavefunction.mol import write_qmcpack_wfn

def gen_wavefunction_and_hamil(tol=0.02):
    mol = gto.M(atom=[['N', (0,0,0)], ['N', (0,0,3.0)]],
                basis='cc-pvdz',
                unit='Bohr')
    nalpha, nbeta = mol.nelec
    rhf = scf.RHF(mol)
    rhf.chkfile = 'scf.chk'
    rhf.kernel()
    # 2. Write CASSCF wavefunction.
    # a. Perform CASSCF calculation, this replicates the calculations from
    # J. Chem. Phys. 127, 144101 (2007).
    # They find a CASSCF energy of -108.916484 Ha, and a ph-AFQMC energy of
    # -109.1975(6) Ha with a 97 determinant CASSCF trial.
    M = 12
    N = 6
    nmo = rhf.mo_coeff.shape[-1]
    mc = mcscf.CASSCF(rhf, M, N)
    mc.chkfile = 'scf.chk'
    mc.kernel()
    nalpha = 3
    nbeta = 3
    # Extract ci expansion in the form of a tuple: (coeff, occ_a, occ_b).
    # Note the tol param which will return wavefunction elements with abs(ci) > tol.
    ci, occa, occb = zip(*fci.addons.large_ci(mc.ci, M, (nalpha,nbeta),
                         tol=tol, return_strs=False))
    # Reinsert frozen core.
    core = [i for i in range(mc.ncore)]
    occa = [numpy.array(core + [o + mc.ncore for o in oa]) for oa in occa]
    occb = [numpy.array(core + [o + mc.ncore for o in ob]) for ob in occb]
    # b. Generate Hamiltonian.
    # Passing 'mccsf' will tell helper function to read mo_coeff/mo_occ from mcscf
    # group in checkpoint file. We will rotate the integrals by mc.mo_coeff.
    scf_data = load_from_pyscf_chk_mol('scf.chk', 'mcscf')
    write_hamil_mol(scf_data, 'afqmc.h5', 1e-5, verbose=True)
    ci = numpy.array(ci, dtype=numpy.complex128)
    # Save CI expansion to scf file. Not necessary for AFQMC but for test
    # reproducibility.
    # try:
        # with h5py.File('scf.chk', 'r+') as fh5:
            # fh5['ci'] = ci
            # fh5['occa'] = occa
            # fh5['occb'] = occb
    # except RuntimeError:
        # print("Reference data already exists in file")
    uhf = True # UHF always true for CI expansions.
    write_qmcpack_wfn('afqmc.h5', (ci, occa, occb), uhf, mol.nelec, nmo)


def gen_test_data():
    scf_data = load_from_pyscf_chk_mol('reference/scf.chk', 'mcscf')
    write_hamil_mol(scf_data, 'afqmc.h5', 1e-5, verbose=True)
    # Save CI expansion to scf file. Not necessary for AFQMC but for test
    # reproducibility.
    with h5py.File('reference/scf.chk', 'r') as fh5:
        ci = fh5['ci'][:]
        occa = fh5['occa'][:]
        occb = fh5['occb'][:]
    uhf = True # UHF always true for CI expansions.
    write_qmcpack_wfn('afqmc.h5', (ci, occa, occb), uhf, (7,7), 28)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        gen_test_data()
    else:
        gen_wavefunction_and_hamil()
