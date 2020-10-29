#! /usr/bin/env python3

import sys
import h5py
import scipy.linalg
import numpy
from pyscf import gto, scf, lib
from afqmctools.utils.pyscf_utils import load_from_pyscf_chk_mol
from afqmctools.utils.linalg import get_ortho_ao_mol
from afqmctools.hamiltonian.mol import write_hamil_mol
from afqmctools.wavefunction.mol import write_qmcpack_wfn

# 1. We will first generate a fake 2 determinant NOMSD trial wavefunction
# expansion made up of the RHF solutions replicated twice. This is nonsense but
# it's just to demonstrate how to write things.
mol = gto.M(atom=[['N', (0,0,0)], ['N', (0,0,3.0)]],
            basis='cc-pvdz',
            unit='Bohr')
nalpha, nbeta = mol.nelec
rhf = scf.RHF(mol)
rhf.chkfile = 'scf.chk'
rhf.kernel()
# Transformation matrix to ortho AO basis.
X = get_ortho_ao_mol(mol.intor('int1e_ovlp_sph'))
Xinv = scipy.linalg.inv(X)
with h5py.File('scf.chk') as fh5:
    fh5['scf/orthoAORot'] = X

# a. Generate Hamiltonian in ortho AO basis.
scf_data = load_from_pyscf_chk_mol('reference/scf_ref.chk')
ortho_ao = True

# b. Fake a two determinant trial wavefunction.
nmo = rhf.mo_coeff.shape[1]
wfn = numpy.zeros((2,nmo,nalpha+nbeta),dtype=numpy.complex128)
Coao = numpy.dot(Xinv, scf_data['mo_coeff'])
wfn[0,:,:nalpha] = Coao[:,:nalpha]
wfn[0,:,nalpha:] = Coao[:,:nbeta]
wfn[1,:,:nalpha] = Coao[:,:nalpha]
wfn[1,:,nalpha:] = Coao[:,:nbeta]
c = numpy.array([1.0/numpy.sqrt(2)+0j])
write_hamil_mol(scf_data, 'afqmc.h5', 1e-5, verbose=True, ortho_ao=ortho_ao)
write_qmcpack_wfn('afqmc.h5', (c,wfn), True, mol.nelec, nmo)
