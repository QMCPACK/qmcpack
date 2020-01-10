#! /usr/bin/env python3

import sys
import os
from pyscf import gto, scf, tools
import numpy

mol = gto.Mole(atom=[['Ne', (0,0,0)]],
            basis='cc-pvdz',
            unit='Angstrom',
            verbose=4)
mol.build()
mf = scf.RHF(mol)
mf.chkfile = 'neon.chk.h5'
mf.kernel()
hcore = mf.get_hcore()
fock = mf.get_veff()
# ascii fcidump is converted using fcidump_to_qmcpack.py.
tools.fcidump.from_scf(mf, 'FCIDUMP')
nmo = mf.mo_coeff.shape[-1]

wfn_header = """&FCI
UHF = 0
CMajor
NCI = 1
TYPE = matrix
/"""
orbmat = numpy.identity(nmo,dtype=numpy.complex128)
with open("wfn_rhf.dat", 'w') as f:
    f.write(wfn_header+'\n')
    f.write('Coefficients: 1.0\n')
    f.write('Determinant: 1\n')
    for i in range(0,mol.nelec[0]):
        occ = '0.0 '*i + '1.0' + ' 0.0'*(nmo-i-1) + '\n'
        f.write(occ)
    # Alternatively add FullMO to the header above and replace for loop above
    # with double loop below.
    # for i in range(0,nmo):
        # for j in range(0,nmo):
            # val = orbmat[i,j]
            # f.write('(%.10e,%.10e) '%(val.real, val.imag))
        # f.write('\n')
