#! /usr/bin/env python3
# Triplet UHF ground state of carbon atom.

from pyscf import gto, scf
import numpy

import h5py

mol = gto.Mole()
mol.basis = 'cc-pvtz'
mol.atom = (('C', 0,0,0),)
mol.spin = 2
mol.verbose = 4
mol.build()

mf = scf.UHF(mol)
mf.chkfile = 'scf.chk'
mf.kernel()
# Check if UHF solution is stable.
mf.stability()
