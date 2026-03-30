#! /usr/bin/env python3
# Triplet UHF ground state of carbon atom.

from pyscf import gto, scf
import numpy

import h5py

$system

mf = scf.UHF(mol)
mf.chkfile = $chkfile
mf.kernel()
# Check if UHF solution is stable.
mf.stability()
