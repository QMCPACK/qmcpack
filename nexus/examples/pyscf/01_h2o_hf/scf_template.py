#! /usr/bin/env python

from pyscf import scf

# Nexus expands this with Mole info
$system

mf = scf.RHF(mol)
mf.kernel()
