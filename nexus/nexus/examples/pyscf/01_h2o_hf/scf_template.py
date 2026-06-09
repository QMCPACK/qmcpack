#! /usr/bin/env python3

from pyscf import scf

# Nexus expands this with Mole info
$system

mf = scf.RHF(mol)
mf.kernel()
