#! /usr/bin/env python3

from ..pyscf import scf

$system

mf = scf.RHF(mol)
mf.kernel()
