#! /usr/bin/env python3

from pyscf import scf

# Nexus expands this with Mole info

### generated system text ###
from pyscf import gto as gto_loc
mol = gto_loc.Mole()
mol.atom     = '''
               O   0.00000000   0.00000000   0.00000000
               H   0.00000000   0.75716000   0.58626000
               H   0.00000000   0.75716000  -0.58626000
               '''
mol.basis    = 'ccpvtz'
mol.unit     = 'A'
mol.charge   = 0
mol.spin     = 0
mol.symmetry = True
mol.build()
### end generated system text ###



mf = scf.RHF(mol)
mf.kernel()
