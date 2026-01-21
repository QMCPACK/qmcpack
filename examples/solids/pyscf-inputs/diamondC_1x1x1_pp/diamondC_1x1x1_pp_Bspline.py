#! /usr/bin/env python3

'''
Gamma point post-HF calculation needs only real integrals.
Methods implemented in finite-size system can be directly used here without
any modification.
'''


import numpy as np
from pyscf.pbc import gto, scf, dft
from pyscf import gto as Mgto
#from mpi4pyscf.pbc import df 
from pyscf.pbc import df 
from pyscf.pbc import ao2mo
from pyscf.pbc import tools
from pyscf.pbc.tools.pbc import super_cell

nmp = [1, 1, 1]

cell = gto.Cell()

cell.a = '''
         3.37316115       3.37316115       0.00000000
         0.00000000       3.37316115       3.37316115
         3.37316115       0.00000000       3.37316115'''
cell.atom = '''  
   C        0.00000000       0.00000000       0.00000000
   C        1.686580575      1.686580575      1.686580575 
            ''' 
cell.basis='bfd-vtz'
cell.ecp = 'bfd'


cell.unit='B'
cell.drop_exponent=0.1

cell.ke_cutoff = 120 # Kinetic energy cut-off. Check the convergence!

cell.verbose = 5


cell.build()


supcell = super_cell(cell, nmp)
supcell.mesh = cell.mesh*nmp // 2 * 2 + 1 # Define FFT mesh for supercell. Do not modify!!
mydf = df.FFTDF(supcell)
mydf.auxbasis = 'weigend'
kpts=[]
mf = dft.RKS(supcell)
mf.xc = 'lda'

mf.exxdiv = 'ewald'
mf.with_df = mydf

e_scf=mf.kernel()
assert mf.converged

ener = open('e_scf','w')
ener.write('%s\n' % (e_scf))

title="C_Diamond"

from PyscfToQmcpack_Spline import pyscf2qmcpackspline

pyscf2qmcpackspline(supcell,mf,title=title,kpts=kpts)
