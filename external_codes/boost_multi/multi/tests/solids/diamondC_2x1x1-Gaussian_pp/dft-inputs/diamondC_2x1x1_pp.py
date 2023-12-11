#! /usr/bin/env python3

'''
Gamma point post-HF calculation needs only real integrals.
Methods implemented in finite-size system can be directly used here without
any modification.
'''


import numpy
from pyscf.pbc import gto, scf, dft
from pyscf import gto as Mgto
#from mpi4pyscf.pbc import df 
from pyscf.pbc import df 
from pyscf.pbc import ao2mo
from pyscf.pbc import tools
from pyscf.pbc.tools.pbc import super_cell


nmp = [2, 1, 1]

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

cell.verbose = 5


cell.build()


supcell = super_cell(cell, nmp)
mydf = df.FFTDF(supcell)
mydf.auxbasis = 'weigend'
kpts=[]
mf = dft.RKS(supcell)
mf.xc = 'lda'

mf.exxdiv = 'ewald'
mf.with_df = mydf

e_scf=mf.kernel()


ener = open('e_scf','w')
ener.write('%s\n' % (e_scf))
print 'e_scf',e_scf


title="C_Diamond"

from PyscfToQmcpack import savetoqmcpack
savetoqmcpack(supcell,mf,title=title,kpts=kpts)
