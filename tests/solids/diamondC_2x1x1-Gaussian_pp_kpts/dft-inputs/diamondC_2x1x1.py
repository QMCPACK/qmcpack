#! /usr/bin/env python3

'''
Gamma point post-HF calculation needs only real integrals.
Methods implemented in finite-size system can be directly used here without
any modification.
'''

import numpy as np
from pyscf import lib
from pyscf.pbc import gto, scf, dft
from pyscf import gto as Mgto
from pyscf.pbc import df 
from pyscf.pbc import ao2mo
from pyscf.pbc import tools
from pyscf.pbc.tools.pbc import super_cell
from functools import reduce
import scipy.linalg as la
import os

restart = False


### generated system text ###
from numpy import array
from pyscf.pbc import gto as gto_loc
cell = gto_loc.Cell()
cell.a             = '''
                     3.37316115   3.37316115    0.000000
                     0.000000   3.37316115   3.37316115
                     3.37316115   0.000000   3.37316115
                     '''
cell.dimension     = 3
cell.basis         = 'bfd-vdz'
cell.ecp           = 'bfd'
cell.unit          = 'B'
cell.atom          = '''
                     C     0.000000   0.00000000   0.00000000
                     C    1.68658058   1.68658058   1.68658058
                     '''
cell.drop_exponent = 0.1
cell.verbose       = 5
cell.charge        = 0
cell.spin          = 0
cell.build()


### end generated system text ###

sp_twist=[0.0, 0.0, 0.0]

kpts = [[ 0.0,  0.0, 0.0],
 [0.46567485,  0.46567485, -0.46567485]]

supcell=cell
mydf = df.FFTDF(supcell,kpts)
#mydf = df.GDF(supcell,kpts)
mydf.auxbasis = 'weigend'
mf = scf.KRHF(supcell,kpts).density_fit()
#mf.tol           = 1e-10 

mf.exxdiv = 'ewald'
mf.with_df = mydf
e_scf=mf.kernel()                      # new

ener = open('e_scf','w')
ener.write('%s\n' % (e_scf))
print('e_scf',e_scf)
ener.close()

title="S2-twist1"
### generated conversion text ###
from PyscfToQmcpack import savetoqmcpack
savetoqmcpack(cell,mf,title=title,kpts=kpts,sp_twist=sp_twist)

 


