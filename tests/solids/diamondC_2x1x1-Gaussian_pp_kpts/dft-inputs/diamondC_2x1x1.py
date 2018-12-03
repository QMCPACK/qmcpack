#!/usr/bin/env python

import numpy
from pyscf.pbc import gto, scf, dft
from pyscf import gto as Mgto
from pyscf.pbc import df 
from pyscf.pbc import tools
from pyscf.pbc.tools.pbc import super_cell


kmesh = [2, 1, 1]

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



kpts = cell.make_kpts(kmesh)
kpts -= kpts[0]

supcell=cell
mydf = df.GDF(supcell,kpts)
mydf.auxbasis = 'weigend'
mf = scf.KRHF(supcell,kpts).density_fit()

mf.exxdiv = 'ewald'
mf.with_df = mydf
e_scf=mf.kernel()


title="C_Diamond-211"

from PyscfToQmcpack import savetoqmcpack
savetoqmcpack(supcell,mf,title=title,kpts=kpts,kmesh=kmesh)

