#! /usr/bin/env python3
import numpy
import h5py
from pyscf.pbc import gto, scf, dft, df
from pyscf.pbc import df 

cell = gto.Cell()
cell.a             = '''
         3.37316115       3.37316115       0.00000000
         0.00000000       3.37316115       3.37316115
         3.37316115       0.00000000       3.37316115'''
cell.atom = '''  
   C        0.00000000       0.00000000       0.00000000
   C        1.686580575      1.686580575      1.686580575 
            ''' 
cell.basis         = 'bfd-vdz'
cell.ecp           = 'bfd'
cell.unit          = 'B'
cell.drop_exponent = 0.1
cell.verbose       = 5
cell.charge        = 0
cell.spin          = 0
cell.build()


sp_twist=[0.07761248, 0.07761248, -0.07761248]

kmesh=[2,1,1]
kpts=[[ 0.07761248,  0.07761248, -0.07761248],[ 0.54328733,  0.54328733, -0.54328733]]


mf = scf.KRHF(cell,kpts)
mf.exxdiv = 'ewald'
mf.max_cycle = 200

e_scf=mf.kernel()

ener = open('e_scf','w')
ener.write('%s\n' % (e_scf))
print('e_scf',e_scf)
ener.close()

title="C_diamond-tiled-cplx"
from PyscfToQmcpack import savetoqmcpack
savetoqmcpack(cell,mf,title=title,kmesh=kmesh,kpts=kpts,sp_twist=sp_twist)
