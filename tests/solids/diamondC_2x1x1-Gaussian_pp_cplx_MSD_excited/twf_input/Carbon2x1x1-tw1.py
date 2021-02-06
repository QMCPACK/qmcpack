#! /usr/bin/env python3
import h5py
from pyscf.pbc import gto, scf, df
from numpy import array

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
kpts=array([[ 0.07761248,  0.07761248, -0.07761248],[ 0.54328733,  0.54328733, -0.54328733]])




mydf = df.GDF(cell,kpts)
mydf.auxbasis = 'weigend'
mydf._cderi_to_save = 'df_ints.h5'
mydf.build()                     
mf = scf.KROHF(cell,kpts).density_fit()
mf.exxdiv = 'ewald'
mf.max_cycle = 200
mf.with_df = mydf
mf.chkfile ='diamond-scf.chk'
mf.with_df._cderi = 'df_ints.h5'
#dm=mf.from_chk('diamond-scf.chk')
#e_scf=mf.kernel(dm)
e_scf=mf.kernel()

ener = open('e_scf','w')
ener.write('%s\n' % (e_scf))
print('e_scf',e_scf)
ener.close()

title="C_diamond-twist"
from PyscfToQmcpack import savetoqmcpack
savetoqmcpack(cell,mf,title=title,kmesh=kmesh,kpts=kpts,sp_twist=sp_twist)

from MolPyscfToQPkpts import pyscf2QP2
pyscf2QP2(cell,mf,kpts=kpts,kmesh=kmesh,sp_twist=sp_twist,int_threshold = 1E-15)
