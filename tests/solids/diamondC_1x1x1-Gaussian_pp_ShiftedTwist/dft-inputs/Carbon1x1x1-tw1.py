#!/usr/bin/env python

from pyscf.pbc import gto, scf, df

cell = gto.Cell()
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

sp_twist=[0.333, 0.333, 0.333]
kpts = [[ 0.333,  0.333, 0.333]]

supcell=cell
mydf = df.FFTDF(supcell,kpts)
mydf.auxbasis = 'weigend'

mf = scf.KRHF(supcell,kpts)
mf.tol           = 1e-10 

mf.exxdiv = 'ewald'
mf.with_df = mydf
mf.chkfile ='diamond-scf.chk'
mf.with_df._cderi = 'df_ints.h5'
e_scf=mf.kernel()                      

ener = open('e_scf','w')
ener.write('%s\n' % (e_scf))
print('e_scf',e_scf)
ener.close()

title="S1-twist1"
from PyscfToQmcpack import savetoqmcpack
savetoqmcpack(cell,mf,title=title,kpts=kpts,sp_twist=sp_twist)

 


