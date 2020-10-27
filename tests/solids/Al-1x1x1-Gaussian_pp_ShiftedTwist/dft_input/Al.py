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


kmesh = [1, 1, 1]

#cell = gto.M(
#    atom = '''He     0.  0.  0. ''',
cell = gto.Cell()
# .a is a matrix for lattice vectors.





cell.a = '''
  4.03893   0.00000   0.00000
 -0.00000   4.03893   0.00000
 -0.00000  -0.00000   4.03893
'''
cell.atom = '''  
Al   0.00000   0.00000   0.00000
Al   2.01946   2.01946   0.00000
Al   2.01946   0.00000   2.01946
Al   0.00000   2.01946   2.01946
'''
#Basis: Al-ccpvDz             '''
cell.basis= {'Al': Mgto.parse('''
Al S
8.257944 0.003287
4.514245 -0.017168
2.467734 0.069766
1.348998 -0.183475
0.737436 -0.147133
0.403123 0.046882
0.220369 0.308423
0.120466 0.451564
0.065853 0.302904
0.035999 0.079545
Al S
0.236926 1.000000
Al P
1.570603 -0.002645
0.977752 -0.037850
0.608683 0.006636
0.378925 0.089291
0.235893 0.134421
0.146851 0.256105
0.091420 0.238970
0.056912 0.260677
0.035429 0.112350
0.022056 0.052665
Al P
0.202698 1.000000
Al D
0.192882 1.000000
''')
}

cell.ecp = {'Al': Mgto.basis.parse_ecp('''
Al nelec 10
Al ul
1 5.073893 3.000000
3 8.607001 15.221680
2 3.027490 -11.165685
Al S
2 7.863954 14.879513
2 2.061358 20.746863
Al P
2 3.125175 7.786227
2 1.414930 7.109015
''')
}



cell.unit='A'
cell.drop_exponent=0.1

cell.verbose = 5
cell.spin =0

cell.build()
sp_twist=[0.11,0.23,-0.34]
twist = np.asarray(sp_twist) / 1.0
kmesh=[1,1,1]
kpts = cell.make_kpts((1,1,1), with_gamma_point=False,  wrap_around=True, scaled_center=twist)

supcell=cell
mydf = df.GDF(supcell,kpts)
mydf.auxbasis = 'weigend'
mydf._cderi_to_save = 'df_ints.h5' 
mydf.build()         
mf = scf.KRHF(supcell,kpts).density_fit()



mf.with_df._cderi = 'df_ints.h5'
mf.exxdiv = 'ewald'
mf.with_df = mydf
mf.chkfile ='Al-TZ.chk'

e_scf=mf.kernel()    


ener = open('e_scf','w')
ener.write('%s\n' % (e_scf))
print('e_scf',e_scf)
ener.close()

title="Al-DZ"

from PyscfToQmcpack import savetoqmcpack
savetoqmcpack(cell,mf,title=title,kmesh=kmesh,kpts=kpts,sp_twist=kpts)



