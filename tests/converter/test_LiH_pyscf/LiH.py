#!/usr/bin/env python


# Note import path which is different to molecule code
#from pyscf.pbc import gto, scf, df, dft
from pyscf import gto, scf, df, dft
import  numpy


cell = gto.M(
   atom ='''Li  0.0 0.0 0.0
            H   0.0 0.0 3.0139239778''',
   basis ='cc-pv5z',
   unit="bohr",
   spin=0,
   verbose = 5,
   cart=False,
)



mf = scf.ROHF(cell)
mf.kernel()

title='LiH'


from PyscfToQmcpack import savetoqmcpack
savetoqmcpack(cell,mf,title)
