#! /usr/bin/env python3


# Note import path which is different to molecule code
#from pyscf.pbc import gto, scf, df, dft
from pyscf import gto, scf, df, dft, mcscf
import  numpy


cell = gto.M(
   atom ='''Li  0.0 0.0 0.0
            H   0.0 0.0 3.0139239778''',
   basis ='cc-pvdz',
   unit="bohr",
   spin=0,
   verbose = 5,
   cart=False,
)



mf = scf.RHF(cell).run()
mycas = mcscf.CASCI(mf, 4, 2).run()

title='LiH'


from PyscfToQmcpack import savetoqmcpack
savetoqmcpack(cell,mycas,title)
