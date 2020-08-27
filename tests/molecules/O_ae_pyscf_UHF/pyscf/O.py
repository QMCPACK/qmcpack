#! /usr/bin/env python3


# Note import path which is different to molecule code
#from pyscf.pbc import gto, scf, df, dft
from pyscf import gto, scf, df, dft
import  numpy


cell = gto.M(
   atom ='''O  0.0 0.0 0.0''',
   basis ='cc-pcvdz',
   unit ='A',
   spin=2,
   verbose = 5,
   cart=False,
)



mf = scf.UHF(cell)
mf.kernel()

title='O-UHF-Triplet'


from PyscfToQmcpack import savetoqmcpack
savetoqmcpack(cell,mf,title)
