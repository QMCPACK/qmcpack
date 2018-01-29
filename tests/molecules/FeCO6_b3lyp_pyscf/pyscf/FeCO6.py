#!/usr/bin/env python

'''
Gamma point Hartree-Fock/DFT
The 2-electron integrals are computed using Poisson solver with FFT by default.
In most scenario, it should be used with pseudo potential.
'''

# Note import path which is different to molecule code
#from pyscf.pbc import gto, scf, df, dft
from pyscf import gto, dft
import  numpy





cell = gto.M(
   atom ='''
             Fe    0.0000013   -0.0000001   -0.0000001
             C     0.0000037   -0.0072194    2.3405884
             O    -0.0000299   -0.0415619    3.4549199
             C     2.3029923   -0.0000137   -0.0000135
             O     3.4182461   -0.0000170   -0.0000168
             C    -2.3029877    0.0000133    0.0000131
             O    -3.4182415    0.0000167    0.0000165
             C    -0.0000019    0.0072196   -2.3405883
             O     0.0000250    0.0415622   -3.4549197
             C    -0.0000028   -2.3369861    0.0065791
             O     0.0000240   -3.4513475    0.0412328
             C     0.0000037    2.3369862   -0.0065788
             O    -0.0000326    3.4513476   -0.0412324
''',
   basis ='bfd-vtz',
   ecp='bfd',
   unit="angstrom",
   spin=4,
   verbose = 5,
   cart=True,
)




#mf = dft.RKS(cell).density_fit()
mf = dft.RKS(cell)
mf.xc = 'b3lyp'
mf.kernel()

title='FeCO6'


from PyscfToQmcpack import savetoqmcpack
savetoqmcpack(cell,mf,title)
