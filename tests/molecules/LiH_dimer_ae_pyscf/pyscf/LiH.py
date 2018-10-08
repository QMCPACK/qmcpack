#!/usr/bin/env python


# Note import path which is different to molecule code
#from pyscf.pbc import gto, scf, df, dft
from pyscf import gto, scf, df, dft
import  numpy


cell = gto.M(
   atom ='''Li  0.0 0.0 0.0
            H   0.0 0.0 3.0139239778''',
   basis ={'H':'cc-pv5z','Li':gto.basis.parse('''
#BASIS SET: (14s,7p,4d,3f,2g,1h) -> [6s,5p,4d,3f,2g,1h]
Li    S
  29493.0000000              0.0000180             -0.0000030
   4417.1010000              0.0001410             -0.0000220
   1005.2230000              0.0007390             -0.0001150
    284.7009000              0.0031070             -0.0004870
     92.8654300              0.0111350             -0.0017460
     33.5117900              0.0346700             -0.0055200
     13.0418000              0.0921710             -0.0149280
      5.3575360              0.1995760             -0.0342060
      2.2793380              0.3288360             -0.0621550
      0.9939900              0.3459750             -0.0959020
Li    S
      0.4334710              1.0000000
Li    S
      0.0955660              1.0000000
Li    S
      0.0446570              1.0000000
Li    S
      0.0206330              1.0000000
Li    P
     11.2500000              0.0013120
      2.5000000              0.0099180
      0.6500000              0.0375420
Li    P
      0.2500000              1.0000000
Li    P
      0.1000000              1.0000000
Li    P
      0.0390000              1.0000000
Li    P
      0.0170000              1.0000000
Li    D
      0.5500000              1.0000000
Li    D
      0.2900000              1.0000000
Li    D
      0.1400000              1.0000000
Li    D
      0.0610000              1.0000000
Li    F
      0.3500000              1.0000000
Li    F
      0.2200000              1.0000000
Li    F
      0.1100000              1.0000000
Li    G
      0.3200000              1.0000000
Li    G
      0.1600000              1.0000000
''')},
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
