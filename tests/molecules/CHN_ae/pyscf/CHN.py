#! /usr/bin/env python3

from pyscf import gto
from pyscf import scf, dft, df

mol = gto.Mole()
mol.verbose = 5
mol.atom =''' 
  C   0.97091888623401     -0.07609978449644     -0.06440023695551
  N   -0.17562187186916     -0.07610010775030     -0.06439988152408
  H   2.03800298563516     -0.07610010775326     -0.06439988152041
 '''
mol.unit='A'
mol.basis = 'cc-pvtz'
#mol.cart=True
mol.build()




mf = dft.RKS(mol).density_fit()
mf.xc ='B3LYP'

mf.chkfile ='CHN.chk'
e_scf=mf.kernel()


title="CHN"
kpts=[]
from PyscfToQmcpack import savetoqmcpack
savetoqmcpack(mol,mf,title=title,kpts=kpts)

