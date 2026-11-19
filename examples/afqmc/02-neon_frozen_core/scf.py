#! /usr/bin/env python3

from pyscf import cc, gto, scf
from pyscf.cc import ccsd_t

mol = gto.Mole()
mol.basis = "aug-cc-pvdz"
mol.atom = (("Ne", 0, 0, 0),)
mol.verbose = 4
mol.build()

mf = scf.RHF(mol)
mf.chkfile = "scf.chk"
ehf = mf.kernel()

ccsd = cc.CCSD(mf)
eccsd = ccsd.kernel()[0]
ecorr_ccsdt = ccsd_t.kernel(ccsd, ccsd.ao2mo())
print(f"E(CCSD(T)) = {ehf + eccsd + ecorr_ccsdt}")
