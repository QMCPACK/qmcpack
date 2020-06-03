#! /usr/bin/env python3

import numpy
from pyscf.pbc import gto, scf
from pyscf.pbc import tools as pbctools
import h5py
from pyscftools.integrals_from_chkfile import eri_to_h5
from pyscftools.integrals_from_chkfile_kpfftdf import (
        eri_to_h5, getOrthoAORotation
)


alat0 = 3.6

cell = gto.Cell()
cell.a = numpy.eye(3)*alat0
cell.atom = [('He',0,0,0)]
cell.basis = 'gth-dzv'
cell.pseudo = 'gth-pade'
cell.mesh = [71]*3
cell.verbose = 5
cell.build()

nk = [2,2,2]
kpts = cell.make_kpts(nk)

mf = scf.KRHF(cell,kpts=kpts)
mf.chkfile = 'scf.dump'
ehf = mf.kernel()

hcore = mf.get_hcore()
fock = (hcore + mf.get_veff())
X,nmo_per_kpt = getOrthoAORotation(cell,kpts,1e-8)
with h5py.File(mf.chkfile) as fh5:
  fh5['scf/hcore'] = hcore
  fh5['scf/fock'] = fock
  fh5['scf/orthoAORot'] = X
  fh5['scf/nmo_per_kpt'] = nmo_per_kpt

eri_to_h5("ham_chol_uc", "./scf.dump", orthoAO=True, gtol=1e-5)
