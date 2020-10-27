#! /usr/bin/env python3

import numpy
from functools import reduce
from pyscf.pbc import gto, scf
from pyscf.pbc import tools as pbctools

alat0 = 3.6

cell = gto.Cell()
cell.a = (numpy.ones((3,3))-numpy.eye(3))*alat0/2.0
cell.atom = (('C',0,0,0),('C',numpy.array([0.25,0.25,0.25])*alat0))
cell.basis = 'gth-dzvp'
cell.pseudo = 'gth-pade'
cell.gs = [10]*3
cell.verbose = 5
cell.build()

nk = [1,1,1]
kpts = cell.make_kpts(nk)

mf = scf.KRHF(cell)
mf.chkfile = 'scf.dump'
ehf = mf.kernel()

import h5py

from pyscftools import integrals_from_chkfile
hcore = mf.get_hcore()                                   # obtain and store core hamiltonian
fock = (hcore + mf.get_veff())                           # store fock matrix (required with orthoAO)
X,nmo_per_kpt = integrals_from_chkfile.getOrthoAORotation(cell,kpts,1e-8)      # store rotation to orthogonal PAO basis
with h5py.File(mf.chkfile) as fh5:
  fh5['scf/hcore'] = hcore
  fh5['scf/fock'] = fock
  fh5['scf/orthoAORot'] = X
  fh5['scf/nmo_per_kpt'] = nmo_per_kpt

integrals_from_chkfile.eri_to_h5("choldump", "./scf.dump", orthoAO=False, gtol=1e-5, gtol_chol=1e-5)

