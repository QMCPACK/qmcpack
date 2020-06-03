#! /usr/bin/env python3

import sys
import os
from pyscf.pbc import scf, gto, tools
import numpy
import time
import h5py

alat0 = 3.6
nks = 1
ngs = 12

cell = gto.Cell()
cell.a = (numpy.ones((3,3))-numpy.eye(3))*alat0/2.0
cell.atom = (('C',0,0,0),('C',numpy.array([0.25,0.25,0.25])*alat0))
cell.basis = 'gth-szv'
cell.pseudo = 'gth-pade'
cell.gs = [ngs]*3  # 10 grids on postive x direction, => 21^3 grids in total
cell.verbose = 5
cell.build()

nk = numpy.array([nks,nks,nks])
kpt = cell.make_kpts(nk)
mf = scf.KRHF(cell,kpt)
mf.chkfile = "scf.dump"
ehf = mf.kernel()

from qmctools import orthoAO
hcore = mf.get_hcore(kpts=kpt)                                   # obtain and store core hamiltonian
fock = (hcore + mf.get_veff(kpts=kpt))                           # store fock matrix (required with orthoAO)
LINDEP = 1e-8
X,nmo_per_kpt = orthoAO.getOrthoAORotation(cell,kpt,LINDEP)      # store rotation to orthogonal PAO basis
with h5py.File(mf.chkfile) as fh5:
  fh5['scf/hcore'] = hcore
  fh5['scf/fock'] = fock
  fh5['scf/orthoAORot'] = X
  fh5['scf/nmo_per_kpt'] = nmo_per_kpt
