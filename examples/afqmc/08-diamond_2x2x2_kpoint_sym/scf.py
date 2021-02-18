#!/usr/bin/env python3

import numpy
import pyscf.pbc.gto as gto
from pyscf.pbc import scf, dft
import h5py
import sys

cell = gto.Cell()
cell.verbose = 5
alat0 = 3.6
cell.a = (numpy.ones((3,3))-numpy.eye(3))*alat0 / 2.0
cell.atom = (('C',0,0,0),('C',numpy.array([0.25,0.25,0.25])*alat0))
cell.basis = 'gth-szv'
cell.pseudo = 'gth-pade'
cell.mesh = [25,25,25]
cell.build()
nk = [2,2,2]
kpts = cell.make_kpts(nk)

mf = scf.KRHF(cell, kpts=kpts)
mf.chkfile = 'scf.chk'
mf.kernel()

from afqmctools.utils.linalg import get_ortho_ao
hcore = mf.get_hcore()
fock = (hcore + mf.get_veff())
X, nmo_per_kpt = get_ortho_ao(cell,kpts,1e-14)
with h5py.File(mf.chkfile, 'r+') as fh5:
  fh5['scf/hcore'] = hcore
  fh5['scf/fock'] = fock
  fh5['scf/orthoAORot'] = X
  fh5['scf/nmo_per_kpt'] = nmo_per_kpt
