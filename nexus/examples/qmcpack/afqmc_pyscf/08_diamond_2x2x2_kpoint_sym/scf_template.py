import pyscf.pbc.gto as gto
from pyscf.pbc import scf, dft
import h5py
import sys

$system

mf = scf.KRHF(cell, kpts=kpts)
mf.chkfile = $chkfile
mf.kernel()

from afqmctools.utils.linalg import get_ortho_ao
hcore = mf.get_hcore()
fock = (hcore + mf.get_veff())
X, nmo_per_kpt = get_ortho_ao(cell,kpts,1e-14)
with h5py.File(mf.chkfile) as fh5:
  fh5['scf/hcore'] = hcore
  fh5['scf/fock'] = fock
  fh5['scf/orthoAORot'] = X
  fh5['scf/nmo_per_kpt'] = nmo_per_kpt
