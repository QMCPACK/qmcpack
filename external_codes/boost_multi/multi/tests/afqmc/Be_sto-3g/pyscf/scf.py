import numpy
from pyscf import scf, fci, gto, ao2mo
from afqmctools.hamiltonian.mol import write_hamil_mol
from afqmctools.wavefunction.mol import write_qmcpack_wfn
from afqmctools.utils.pyscf_utils import load_from_pyscf_chk_mol

mol = gto.M(atom=[('Be', 0, 0, 0)], basis='sto-3g', verbose=0)
mf = scf.RHF(mol)
mf.chkfile = 'scf.chk'
ehf = mf.kernel()
cisolver = fci.direct_spin1.FCI(mol)
hcore = mf.get_hcore()
h1e = numpy.dot(mf.mo_coeff.T, numpy.dot(hcore, mf.mo_coeff))
eri = ao2mo.kernel(mol, mf.mo_coeff, aosym=1)
e_fci, ci_fci = cisolver.kernel(h1e, eri, h1e.shape[1], mol.nelec,
                                ecore=mol.energy_nuc())
coeff, oa, ob = zip(*fci.addons.large_ci(ci_fci, mf.mo_coeff.shape[0],
                    mol.nelec, tol=1e-8, return_strs=False))
scf_data = load_from_pyscf_chk_mol(mf.chkfile)
numpy.random.seed(7)
nmo = h1e.shape[1]
nalpha = mol.nelec[0]
init = numpy.array(numpy.random.rand(nmo,nalpha),dtype=numpy.complex128)
write_hamil_mol(scf_data, 'rqmcpack.h5', 1e-5, real_chol=True)
write_qmcpack_wfn('rqmcpack.h5', (numpy.array(coeff),oa,ob), True,
                  mol.nelec, h1e.shape[1], init=[init,init])
