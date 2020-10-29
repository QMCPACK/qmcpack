import ast
import h5py
import numpy
import scipy.linalg
import sys
from pyscf import fci
from afqmctools.utils.io import to_qmcpack_complex, add_group, add_dataset

def write_wfn_mol(scf_data, ortho_ao, filename, wfn=None,
                  init=None, verbose=False):
    """Generate QMCPACK trial wavefunction.

    Parameters
    ----------
    scf_data : dict
        Dictionary containing scf data extracted from pyscf checkpoint file.
    ortho_ao : bool
        Whether we are working in orthogonalised AO basis or not.
    filename : string
        HDF5 file path to store wavefunction to.
    wfn : tuple
        User defined wavefunction. Not fully supported. Default None.

    Returns
    -------
    wfn : :class:`numpy.ndarray`
        Wavefunction as numpy array. Format depends on wavefunction.
    """
    ghf = False
    mol = scf_data['mol']
    nelec = mol.nelec
    nalpha, nbeta = nelec
    C = scf_data['mo_coeff']
    X = scf_data['X']
    rohf = scf_data.get('rohf', False)
    uhf = scf_data['isUHF']
    # For RHF only nalpha entries will be filled.
    if uhf:
        norb = C[0].shape[-1]
    else:
        norb = C.shape[-1]
    if wfn is None:
        wfn = numpy.zeros((1,norb,nalpha+nbeta), dtype=numpy.complex128)
        wfn_type = 'NOMSD'
        coeffs = numpy.array([1.0+0j])
        if ortho_ao:
            if X.shape[0] != X.shape[1]:
                Xinv = scipy.linalg.pinv(X)
            else:
                Xinv = scipy.linalg.inv(X)
            if uhf:
                # We are assuming C matrix is energy ordered.
                wfn[0,:,:nalpha] = numpy.dot(Xinv, C[0])[:,:nalpha]
                wfn[0,:,nalpha:] = numpy.dot(Xinv, C[1])[:,:nbeta]
            elif rohf:
                wfn[0,:,:nalpha] = numpy.dot(Xinv, C)[:,:nalpha]
                wfn[0,:,nalpha:] = numpy.dot(Xinv, C)[:,:nbeta]
            else:
                wfn[0,:,:nalpha] = numpy.dot(Xinv, C)[:,:nalpha]
        else:
            # Assuming we are working in MO basis, only works for RHF, ROHF trials.
            I = numpy.identity(C.shape[-1], dtype=numpy.float64)
            wfn[0,:,:nalpha] = I[:,:nalpha]
            if rohf:
                wfn[0,:,nalpha:] = I[:,:nbeta]
            if uhf:
                print(" # Warning: UHF trial wavefunction can only be used of "
                      "working in ortho AO basis.")
    write_qmcpack_wfn(filename, (numpy.array([1.0+0j]),wfn), uhf or rohf,
                      nelec, norb, verbose=verbose)
    return nelec

#def write_qmcpack_wfn(filename, wfn, walker_type, nelec, norb, init=None,
def write_qmcpack_wfn(filename, wfn, uhf, nelec, norb, init=None,
                      orbmat=None, verbose=False):
    # User defined wavefunction.
    # PHMSD is a list of tuple of (ci, occa, occb).
    # NOMSD is a tuple of (list, numpy.ndarray).
    if len(wfn) == 3:
        coeffs, occa, occb = wfn
        wfn_type = 'PHMSD'
    elif len(wfn) == 2:
        coeffs, wfn = wfn
        wfn_type = 'NOMSD'
    else:
        print("Unknown wavefunction type passed.")
        sys.exit()

    fh5 = h5py.File(filename, 'a')
    nalpha, nbeta = nelec
    # TODO: FIX for GHF eventually.
#    if walker_type == 'ghf':
#        walker_type = 3
#    elif walker_type == 'uhf':
#        walker_type = 2
#        uhf = True
#    else:
#        walker_type = 1
#        uhf = False
    if uhf:
        walker_type = 2
    else:
        walker_type = 1
    if wfn_type == 'PHMSD':
        walker_type = 2
    if wfn_type == 'NOMSD':
        wfn_group = add_group(fh5, 'Wavefunction/NOMSD')
        write_nomsd(wfn_group, wfn, uhf, nelec, init=init)
    else:
        wfn_group = add_group(fh5, 'Wavefunction/PHMSD')
        write_phmsd(wfn_group, occa, occb, nelec, norb,
                    init=init, orbmat=orbmat)
    if coeffs.dtype == float:
        if verbose:
            print(" # Found real MSD coefficients. Converting to complex.")
        coeffs = numpy.array(coeffs, dtype=numpy.complex128)
    wfn_group['ci_coeffs'] = to_qmcpack_complex(coeffs)
    dims = [norb, nalpha, nbeta, walker_type, len(coeffs)]
    wfn_group['dims'] = numpy.array(dims, dtype=numpy.int32)
    fh5.close()

def write_nomsd(fh5, wfn, uhf, nelec, thresh=1e-8, init=None):
    """Write NOMSD to HDF.

    Parameters
    ----------
    fh5 : h5py group
        Wavefunction group to write to file.
    wfn : :class:`numpy.ndarray`
        NOMSD trial wavefunctions.
    uhf : bool
        UHF style wavefunction.
    nelec : tuple
        Number of alpha and beta electrons.
    thresh : float
        Threshold for writing wavefunction elements.
    """
    nalpha, nbeta = nelec
    nmo = wfn.shape[1]
    wfn[abs(wfn) < thresh] = 0.0
    if init is not None:
        add_dataset(fh5, 'Psi0_alpha', to_qmcpack_complex(init[0]))
        add_dataset(fh5, 'Psi0_beta', to_qmcpack_complex(init[1]))
    else:
        add_dataset(fh5, 'Psi0_alpha',
                    to_qmcpack_complex(wfn[0,:,:nalpha].copy()))
        if uhf:
            add_dataset(fh5, 'Psi0_beta',
                        to_qmcpack_complex(wfn[0,:,nalpha:].copy()))
    for idet, w in enumerate(wfn):
        # QMCPACK stores this internally as a csr matrix, so first convert.
        ix = 2*idet if uhf else idet
        psia = scipy.sparse.csr_matrix(w[:,:nalpha].conj().T)
        write_nomsd_single(fh5, psia, ix)
        if uhf:
            ix = 2*idet + 1
            psib = scipy.sparse.csr_matrix(w[:,nalpha:].conj().T)
            write_nomsd_single(fh5, psib, ix)

def write_nomsd_single(fh5, psi, idet):
    """Write single component of NOMSD to hdf.

    Parameters
    ----------
    fh5 : h5py group
        Wavefunction group to write to file.
    psi : :class:`scipy.sparse.csr_matrix`
        Sparse representation of trial wavefunction.
    idet : int
        Determinant number.
    """
    base = 'PsiT_{:d}/'.format(idet)
    dims = [psi.shape[0], psi.shape[1], psi.nnz]
    fh5[base+'dims'] = numpy.array(dims, dtype=numpy.int32)
    fh5[base+'data_'] = to_qmcpack_complex(psi.data)
    fh5[base+'jdata_'] = psi.indices
    fh5[base+'pointers_begin_'] = psi.indptr[:-1]
    fh5[base+'pointers_end_'] = psi.indptr[1:]

def write_phmsd(fh5, occa, occb, nelec, norb, init=None, orbmat=None):
    """Write NOMSD to HDF.

    Parameters
    ----------
    fh5 : h5py group
        Wavefunction group to write to file.
    nelec : tuple
        Number of alpha and beta electrons.
    """
    # TODO: Update if we ever wanted "mixed" phmsd type wavefunctions.
    na, nb = nelec
    if init is not None:
        add_dataset(fh5, 'Psi0_alpha', to_qmcpack_complex(init[0]))
        add_dataset(fh5, 'Psi0_beta', to_qmcpack_complex(init[1]))
    else:
        init = numpy.eye(norb, dtype=numpy.complex128)
        add_dataset(fh5, 'Psi0_alpha',
                    to_qmcpack_complex(init[:,occa[0]].copy()))
        add_dataset(fh5, 'Psi0_beta',
                    to_qmcpack_complex(init[:,occb[0]].copy()))
    if orbmat is not None:
        fh5['type'] = 1
        # Expects conjugate transpose.
        oa = scipy.sparse.csr_matrix(orbmat[0].conj().T)
        write_nomsd_single(fh5, oa, 0)
        ob = scipy.sparse.csr_matrix(orbmat[1].conj().T)
        write_nomsd_single(fh5, ob, 1)
    else:
        fh5['type'] = 0
    occs = numpy.zeros((len(occa), na+nb), dtype=numpy.int32)
    occs[:,:na] = numpy.array(occa)
    occs[:,na:] = norb+numpy.array(occb)
    # Reading 1D array currently in qmcpack.
    fh5['occs'] = occs.ravel()

#
# Graveyard. Old QMCPACK wavefunction plain text format.
# Keep around for backwards compatability.
#
def write_nomsd_wfn(filename, wfn, nalpha, uhf, coeffs=[1.0]):
    if len(wfn.shape) == 2:
        wfn = wfn.reshape((1,wfn.shape[0],wfn.shape[1]))
    namelist = qmcpack_wfn_namelist(wfn.shape[0], uhf)
    with open(filename, 'a') as f:
        f.write(namelist)
        f.write('Coefficients: ' + ' '.join(str(c) for c in coeffs) +'\n')
        for (i,d) in enumerate(wfn):
            f.write('Determinant: {}\n'.format(i+1))
            if uhf:
                write_single(f, d[:,:nalpha])
                write_single(f, d[:,nalpha:])
            else:
                write_single(f, d[:,:nalpha])


def qmcpack_wfn_namelist(nci, uhf):
    return ("&FCI\n UHF = {}\n CMajor\n "
            "NCI = {}\n TYPE = matrix\n/\n".format(int(uhf),nci))

def write_single(out, mos):
    for j in range(0, mos.shape[1]):
        for i in range(0, mos.shape[0]):
            val = mos[i,j]
            out.write('(%.10e,%.10e) '%(val.real, val.imag))
        out.write('\n')

def gen_multi_det_wavefunction(mc, weight_cutoff=0.95, verbose=False,
                               max_ndets=1e5, norb=None,
                               filename=None):
    """Generate multi determinant particle-hole trial wavefunction.

    Format adopted to be compatable with QMCPACK PHMSD type wavefunction.

    Parameters
    ----------
    mc : pyscf CI solver type object
        Input object containing multi determinant coefficients.
    weight_cutoff : float, optional
        Print determinants until accumulated weight equals weight_cutoff.
        Default 0.95.
    verbose : bool
        Print information about process. Default False.
    max_ndets : int
        Max number of determinants to print out. Default 1e5.
    norb : int or None, optional
        Total number of orbitals in simulation. Used if we want to run CI within
        active space but QMC in full space. Deault None.
    filename : string
        Output filename. Default "multi_det.dat"
    """
    occlists = fci.cistring._gen_occslst(range(mc.ncas), mc.nelecas[0])

    ci_coeffs = mc.ci.ravel()
    # Sort coefficients in terms of increasing absolute weight.
    ix_sort = numpy.argsort(numpy.abs(ci_coeffs))[::-1]
    cweight = numpy.cumsum(ci_coeffs[ix_sort]**2)
    max_det = numpy.searchsorted(cweight, weight_cutoff)
    ci_coeffs = ci_coeffs[ix_sort]
    ndets = min(max_det,max_ndets)
    if verbose:
        print(" # Number of dets in CI expansion: {:d}".format(ndets))

    output = open(filename, 'w')
    namelist = "&FCI\n UHF = 0\n NCI = %d\n TYPE = occ\n&END" % ndets
    output.write(namelist+'\n')
    output.write("Configurations:"+'\n')
    if norb is None:
        norb = mc.ncas

    occups = []
    occdns = []
    coeffs = []
    for idet in range(min(max_det,max_ndets)):
        if mc.ncore > 0:
            ocore_up = ' '.join('{:d}'.format(x+1) for x in range(mc.ncore))
            ocore_dn = ' '.join('{:d}'.format(x+1+norb) for x in range(mc.ncore))
        else:
            ocore_up = ' '
            ocore_dn = ' '
        coeff = '%.13f'%ci_coeffs[idet]
        coeffs.append(ci_coeffs[idet])
        ix_alpha = ix_sort[idet] // len(occlists)
        ix_beta = ix_sort[idet] % len(occlists)
        ia = occlists[ix_alpha]
        ib = occlists[ix_beta]
        oup = ' '.join('{:d}'.format(x+1+mc.ncore) for x in ia)
        odown = ' '.join('{:d}'.format(x+norb+1+mc.ncore) for x in ib)
        occups.append([int(o) for o in oup.split()])
        occdns.append([int(o) for o in odown.split()])
        output.write(coeff+' '+ocore_up+' '+oup+' '+ocore_dn+' '+odown+'\n')
    return coeffs, [occups,occdns]


def write_phmsd_wfn(filename, occs, nmo, ncore=0):
    output = open(filename, 'w')
    namelist = "&FCI\n UHF = 0\n NCI = %d\n TYPE = occ\n&END" % len(occs)
    output.write(namelist+'\n')
    output.write("Configurations:"+'\n')
    corea = [i + 1 for i in range(ncore)]
    coreb = [i + nmo + 1 for i in range(ncore)]
    for c, occup, occdn in occs:
        # occup = corea + [ncore + oa + 1 for oa in da.tolist()]
        # occdn = coreb + [ncore + nmo + ob + 1 for ob in db.tolist()]
        # print(occup, occdn)
        occstra = ' '.join('{:d} '.format(x+1) for x in occup)
        occstrb = ' '.join('{:d}'.format(x+1) for x in occdn)
        output.write('%13.8e '%c + occstra + occstrb + '\n')
