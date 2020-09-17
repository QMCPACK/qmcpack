import ast
import h5py
import numpy
import scipy.sparse
import struct
from afqmctools.utils.io import from_qmcpack_complex

def read_qmcpack_ascii_wavefunction(filename, nmo, nelec):
    na, nb = nelec
    cmajor = False
    nci = 1
    fullmo = False
    with open(filename) as f:
        line = f.readline()
        cnt = 1
        while line:
            line = f.readline().strip().split('=')
            if len(line) == 2:
                if line[0].strip() == 'NCI':
                    nci = int(line[1])
                elif line[0].strip() == 'TYPE':
                    wfn_type = line[1].strip()
                elif line[0].strip() == 'UHF':
                    uhf = int(line[1])
            else:
                if line[0] == 'CMajor':
                    cmajor = True
                elif line[0] == 'FullMO':
                    fullmo = True
                elif line[0] == '&END' or line[0] == '/':
                    break
                elif len(line) == 0:
                    break
        if wfn_type == 'occ':
            uhf = True
            wfn = read_phmsd(f, na, nb, nmo)
        else:
            wfn = read_nomsd(f, nmo, na, nb, nci, uhf, fullmo, cmajor)

    return wfn, True if uhf else False

def read_phmsd(f, na, nb, nmo):
    line = f.readline()
    coeffs = []
    occa = []
    occb = []
    while line:
        line = f.readline().split()
        if len(line) > 0:
            coeffs.append(convert_string(line[0]))
            occa.append([int(i)-1 for i in line[1:na+1]])
            occb.append([int(i)-1-nmo  for i in line[na+1:]])

    return numpy.array(coeffs), numpy.array(occa), numpy.array(occb)

def read_nomsd(f, nmo, na, nb, nci, uhf, fullmo, cmajor):
    wfn = numpy.zeros((nci,nmo,na+nb), dtype=numpy.complex128)
    cnt = 0
    coeffs = []
    while True:
        line = f.readline().split()
        if len(line) == 1 and line[0] == 'Coefficients:':
            continue
        elif len(line) > 1 and line[0] == 'Coefficients:':
            coeffs = [convert_string(s) for s in line[1:]]
        elif 'Determinant' in line[0]:
            break
        else:
            coeffs += [convert_string(s) for s in line]
    assert nci == len(coeffs)
    data = []
    while True:
        line = f.readline().split()
        if len(line) > 0:
            if 'Determinant' in line[0]:
                continue
            for v in line:
                val = convert_string(v)
                data.append(val)
        else:
            break
    if fullmo:
        if uhf:
            nvals = 2*nmo*nmo
        else:
            nvals = nmo*nmo
        noa = nmo*nmo
        nob = nmo*nmo
        shapea = (nmo,nmo)
        shapeb = (nmo,nmo)
    else:
        if uhf:
            nvals = nmo*(na+nb)
            nvals = nmo*(na+nb)
        else:
            nvals = nmo*na
        noa = nmo*na
        nob = nmo*nb
        shapea = (nmo,na)
        shapeb = (nmo,nb)
    assert len(data) == nvals*nci
    nspin = 2 if uhf else 1
    order = 'F' if cmajor else 'C'
    for i in range(nci):
        orbs = data[nspin*i*noa:(nspin*i+1)*noa]
        wfn[i,:,:na] = numpy.array(orbs).reshape(shapea, order=order)[:,:na]
        if uhf:
            orbs = data[(nspin*i+1)*noa:(nspin*i+2)*noa]
            wfn[i,:,na:] = numpy.array(orbs).reshape(shapeb, order=order)[:,:nb]

    return (numpy.array(coeffs), wfn)

def convert_string(s):
    try:
        c = complex(s)
    except ValueError:
        c = ast.literal_eval(s)
        c = c[0] + 1j*c[1]
    return c

def read_orbitals():
    with open(filename) as f:
        content = f.readlines()[nskip:]
    useable = numpy.array([c.split() for c in content]).flatten()
    tuples = [ast.literal_eval(u) for u in useable]
    orbs = [complex(t[0], t[1]) for t in tuples]
    return numpy.array(orbs)

def get_occupied(det, nel, nmo):
    nset = 0
    pos = numpy.uint64(0)
    occs = []
    shift = 0
    one = numpy.uint64(1)
    all_found = False
    for d in det:
        while pos < min(nmo,64):
            if d & (one<<pos):
                nset += 1
                occs.append(int(pos+shift))
            if nset == nel:
                all_found = True
                break
            pos += numpy.uint64(1)
        # Assuming 64 bit integers
        pos = 0
        if all_found:
            break
        shift += 64
    return occs

def read_dmc_ci_wavefunction(input_file, nelec, nmo, ndets=None):
    if ndets is None:
        ndets = -1
    with h5py.File(input_file) as fh5:
        try:
            nmo_read = fh5['parameters/numMO'][:][0]
            na = fh5['parameters/NbAlpha'][:][0]
            nb = fh5['parameters/NbBeta'][:][0]
            if nelec is not None:
                assert na == nelec[0]
                assert nb == nelec[1]
            if nmo is not None:
                assert nmo_read == nmo
            else:
                nmo = nmo_read
            nelec = (na, nb)
        except KeyError:
            pass
        assert nelec is not None
        assert nmo is not None
        ci_a = numpy.array(fh5['MultiDet/CI_Alpha'][:], dtype=numpy.uint64)
        ci_b = numpy.array(fh5['MultiDet/CI_Beta'][:], dtype=numpy.uint64)
        nbs = fh5['MultiDet/Nbits'][()]
        coeffs = fh5['MultiDet/Coeff'][:][:ndets]
        try:
            coeffs_imag = fh5['MultiDet/Coeff_imag'][:][:ndets]
        except KeyError:
            coeffs_imag = 0j
        coeffs = coeffs + 1j*coeffs_imag
        occa = []
        occb = []
        for i, (ca, cb) in enumerate(zip(ci_a[:ndets], ci_b[:ndets])):
            occa.append(get_occupied(ca, na, nmo))
            occb.append(get_occupied(cb, nb, nmo))
    wfn = (coeffs, numpy.array(occa), numpy.array(occb))
    return wfn, True, nmo, (na,nb)

def write_phf_rhf(fname, cf):
    assert cf.dtype == numpy.dtype(float) or cf.dtype == numpy.dtype(complex)
    assert len(cf.shape) == 2
    nb, no = cf.shape

    ia = numpy.array([1,nb,no,1], dtype=int)
    i1_ = struct.pack('i',4*4)
    i2_ = struct.pack('i',16*nb*no)
    ia_ = struct.pack('i'*4,*ia)

    cfx = numpy.ascontiguousarray(numpy.zeros((2*nb*no,),
                                              dtype=float,
                                              order='C'))
    if numpy.iscomplexobj(cf):
        cfx[0::2] = cf.real.reshape((nb*no,), order='F')[:]
        cfx[1::2] = cf.imag.reshape((nb*no,), order='F')[:]
    else:
        cfx[0::2] = cf.reshape((nb*no,), order='F')[:]
    cf_ = struct.pack('d'*(nb*no*2), *cfx)

    with open(fname, 'wb') as f:
        f.write(i1_)
        f.write(ia_)
        f.write(i1_)
        f.write(i2_)
        f.write(cf_)
        f.write(i2_)

def write_phf_uhf(fname, cf):
    assert len(cf.shape) == 3
    assert cf[0].dtype == numpy.dtype(float) or cf[0].dtype == numpy.dtype(complex)
    assert cf[1].dtype == numpy.dtype(float) or cf[1].dtype == numpy.dtype(complex)
    assert cf[0].shape == cf[1].shape
    nb, no = cf[0].shape

    ia = numpy.array([2,nb,no,1], dtype=int)
    i1_ = struct.pack('i',4*4)
    i2_ = struct.pack('i',32*nb*no)
    ia_ = struct.pack('i'*4,*ia)

    cfx = numpy.ascontiguousarray(numpy.zeros((4*nb*no,),
                                              dtype=float,
                                              order='C'))
    if numpy.iscomplexobj(cf[0]):
        cfx[0:2*nb*no:2] = cf[0].real.reshape((nb*no,), order='F')[:]
        cfx[1:2*nb*no:2] = cf[0].imag.reshape((nb*no,), order='F')[:]
    else:
        cfx[0:2*nb*no:2] = cf[0].reshape((nb*no,), order='F')[:]
    if numpy.iscomplexobj(cf[1]):
        cfx[2*nb*no+0::2] = cf[1].real.reshape((nb*no,), order='F')[:]
        cfx[2*nb*no+1::2] = cf[1].imag.reshape((nb*no,), order='F')[:]
    else:
        cfx[2*nb*no+0::2] = cf[1].reshape((nb*no,), order='F')[:]
    cf_ = struct.pack('d'*(nb*no*4), *cfx)

    with open(fname, 'wb') as f:
        f.write(i1_)
        f.write(ia_)
        f.write(i1_)
        f.write(i2_)
        f.write(cf_)
        f.write(i2_)

def read_qmcpack_wavefunction(filename):
    try:
        with h5py.File(filename, 'r') as fh5:
            wgroup = fh5['Wavefunction/NOMSD']
            wfn, psi0, nelec = read_qmcpack_nomsd_hdf5(wgroup)
    except KeyError:
        with h5py.File(filename, 'r') as fh5:
            wgroup = fh5['Wavefunction/PHMSD']
            wfn, psi0, nelec = read_qmcpack_phmsd_hdf5(wgroup)
    except KeyError:
        print("Wavefunction not found.")
        sys.exit()
    return wfn, psi0, nelec

def read_qmcpack_nomsd_hdf5(wgroup):
    dims = wgroup['dims']
    nmo = dims[0]
    na = dims[1]
    nb = dims[2]
    walker_type = dims[3]
    if walker_type == 2:
        uhf = True
    else:
        uhf = False
    nci = dims[4]
    coeffs = from_qmcpack_complex(wgroup['ci_coeffs'][:], (nci,))
    psi0a = from_qmcpack_complex(wgroup['Psi0_alpha'][:], (nmo,na))
    if uhf:
        psi0b = from_qmcpack_complex(wgroup['Psi0_beta'][:], (nmo,nb))
    psi0 = numpy.zeros((nmo,na+nb),dtype=numpy.complex128)
    psi0[:,:na] = psi0a.copy()
    if uhf:
        psi0[:,na:] = psi0b.copy()
    else:
        psi0[:,na:] = psi0a[:,:nb].copy()
    wfn = numpy.zeros((nci,nmo,na+nb), dtype=numpy.complex128)
    for idet in range(nci):
        ix = 2*idet if uhf else idet
        pa = orbs_from_dset(wgroup['PsiT_{:d}/'.format(idet)])
        wfn[idet,:,:na] = pa
        if uhf:
            ix = 2*idet + 1
            wfn[idet,:,na:] = orbs_from_dset(wgroup['PsiT_{:d}/'.format(ix)])
        else:
            wfn[idet,:,na:] = pa[:,:nb]
    return (coeffs,wfn), psi0, (na, nb)

def read_qmcpack_phmsd_hdf5(wgroup):
    dims = wgroup['dims']
    nmo = dims[0]
    na = dims[1]
    nb = dims[2]
    walker_type = dims[3]
    if walker_type == 2:
        uhf = True
    else:
        uhf = False
    nci = dims[4]
    coeffs = from_qmcpack_complex(wgroup['ci_coeffs'][:], (nci,))
    occs = wgroup['occs'][:].reshape((nci,na+nb))
    occa = occs[:,:na]
    occb = occs[:,na:]-nmo
    wfn = (coeffs, occa, occb)
    psi0a = from_qmcpack_complex(wgroup['Psi0_alpha'][:], (nmo,na))
    if uhf:
        psi0b = from_qmcpack_complex(wgroup['Psi0_beta'][:], (nmo,nb))
    psi0 = numpy.zeros((nmo,na+nb),dtype=numpy.complex128)
    psi0[:,:na] = psi0a.copy()
    if uhf:
        psi0[:,na:] = psi0b.copy()
    else:
        psi0[:,na:] = psi0a.copy()
    return wfn, psi0, (na,nb)

def orbs_from_dset(dset):
    """Will read actually A^{H} but return A.
    """
    dims = dset['dims'][:]
    wfn_shape = (dims[0],dims[1])
    nnz = dims[2]
    data = from_qmcpack_complex(dset['data_'][:],(nnz,))
    indices = dset['jdata_'][:]
    pbb = dset['pointers_begin_'][:]
    pbe = dset['pointers_end_'][:]
    indptr = numpy.zeros(dims[0]+1)
    indptr[:-1] = pbb
    indptr[-1] = pbe[-1]
    wfn = scipy.sparse.csr_matrix((data,indices,indptr),shape=wfn_shape)
    return wfn.toarray().conj().T.copy()
