import ast
import h5py
import numpy

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
    pos = numpy.int64(0)
    occs = []
    shift = 0
    all_found = False
    for d in det:
        while pos < min(nmo,64):
            if d & (1<<pos):
                nset += 1
                occs.append(pos+shift)
            if nset == nel:
                all_found = True
                break
            pos += 1
        # Assuming 64 bit integers
        pos = 0
        if all_found:
            break
        shift += 64
    return occs

def read_qmcpack_ci_wavefunction(input_file, nelec, nmo, ndets=None):
    if ndets is None:
        ndets = -1
    na, nb = nelec
    with h5py.File(input_file) as fh5:
        try:
            nmo = fh5['parameters/numMO'][:][0]
            na = fh5['parameters/NbAlpha'][:][0]
            nb = fh5['parameters/NbBeta'][:][0]
        except KeyError:
            pass
        ci_a = fh5['MultiDet/CI_Alpha'][:]
        ci_b = fh5['MultiDet/CI_Beta'][:]
        nbs = fh5['MultiDet/Nbits'][()]
        coeffs = fh5['MultiDet/Coeff'][:][:ndets]
        coeffs_imag = fh5['MultiDet/Coeff_imag'][:][:ndets]
        coeffs = coeffs + 1j*coeffs_imag
        occa = []
        occb = []
        for i, (ca, cb) in enumerate(zip(ci_a[:ndets], ci_b[:ndets])):
            occa.append(get_occupied(ca, na, nmo))
            occb.append(get_occupied(cb, nb, nmo))
    wfn = (coeffs, numpy.array(occa), numpy.array(occb))
    return wfn, True, nmo, (na,nb)
