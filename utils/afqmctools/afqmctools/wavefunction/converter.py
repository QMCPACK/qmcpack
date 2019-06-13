import numpy
import ast

def read_qmcpack_ascii_wavefunction(filename, nmo, nelec):
    na, nb = nelec
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
            wfn = read_phmsd(f, na, nb)
        else:
            wfn = read_nomsd(f, nmo, na, nb, nci, uhf, fullmo)

    return wfn, 'uhf' if uhf else 'rhf'

def read_phmsd(f, na, nb):
    line = f.readline()
    coeffs = []
    occa = []
    occb = []
    while line:
        line = f.readline().split()
        if len(line) > 0:
            coeffs.append(convert_string(line[0]))
            occa.append([int(i) for i in line[1:na+1]])
            occb.append([int(i) for i in line[na+1:]])

    return numpy.array(coeffs), numpy.array(occa), numpy.array(occb)

def read_nomsd(f, nmo, na, nb, nci, uhf, fullmo):
    wfn = numpy.zeros((nci,nmo,na+nb), dtype=numpy.complex128)
    line = f.readline().split()
    coeffs = [convert_string(s) for s in line[1:]]
    line = f.readline().split()
    data = []
    while True:
        line = f.readline().split()
        if len(line) > 0:
            for v in line:
                try:
                    val = convert_string(v)
                except SyntaxError:
                    pass
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
    for i in range(nci):
        orbs = data[i*noa:(i+1)*noa]
        wfn[i,:,:na] = numpy.array(orbs).reshape(shapea)[:,:na]
        if uhf:
            orbs = data[(i+1)*noa:(i+2)*noa]
            wfn[i,:,na:] = numpy.array(orbs).reshape(shapeb)[:,:nb]

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

def read_phfmol_wfn(filename, nmo):
    with open(filename) as f:
        content = f.read().split()
    start = False
    idet = 0
    data = []
    for (i,f) in enumerate(content):
        if 'NCI' in f:
            try:
                ndets = int(content[i+1])
            except ValueError:
                ndets = int(content[i+2])
            dets = numpy.zeros((ndets,nmo,nmo), dtype=numpy.complex128)
        # print(f,start,data)
        # print(len(data),f)
        if 'Coefficients' in f:
            string_coeffs = content[i+1:i+1+ndets]
        if 'Determinant' in f:
            break
    start = i + 2
    coeffs = []
    for c in string_coeffs:
        v = ast.literal_eval(c)
        coeffs.append(complex(v[0],v[1]))

    for idet in range(ndets):
        end = start+nmo*nmo
        data = []
        for line in content[start:end]:
            v = ast.literal_eval(line)
            data.append(complex(v[0],v[1]))
        C = numpy.copy(numpy.array(data).reshape(nmo,nmo).T)
        dets[idet] = C
        dets[idet] = C
        start = end + 2
    return numpy.array(coeffs), dets
