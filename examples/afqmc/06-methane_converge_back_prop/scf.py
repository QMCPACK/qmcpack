#! /usr/bin/env python3

from math import cos, sin, pi, acos
import numpy
from pyscf import gto, scf

# Bond angle HCH
alpha = acos(-1.0/3.0)
# Rotation angle about z.
beta = 2*pi / 3.0
c = cos(beta)
s = sin(beta)
# Rotation matrix about z axis.
R = numpy.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
# C-H bond length in angstrom
rch = 1.1085
# Rotation matrix for geometry
C = numpy.array([0, 0, 0])
H1 = numpy.array([0, 0, rch])
# Arbitrary position of first HCH molecule
H2 = numpy.array([rch*cos(alpha-pi/2.), 0, -rch*sin(alpha-pi/2.)])
H3 = R.dot(H2)
H4 = R.dot(H3)

mol = gto.M(atom=[['C', C], ['H', H1], ['H', H2], ['H', H3], ['H', H4]],
            basis='sto-3g',
            unit='Angstrom')

mf = scf.RHF(mol)
mf.chkfile = 'scf.chk'
ehf = mf.scf()
