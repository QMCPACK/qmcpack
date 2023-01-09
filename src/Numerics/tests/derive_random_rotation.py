
# Generate a 3D rotation matrix used to test the code in RotationMatrix3D.h.
# Output of the script is used in test_RotationMatrix3D.cpp.
# The rotation matrix is used for random rotations of the spherical integration grid
# for the non-local pseudopotential.

from sympy import *
from math import pi,acos

# Make a 2x2 rotation matrix.
#  c - cos of value
#  s - sin of value
def make_rotation(c, s):
    return Matrix([[c, -s],[s, c]])

# Enlarge a 2x2 matrix to 3x3, inserting identity elements as appropriate
def enlarge_matrix(m, idx):
    if idx == 2:
        mb = m.row_insert(2, Matrix([[0, 0]]))
        mc = mb.col_insert(2, Matrix([[0],[0],[1]]))
    if idx == 1:
        mb = m.row_insert(1, Matrix([[0, 0]]))
        mc = mb.col_insert(1, Matrix([[0],[1],[0]]))
    if idx == 0:
        mb = m.row_insert(0, Matrix([[0, 0]]))
        mc = mb.col_insert(0, Matrix([[1],[0],[0]]))
    return mc

theta = Symbol('theta')
phi = Symbol('phi')
psi = Symbol('psi')

cth = cos(theta)
sth = sin(theta)
sph = sin(phi)
cph = cos(phi)
sps = sin(psi)
cps = cos(psi)


def get_rotation_matrix():
    # Based on comments from the code in
    # https://github.com/QMCPACK/qmcpack/commit/a9b5aadbe8b37609c7f65222df8d1daec5f97f93#diff-0d7403fee451dfd57cc1cc6b2509e25d2d7e52bc1877d1c547176e63cc2e3b71R218-R241

    r1 = make_rotation(cph, -sph) # Counterclockwise?
    r2 = make_rotation(cth, sth)  # Clockwise?
    r3 = make_rotation(cps, -sps) # Counterclockwise?

    r1a = enlarge_matrix(r1, 2)  # Rotate around z
    r2a = enlarge_matrix(r2, 1)  # Rotate around y
    r3a = enlarge_matrix(r3, 2)  # Rotate around z
    print('r1',r1a)

    m = r1a * r2a * r3a
    print(m)

    return m


# Taken from QMCWaveFunctions/tests/gen_matrix_ops.py
# TODO - share python code between directories
# Print a check for each matrix entry
def print_matrix_for_check(m, matrix_name):
    for i in range(m.rows):
        for j in range(m.cols):
            print("  CHECK({matrix_name}({row}, {col}) == Approx({val:15g}));".format(
                matrix_name=matrix_name, row=i, col=j, val=m[i,j]))

# To get uniform covering of the sphere,
# psi and phi should sampled from 0 to 2*pi
# and cos(theta) should be sampled from from -1 to 1

rmat = get_rotation_matrix()

s1 = {'phi':0.0, 'p2i':0.0, 'theta':0.0}
print('# ',s1)
m1 = rmat.subs(s1)
print(m1)

rng2 = [0.1, 0.2, 0.3]
s2 = {'psi':2*pi*rng2[0], 'phi':2*pi*rng2[1], 'theta':acos(1 - 2*rng2[2])}

print('\n# rng inputs: ',rng2)
print('# angles: ',s2)
m2 = rmat.subs(s2)
print_matrix_for_check(m2,"rmat2")

rng3 = [0.9, 0.5, 0.8]
s3 = {'psi':2*pi*rng3[0], 'phi':2*pi*rng3[1], 'theta':acos(1 - 2*rng3[2])}

print('\n# rng inputs: ',rng3)
print('#  angles: ',s3)
m3 = rmat.subs(s3)
print_matrix_for_check(m3,"rmat3")
