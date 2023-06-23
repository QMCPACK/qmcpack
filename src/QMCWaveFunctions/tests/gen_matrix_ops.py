# Use sympy to verify the matrix exponential for some concrete numerical cases
# Output used in test_RotatedSPOs.cpp

from sympy import *
import scipy.linalg
import numpy as np

# Create 2x2 skew symmetric matrix
def create2x2_matrix(k1):
    return Matrix([[0.0, -k1],
                   [k1, 0.0]])

# Create 3x3 skew symmetric matrix
# The pattern of plus/minus reflects how the matrix is constructed in QMCPACK.
# The standard presentation for a 3x3 matrix (that creates a rotation matrix)
# numbers the variables differently, and flips the sign on k2
def create3x3_matrix(k1, k2, k3):
    return Matrix([[0.0, -k1, -k2],
                   [k1,  0.0, -k3],
                   [k2,   k3, 0.0]])

def extract3x3parameters(m3):
    k1 = m3[1,0]
    k2 = m3[2,0]
    k3 = m3[1,2]

    # Test that it is really an antisymmetric matrix
    tol = 1e-12
    for i in range(3):
        assert(abs(m3[i,i]) < tol)

    for (i,j) in [(1,0), (2,0), (1,2)]:
        assert(abs(m3[i,j] + m3[j,i]) < tol)

    return k1, k2, k3

def create4x4_matrix(k1, k2, k3, k4, k5, k6):
    return np.array([[0.0, -k1, -k2, -k3],
                     [k1,  0.0, -k4, -k5],
                     [k2,   k4, 0.0, -k6],
                     [k3,   k5,  k6, 0.0]])

def extract4x4parameters(m4):
    k1 = m4[1,0]
    k2 = m4[2,0]
    k3 = m4[3,0]
    k4 = m4[2,1]
    k5 = m4[3,1]
    k6 = m4[3,2]

    # Test that it is really an antisymmetric matrix
    tol = 1e-12
    for i in range(4):
        assert(abs(m4[i,i]) < tol)

    for (i,j) in [(1,0), (2,0), (3,0), (1,2), (1,3), (2,3)]:
        assert(abs(m4[i,j] + m4[j,i]) < tol)

    return k1, k2, k3, k4, k5, k6

# Print a check for each matrix entry
def print_matrix_for_check(m, matrix_name):
    for i in range(m.rows):
        for j in range(m.cols):
            print("  CHECK({matrix_name}({row},{col}) == ValueApprox({val:15g}));".format(
                matrix_name=matrix_name, row=i, col=j, val=m[i,j]))


# Print matrix as initializer list (for use in checkMatrix)
# The entries are printed nicely to appear as a matrix, but form a single list in row-major ordering, which is the same as used by OhmmsMatrix.
def print_matrix_as_initializer_list(m):
    try:
        # Sympy matrices
        rows = m.rows
        cols = m.cols
    except AttributeError:
        # Numpy arrays
        rows = m.shape[0]
        cols = m.shape[1]

    print("{", end="")
    comma = ","
    for i in range(rows):
        for j in range(cols):
            # Skip comma after the last entry
            if i == rows-1 and j == cols-1:
                comma = ""
            print(" {:>18.15g}{comma}".format(float(m[i,j]), comma=comma), end="")
        if i == rows-1:
            print(" };")
        else:
            print()
            print(" ", end="")


# Only have one choice for antisymmetric 1x1 matrix
print("1x1 matrix")
m1 = Matrix([0.0])
m1exp = exp(-m1)

print(m1exp)
print()

print("2x2 matrix")
m2 = create2x2_matrix(0.1)
m2exp = exp(m2)

print("Input")
print_matrix_as_initializer_list(m2)
print("\nExp(Input)")
print_matrix_as_initializer_list(m2exp)
print()

print("3x3 matrix")
m3 = create3x3_matrix(0.3, 0.1, 0.2)
m3exp = exp(m3)

print("Input")
print_matrix_as_initializer_list(m3)
print("\nExp(Input)")
print_matrix_as_initializer_list(m3exp)


# Application of rotation matrix in the test "RotatedSPOs construct delta matrix" in test_RotatedSPOs.cpp
# Sympy operations are very slow for a 4x4 matrix, and it doesn't have a matrix log.
# Use the Scipy versions instead.
print()
print("4x4 rotation matrix in 'RotatedSPOs construct delta matrix' test")
print()
# Be aware of the order of the indices when comparing with qmcpack.
# This system has 2 electrons and 4 orbitals.
# The 0,1 (core->core) and 2,3 (unoccupied->unoccupied) are parameters k1 and k6, respectively,
#  and are put at the end of the full rotation indices list.
m4old = create4x4_matrix(-1.1, 1.5, 0.2, -0.15, 0.03, 0.05)
m4delta = create4x4_matrix(0.0, 0.1, 0.3, 0.2, -0.1, 0.0)

# Why is the order reversed compared with constructDeltaRotations?
#  Probably has to do with the data ordering of the matrices?
m4new = np.dot(scipy.linalg.expm(m4old), scipy.linalg.expm(m4delta))
print("New rotation matrix")
print_matrix_as_initializer_list(m4new)

param4 = scipy.linalg.logm(m4new)
print("Corresponding antisymmetric matrix")
print(param4)
ks = extract4x4parameters(param4)
print('Extracted parameters (Reminder: ordering differs from QMCPACK)')
print(ks)


# Application of rotation matrix in test_RotatedSPOs_LCAO.cpp "Rotated LCAO global rotation consistency"
print()
print("3x3 rotation matrix in 'Rotated LCAO global rotation consistency' test")
print()

# For this test, nel = 1 nmo = 3
# The core-unoccupied rotations come first (0,1),(0,2) followed by the unoccupied-unoccupied rotations (1,2)
# In this case they fall in the natural order for k1, k2, and k3.
m3old = np.array(create3x3_matrix(0.1, 0.2, 0.0))
m3next = np.array(create3x3_matrix(0.3, 0.15, 0.0))

print('Original rotation matrix')
print_matrix_as_initializer_list(scipy.linalg.expm(m3old).T)

m3new = np.dot(scipy.linalg.expm(m3next).T, scipy.linalg.expm(m3old).T)
print('After second rotation')
print_matrix_as_initializer_list(m3new)

param3 = scipy.linalg.logm(m3new)

ks3 = extract3x3parameters(param3.T)
print('Extracted parameters after second rotation')
print(ks3)
