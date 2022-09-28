# Use sympy to verify the matrix exponential for some concrete numerical cases
# Output used in test_RotatedSPOs.cpp

from sympy import *

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


# Print a check for each matrix entry
def print_matrix_for_check(m, matrix_name):
    for i in range(m.rows):
        for j in range(m.cols):
            print("  CHECK({matrix_name}({row},{col}) == ValueApprox({val:15g}));".format(
                matrix_name=matrix_name, row=i, col=j, val=m[i,j]))


# Print matrix as initializer list (for use in checkMatrix)
def print_matrix_as_intializer_list(m):
    print("{", end="")
    comma = ","
    for i in range(m.rows):
        for j in range(m.cols):
            # Skip comma after the last entry
            if i == m.rows-1 and j == m.cols-1:
                comma = ""
            print(" {:>18.15g}{comma}".format(m[i,j], comma=comma), end="")
        if i == m.rows-1:
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
print_matrix_as_intializer_list(m2)
print("\nExp(Input)")
print_matrix_as_intializer_list(m2exp)
print()

print("3x3 matrix")
m3 = create3x3_matrix(0.3, 0.1, 0.2)
m3exp = exp(m3)

print("Input")
print_matrix_as_intializer_list(m3)
print("\nExp(Input)")
print_matrix_as_intializer_list(m3exp)
