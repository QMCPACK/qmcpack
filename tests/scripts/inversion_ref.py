# no shebang on purpose run with $: python3 inversion_ref.py

# //////////////////////////////////////////////////////////////////////////////////////
# // This file is distributed under the University of Illinois/NCSA Open Source License.
# // See LICENSE file in top directory for details.
# //
# // Copyright (c) 2021 QMCPACK developers.
# //
# // File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
# //
# // File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
# //////////////////////////////////////////////////////////////////////////////////////

# this python code was used to generate the reference data used in
# test_cuBLAS_LU.cpp

import numpy as np
import scipy as sp
import scipy.linalg as sl
import cmath
import math
import itertools as it

def complex_det_log(lu_diag, pivots):
    log_value = complex(0.0, 0.0)
    for i in range(4):
        log_value += cmath.log(lu_diag[i] * (1 - (2 * (pivot[i] - 1 != i))))

def complex_cplusplus_format(a_mat):
    for i in range (a_mat.shape[0]):
        for j in range (a_mat.shape[1]):
            print("a({},{}) = {{ {}, {} }};".format(i,j,a_mat[i,j].real,a_mat[i,j].imag))

def double_cplusplus_format(a_mat):
    for i in range (a_mat.shape[0]):
        for j in range (a_mat.shape[1]):
            print("a({},{}) = {};".format(i,j,a_mat[i,j]))

# note that OhmmsMatrix is row major and by default so is numpy
# numpy and scipy cover this up for lapack calls
# as well as translate pivots to python 0 base indexing

# cuBLAS expects col major layouts
# at least some transposes that are unexplained in QMCPACK are about this.
def print_real_col_major(a_mat):
    for j in range (a_mat.shape[1]):
        for i in range(a_mat.shape[0]):
            print("{} ".format(a_mat[i,j]), end="")

def print_complex_col_major(a_mat):
    for j in range (a_mat.shape[1]):
        for i in range(a_mat.shape[0]):
            print("{}, {}, ".format(a_mat[i,j].real,a_mat[i,j].imag), end="")

# lapack and cuBLAS use 1 based indexing for pivots
def print_pivots(pivots):
    print("{{ {} }}".format(", ".join([str(p+1) for p in pivots])))


def math_log_det(lu_diag, pivots):
    log_value = complex(0.0, 0.0)
    for i in range(4):
        lud = lu_diag[i] * (1 -2 * (i != pivots[i]))
        log_value += complex(math.log(abs(lu_diag[i])), 0)
        log_value += complex(0,(lud < 0) * np.pi)
    return log_value

def cmath_log_det(lu_diag, pivot):
    log_value = complex(0.0,0.0)
    for i in range(4):
        log_value += cmath.log(lu_diag[i] * (1 - (2 * (pivot[i] != i))))
    return log_value

def com_cmath_log_det(lu_diag, pivot):
    log_value = complex(0.0,0.0)
    for i in range(4):
        pivot_factor = (1 - (2 * (pivot[i] != i)))
        diag = complex(lu_diag[i].real * pivot_factor, lu_diag[i].imag * pivot_factor)
        log_value += cmath.log(diag)
    return log_value


def com_math_log_det(lu_diag, pivot):
    log_value = complex(0.0,0.0)
    for i in range(4):
        pivot_factor = (1 - (2 * (pivot[i] != i)))
        diag = complex(lu_diag[i].real * pivot_factor, lu_diag[i].imag * pivot_factor)
        log_value += complex(math.log(math.sqrt(diag.real * diag.real + diag.imag * diag.imag)), 0)
        log_value += complex(0, math.atan2(diag.imag, diag.real))
        print(math.atan2(diag.imag, diag.real))
    return log_value


def pairwise(iterable):
    a, b = it.tee(iterable)
    next(b, None)
    return zip(a, b)

####################### Reference for cuBLAS_LU::computeLogDet_real
a_mat = np.array([[2, 5, 8, 7],
                  [5, 2, 2, 8],
                  [7, 5, 6, 6],
                  [5, 4, 4, 8]])

print("For computeLogDet real")
print("a_mat:")
print("{ ", end="")
print_real_col_major(a_mat)
print("}\n")

a_lu, a_piv = sl.lu_factor(a_mat);

print("a_lu_mat:")
print("{ ", end="")
print_real_col_major(a_lu)
print("}\n")

print("pivots:")
print_pivots(a_piv)

print("log_det:")
print(math_log_det(np.diag(a_lu),a_piv))

print("=======")

##################### Reference for cuBLAS_LU::computeLogDet_complex
complex_a_mat = np.array([[2.+0.1j, 5.+0.1j, 7.+0.2j, 5.+0.j ],
                          [5.+0.1j, 2.+0.2j, 5.+1.j , 4.-0.1j],
                          [8.+0.5j, 2.+0.1j, 6.-0.2j, 4.-0.6j],
                          [7.+1.j , 8.+0.5j, 6.-0.2j, 8.-2.j ]])

print("For direct cuBLAS test:")
print("complex_a_mat")
print("{ ", end="")
print_complex_col_major(complex_a_mat)
print("}\n")

com_LU, piv_LU = sl.lu_factor(complex_a_mat)

print("complex_lu_mat")
print("{ ", end="")
print_complex_col_major(com_LU)
print("}\n")

print("pivots:")
print_pivots(piv_LU)

print("log det:")
print(cmath_log_det(np.diag(com_LU),piv_LU))

print("=======")

####################### Reference for cuBLAS_LU::computeLogDet(batch=2)
print("batch=2 test")
print()

complex_a_mat = np.array([[2.+0.1j, 5.+0.1j, 7.+0.2j, 5.+0.j ],
                          [5.+0.1j, 2.+0.2j, 5.+1.j , 4.-0.1j],
                          [8.+0.5j, 2.+0.1j, 6.-0.2j, 4.-0.6j],
                          [7.+1.j , 8.+0.5j, 6.-0.2j, 8.-2.j ]])

complex_b_mat = np.array([[2.+0.1j, 2.+0.1j, 3.+0.2j, 5.+0.j ],
                          [4.+0.1j, 6.+0.2j, 5.+1.j , 4.-0.1j],
                          [8.+0.5j, 2.+0.1j, 3.-0.2j, 4.-0.6j],
                          [7.+1.j , 8.+0.5j, 6.-0.2j, 8.-2.j ]])

print ("complex_b_mat:")
print("{ ", end="")
print_complex_col_major(complex_b_mat)
print("}\n")

com_LU_b, piv_LU_b = sl.lu_factor(complex_b_mat)
p, l, u = sl.lu(complex_b_mat)
print("com b permutation ok?: ", np.allclose(complex_b_mat - p @ l @ u, np.zeros((4,4))))


print("com_LU_2:")
print_complex_col_major(com_LU_b)


lu_diag_1 = np.diag(com_LU);
lu_diag_2 = np.diag(com_LU_b);

upcom = lambda c: [c.real, c.imag]

print("diag1:")
print("{", end="")
[ print(" {}, {} ".format(*upcom(c)), end="") for c in lu_diag_1 ]
print("}")
print("diag2:")
print("{", end="")
[ print(" {}, {} ".format(*upcom(c)), end="") for c in lu_diag_2 ]
print("}")

print("pivots_a:")
print_pivots(piv_LU)
print("pivots_b:")
print_pivots(piv_LU_b)

pivots_1 = piv_LU
pivots_2 = piv_LU_b

print ("cmath complex log")
print ("log_det_a:")
print("{{ {}, {} }}".format(*upcom(com_cmath_log_det(lu_diag_1, pivots_1))))
print ("log_det_b:")
print("{{ {}, {} }}".format(*upcom(com_cmath_log_det(lu_diag_2, pivots_2))))
print()

print("math complex log")
print ("log_det_a:")
print("{{ {}, {} }}".format(*upcom(com_math_log_det(lu_diag_1, pivots_1))))
print ("log_det_b:")
print("{{ {}, {} }}".format(*upcom(com_math_log_det(lu_diag_2, pivots_2))))

      
