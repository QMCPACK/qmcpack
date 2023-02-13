#! /usr/bin/env python3

'''
Provide information about a multideterminant wavefunction stored in HDF5 format.
Can be used as a library or as CLI tool

determinants_tools.py ./tests/molecules/C2_pp/C2.h5    # Summary excitation info only
determinants_tools.py -v ./tests/molecules/C2_pp/C2.h5 # Verbose. Gives details of each determinant
'''

from __future__ import print_function
from builtins import zip

from functools import partial

# A detspin is a list of integers that represent the alpha or beta part of determinant
# A determinant is a pair of detspin alpha and detspin beta
# A CSF is represented as a pair of list of integers coding respectively for the double and single occupied state.

# A bitmask is a `n_state` length binary representation of a list of integers encoded in `bit_kind_size` bits.
#   Note that this representation is 'reversed'. This allows the first bit to represent the first orbital, and so on.
#   For performance and convenience, if `n_state` < `bit_kind_size` * `n_integer`,
#      we assume that the bits coding for excessive state of the last integer are set to 0.
#      To be clearer, if we have `n_state` = 1 with one alpha electron and `bit_kind_size` = 64
#      the only valid integer representation of this determinants is ( (1,), (0,) )

# Assume that each integer is unsigned. If not they should be converted before using the sanitize function.
import numpy as np
from collections import namedtuple
Determinants = namedtuple('Determinants', ['alpha', 'beta'])
CSF = namedtuple('CSF', ['double', 'single'])

def sanitize(det_spin, bit_kind_size):
     '''
     This function transform signed numpy variable to native "positive" python number
     In the case of :
        -positive Numpy number, we just cast it to Python data-type
        -negative Numpy number, we do the two-compelement trick. 
                                Because this value cannot fit on the initial Numpy datatype, 
                                it will silently casted to int.

     >>> import numpy as np
     >>> sanitize([np.int8(1)], 8)
     [1]
     >>> sanitize([np.int8(-1)], 8)
     [255]
     '''
     return [ (s + (1 << bit_kind_size)) if s < 0 else int(s) for s in det_spin]

def int_to_bitmask(s, bit_kind_size):
    '''
    >>> int_to_bitmask(1, 3)
    '100'
    >>> int_to_bitmask(2, 3)
    '010'
    >>> int_to_bitmask(3, 3)
    '110'
    >>> int_to_bitmask(7, 3)
    '111'
    '''
    return '{s:0{width}b}'.format(s=s, width=bit_kind_size)[::-1]


def detspin_to_bitmask(detspin, bit_kind_size, size=None):
    '''
    >>> detspin_to_bitmask([1,7], 3,6)
    '100111'
    >>> detspin_to_bitmask([1,2], 3,6)
    '100010'
    >>> detspin_to_bitmask([1,2], 3,5)
    '10001'
    '''
    return ''.join(int_to_bitmask(i, bit_kind_size) for i in detspin)[:size]

def det_excitation_degree(det_a, det_b):
    '''
    >>> det_excitation_degree( Determinants(alpha=[1],beta=[0] ),  
    ...                        Determinants(alpha=[1],beta=[0] ) )
    0
    >>> det_excitation_degree( Determinants(alpha=[1],beta=[0] ),  
    ...                        Determinants(alpha=[2],beta=[0] ) )
    1
    >>> det_excitation_degree( Determinants(alpha=[1,2],beta=[0] ),  
    ...                        Determinants(alpha=[2,1],beta=[0] ) )
    2
    >>> det_excitation_degree( Determinants(alpha=[7,0],beta=[0] ),  
    ...                        Determinants(alpha=[0,7],beta=[0] ))
    3
    '''
    or_and_popcount = lambda det : sum(bin(i ^ j).count("1") for i, j in zip(*det))
    return sum(map(or_and_popcount,zip(det_a, det_b))) // 2

def det_to_csf(det):
    '''
    >>> det_to_csf( ( ([1,0]),([2,0]) ) )
    CSF(double=(0, 0), single=(3, 0))
    '''
    # This function assume that determinant are zero-padded.
    d = tuple(a & b for a, b in zip(*det))
    s = tuple(a ^ b for a, b in zip(*det))
    return CSF(d, s)


def csf_to_string(csf, bit_kind_size, size):
    '''
    >>> csf_to_string( CSF(double=[1], single=[0]), 3, 3)
    '200'
    >>> csf_to_string( CSF(double=[0], single=[1]),  3, 3)
    '100'
    >>> csf_to_string( CSF(double=[1], single=[2]),  3, 3)
    '210'
    '''
    f = partial(detspin_to_bitmask, bit_kind_size=bit_kind_size, size=size)
    d_str, s_str = map(f, csf)

    r = ''
    for d, s in zip(d_str, s_str):
        if d != '0':
            r += '2'
        elif s != '0':
            r += '1'
        else:
            r += '0'
    return r

if __name__ == '__main__':
    import h5py
    import argparse

    from collections import Counter
    from itertools import chain

    parser = argparse.ArgumentParser(
        description='Provide information about a multideterminant wavefunction stored in HDF5 format.')

    parser.add_argument("h5path", type=str,
        help="Path to a h5file")

    parser.add_argument("-v", "--verbose",
        help="For each determinant, print its bitmask representation, associated CSF and excitation degree",
        action="store_true")

    parser.add_argument("-s", "--silent",
        help="Disable check", action="store_true")

    parser.add_argument("--transform",
        help="Transform a HDF5 that uses signed integers to unsigned integers.", action="store_true")

    args = parser.parse_args()

    # Run the unit test if needed
    if not args.silent:
        import doctest
        r = doctest.testmod(verbose=False)
        if r.failed:
            import sys
            sys.exit(1)

    open_mode = 'r'
    if args.transform:
        open_mode = 'a'
    with h5py.File(args.h5path, open_mode) as f:
        alpha_group = 'MultiDet/CI_Alpha'
        if alpha_group not in f:
            alpha_group = 'MultiDet/CI_0'
        beta_group = 'MultiDet/CI_Beta'
        if beta_group not in f:
            beta_group = 'MultiDet/CI_1'

        bit_kind_size = f[alpha_group].dtype.itemsize * 8

        # We sanitize the determinants. Aka we transform any negative integers, into the proper unsigned one.
        i_alpha= (sanitize(d,bit_kind_size) for d in f[alpha_group])
        i_beta= (sanitize(d,bit_kind_size) for d in f[beta_group])
        if args.transform:
            alpha_dets = np.array([i for i in i_alpha], dtype=np.uint64)
            del f[alpha_group]
            f.create_dataset(alpha_group, data=alpha_dets)
            beta_dets = np.array([i for i in i_beta], dtype=np.uint64)
            del f[beta_group]
            f.create_dataset(beta_group, data=beta_dets)
            # regenerate iterators
            i_alpha= (sanitize(d,bit_kind_size) for d in f[alpha_group])
            i_beta= (sanitize(d,bit_kind_size) for d in f[beta_group])

        if f[alpha_group].dtype == 'int64':
            raise TypeError("This HDF5 file uses signed 64 bit integers (int64) to represent the determinants. QMCPACK expects unsigned integers (uint64) this can cause issues when used in qmcpack.")

        i_det = zip(i_alpha, i_beta)

        # Getting the reference determinant to compute the excitation degree.
        hf_det = next(i_det)
        i_det = chain([hf_det], i_det) # Put back the original

        # Gathering information only needed in verbose mode
        if args.verbose:
            nstate = f['MultiDet/nstate'][0]
            to_string = partial(detspin_to_bitmask, bit_kind_size=bit_kind_size, size=nstate)

        # Aggregate the result and display curent information about the determinant
        c_e = Counter(); c_csf = Counter()
        for i, det in enumerate(i_det, 1):
            e = det_excitation_degree(det, hf_det); c_e[e] += 1
            c = det_to_csf(det); c_csf[c] += 1

            if args.verbose:
                str_alpha, str_beta = map(to_string, det)
                print(i)
                print('alpha ', str_alpha)
                print('beta  ', str_beta)
                print('scf   ', csf_to_string(c, bit_kind_size, nstate))
                print('excitation degree ', e)
                print('')

        # Display computed data
        n_det = sum(c_e.values()); n_csf = len(c_csf)

        print('Summary:')
        for e, i in sorted(c_e.items()):
            print('excitation degree', e, 'count:', i)

        print('')
        print('n_det', n_det)
        print('n_csf', n_csf)
