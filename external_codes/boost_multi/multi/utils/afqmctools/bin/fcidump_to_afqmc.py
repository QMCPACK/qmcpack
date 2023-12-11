#! /usr/bin/env python3

import argparse
import scipy.sparse
import sys
import time
import numpy
from afqmctools.hamiltonian.mol import (
        write_qmcpack_sparse
        )
from afqmctools.hamiltonian.converter import read_fcidump
from afqmctools.utils.linalg import modified_cholesky_direct


def parse_args(args):
    """Parse command-line arguments.

    Parameters
    ----------
    args : list of strings
        command-line arguments.

    Returns
    -------
    options : :class:`argparse.ArgumentParser`
        Command line arguments.
    """

    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument('-i', '--input', dest='input_file', type=str,
                        default=None, help='Input FCIDUMP file.')
    parser.add_argument('-o', '--output', dest='output_file',
                        type=str, default='fcidump.h5',
                        help='Output file name for PAUXY data.')
    parser.add_argument('--write-complex', dest='write_complex',
                        action='store_true', default=False,
                        help='Output integrals in complex format.')
    parser.add_argument('-t', '--cholesky-threshold', dest='thresh',
                        type=float, default=1e-5,
                        help='Cholesky convergence threshold.')
    parser.add_argument('-s', '--symmetry', dest='symm',
                        type=int, default=8,
                        help='Symmetry of integral file (1,4,8).')
    parser.add_argument('-v', '--verbose', dest='verbose',
                        action='store_true', default=False,
                        help='Verbose output.')

    options = parser.parse_args(args)

    if not options.input_file:
        parser.print_help()
        sys.exit(1)

    return options

def main(args):
    """Convert FCIDUMP to QMCPACK readable Hamiltonian format.

    Parameters
    ----------
    args : list of strings
        command-line arguments.
    """
    options = parse_args(args)
    (hcore, eri, ecore, nelec) = read_fcidump(options.input_file,
                                              symmetry=options.symm,
                                              verbose=options.verbose)
    norb = hcore.shape[-1]

    # If the ERIs are complex then we need to form M_{(ik),(lj}} which is
    # Hermitian. For real integrals we will have 8-fold symmetry so trasposing
    # will have no effect.
    eri = numpy.transpose(eri,(0,1,3,2))

    chol = modified_cholesky_direct(eri.reshape(norb**2,norb**2),
                                    options.thresh, options.verbose,
                                    cmax=20).T.copy()
    cplx_chol = options.write_complex or numpy.any(abs(eri.imag)>1e-14)
    write_qmcpack_sparse(hcore, chol, nelec, norb, e0=ecore,
                         real_chol=(not cplx_chol),
                         filename=options.output_file)

if __name__ == '__main__':
    main(sys.argv[1:])
