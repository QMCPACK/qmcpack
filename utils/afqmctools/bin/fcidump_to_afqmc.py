#!/usr/bin/env python

import argparse
import scipy.sparse
import sys
import time
import numpy
from afqmctools.hamiltonian.mol import (
        read_ascii_integrals,
        modified_cholesky_direct,
        write_qmcpack_cholesky
        )


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
    parser.add_argument('-c', '--complex', dest='complex_chol',
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
    (hcore, eri, ecore, nelec) = read_ascii_integrals(options.input_file,
                                                      symmetry=options.symm,  
                                                      verbose=options.verbose)
    nelec = (nelec//2,nelec//2)
    norb = hcore.shape[-1]
    
    if options.symm == 4: # assuming complex
        eri = numpy.transpose(eri,(0,1,3,2))

    chol = modified_cholesky_direct(eri.reshape(norb**2,norb**2),
                                    options.thresh, options.verbose).T
    write_qmcpack_cholesky(hcore, scipy.sparse.csr_matrix(chol),
                           nelec, norb, e0=ecore,
                           real_chol=(not options.complex_chol),
                           filename=options.output_file)

if __name__ == '__main__':

    main(sys.argv[1:])
