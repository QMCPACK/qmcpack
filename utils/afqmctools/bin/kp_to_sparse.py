#! /usr/bin/env python3

import argparse
import sys
from afqmctools.hamiltonian.converter import (
        read_qmcpack_hamiltonian,
        kpoint_to_sparse
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
                        default=None, help='Input kpoint factorized file.')
    parser.add_argument('-o', '--output', dest='output_file',
                        type=str, default='sparse.h5',
                        help='Output file for sparse hamiltonian.')
    parser.add_argument('-r', '--real-chol', dest='real_chol',
                        action='store_true', default=False,
                        help='Dump real integrals.')
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
    kpoint_to_sparse(options.input_file,
                     options.output_file,
                     real_chol=options.real_chol,
                     verbose=options.verbose)


if __name__ == '__main__':
    main(sys.argv[1:])
