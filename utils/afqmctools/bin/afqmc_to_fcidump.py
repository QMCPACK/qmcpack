#! /usr/bin/env python3

import argparse
import sys
from afqmctools.hamiltonian.converter import (
        read_qmcpack_hamiltonian,
        write_fcidump,
        write_fcidump_kpoint
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
                        default=None, help='Input AFQMC hamiltonian file.')
    parser.add_argument('-o', '--output', dest='output_file',
                        type=str, default='FCIDUMP',
                        help='Output file for FCIDUMP.')
    parser.add_argument('-s', '--symmetry', dest='symm',
                        type=int, default=1,
                        help='Symmetry of integral file (1,4,8).')
    parser.add_argument('-t', '--tol', dest='tol',
                        type=float, default=1e-12,
                        help='Cutoff for integrals.')
    parser.add_argument('-c', '--complex', dest='cplx',
                        action='store_true', default=False,
                        help='Whether to write integrals as complex numbers.')
    parser.add_argument('--complex-paren', dest='cplx_paren',
                        action='store_true', default=False,
                        help='Whether to write FORTRAN format complex numbers.')
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
    hamil = read_qmcpack_hamiltonian(options.input_file)
    if hamil is None:
        sys.exit()

    if hamil.get('qk_k2') is None:
        write_fcidump(options.output_file,
                      hamil['hcore'],
                      hamil['chol'],
                      hamil['enuc'],
                      hamil['nmo'],
                      hamil['nelec'],
                      tol=options.tol,
                      sym=options.symm,
                      cplx=options.cplx,
                      paren=options.cplx_paren)
    else:
        write_fcidump_kpoint(options.output_file,
                             hamil['hcore'],
                             hamil['chol'],
                             hamil['enuc'],
                             hamil['nmo'],
                             hamil['nelec'],
                             hamil['nmo_pk'],
                             hamil['qk_k2'],
                             tol=options.tol,
                             sym=options.symm,
                             cplx=options.cplx,
                             paren=options.cplx_paren)

if __name__ == '__main__':
    main(sys.argv[1:])
