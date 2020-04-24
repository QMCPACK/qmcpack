#! /usr/bin/env python3

import argparse
import scipy.sparse
import sys
import time
from afqmctools.wavefunction.mol import write_qmcpack_wfn
from afqmctools.wavefunction.converter import (
        read_qmcpack_ascii_wavefunction,
        read_dmc_ci_wavefunction
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
                        default=None, help='Input wavefunction file.')
    parser.add_argument('-o', '--output', dest='output_file',
                        type=str, default='wfn.h5',
                        help='Output file name for wavefunction.')
    parser.add_argument('-v', '--verbose', dest='verbose',
                        action='store_true', default=False,
                        help='Verbose output.')
    parser.add_argument('--qmcpack', dest='qmcpack',
                        action='store_true', default=False,
                        help='Convert from DMC CI wavefunction.')
    parser.add_argument('--nmo', dest='nmo', type=int, default=None,
                        help='Number of orbitals.')
    parser.add_argument('--nalpha', dest='nalpha', type=int,
                        default=None, help='Number of alpha electrons.')
    parser.add_argument('--nbeta', dest='nbeta', type=int,
                        default=None, help='Number of alpha electrons.')
    parser.add_argument('--ndets', dest='ndets', type=int,
                        default=None, help='Number of determinants to write.')

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
    input_file = options.input_file
    output_file = options.output_file
    nmo = options.nmo
    nalpha = options.nalpha
    nbeta = options.nbeta
    ndets = options.ndets
    nelec = (nalpha, nbeta)
    if options.qmcpack:
        wfn, walker_type, nmo, nelec = (
                read_dmc_ci_wavefunction(input_file, nelec, nmo, options.ndets)
                )
    else:
        wfn, walker_type = read_qmcpack_ascii_wavefunction(input_file, nmo, nelec)
    write_qmcpack_wfn(output_file, wfn, walker_type, nelec, nmo)

if __name__ == '__main__':

    main(sys.argv[1:])
