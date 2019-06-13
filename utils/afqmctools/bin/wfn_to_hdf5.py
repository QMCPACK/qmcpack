#!/usr/bin/env python

import argparse
import scipy.sparse
import sys
import time
from afqmctools.wavefunction.mol import write_qmcpack_wfn
from afqmctools.wavefunction.converter import read_qmcpack_ascii_wavefunction


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
                        type=str, default='wfn.h5',
                        help='Output file name for QMCPACK data.')
    parser.add_argument('-v', '--verbose', dest='verbose',
                        action='store_true', default=False,
                        help='Verbose output.')
    parser.add_argument('--nmo', dest='nmo', type=int, default=None)
    parser.add_argument('--nalpha', dest='nalpha', type=int, default=None)
    parser.add_argument('--nbeta', dest='nbeta', type=int, default=None)

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
    nelec = (nalpha, nbeta)
    wfn, walker_type = read_qmcpack_ascii_wavefunction(input_file, nmo, nelec)
    write_qmcpack_wfn(output_file, wfn, walker_type, nelec, nmo)

if __name__ == '__main__':

    main(sys.argv[1:])
