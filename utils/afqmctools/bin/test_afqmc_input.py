#!/usr/bin/env python3

import argparse
import sys

from afqmctools.hamiltonian.converter import read_qmcpack_hamiltonian

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
    parser.add_argument('-i', '--input', dest='input', type=str,
                        default=None, help='Input afqmc file.')

    options = parser.parse_args(args)

    if not options.input:
        parser.print_help()
        sys.exit()

    return options

def main(args):
    """Sanity check for afqmc Hamiltonian.

    Parameters
    ----------
    args : list of strings
        command-line arguments.
    """
    options = parse_args(args)
    hamil = read_qmcpack_hamiltonian(options.input)
    if hamil is None:
        sys.exit(1)
    nerror = 0
    for k, v in hamil.items():
        if v is None:
            nerror += 1
    if nerror > 0:
        print("Found {:} non fatal error reading Hamiltonian file.".format(nerror))
        sys.exit(1)
    else:
        sys.exit(0)

if __name__ == '__main__':
    main(sys.argv[1:])
