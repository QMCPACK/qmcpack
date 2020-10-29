#! /usr/bin/env python3

import argparse
import h5py
import json
from mpi4py import MPI
import numpy
import os
import scipy.sparse
import sys
import time
from afqmctools.inputs.from_pyscf import write_qmcpack
from afqmctools.utils.misc import get_git_hash

def parse_args(args, comm):
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

    if comm.rank == 0:
        parser = argparse.ArgumentParser(description = __doc__)
        parser.add_argument('-i', '--input', dest='chk_file', type=str,
                            default=None, help='Input pyscf .chk file.')
        parser.add_argument('-o', '--output', dest='hamil_file',
                            type=str, default='hamiltonian.h5',
                            help='Output file name for QMCPACK hamiltonian.')
        parser.add_argument('-w', '--wavefunction', dest='wfn_file',
                            type=str, default=None,
                            help='Output file name for QMCPACK wavefunction. '
                                 'By default will write to hamil_file.')
        parser.add_argument('-q', '--qmcpack-input', dest='qmc_input',
                            type=str, default=None,
                            help='Generate skeleton QMCPACK input xml file.')
        parser.add_argument('-t', '--cholesky-threshold', dest='thresh',
                            type=float, default=1e-5,
                            help='Cholesky convergence threshold.')
        parser.add_argument('-k', '--kpoint', dest='kpoint_sym',
                            action='store_true', default=False,
                            help='Generate explicit kpoint dependent integrals.')
        parser.add_argument('--density-fit', dest='df',
                            action='store_true', default=False,
                            help='Use density fitting integrals stored in '
                            'input pyscf chkpoint file.')
        parser.add_argument('-a', '--ao', '--ortho-ao', dest='ortho_ao',
                            action='store_true', default=False,
                            help='Transform to ortho AO basis. Default assumes '
                            'we work in MO basis')
        parser.add_argument('-c', '--cas',
                            help='Specify a CAS in the form of N,M.',
                            type=lambda s: [int(item) for item in s.split(',')],
                            default=None)
        parser.add_argument('-d', '--disable-ham', dest='disable_ham',
                            action='store_true', default=False,
                            help='Disable hamiltonian generation.')
        parser.add_argument('-n', '--num-dets', dest='ndet_max',
                            type=int, default=1,
                            help='Set upper limit on number of determinants to '
                            'generate.')
        parser.add_argument('-r', '--real-ham', dest='real_chol',
                            action='store_true', default=False,
                            help='Write integrals as real numbers.')
        parser.add_argument('-p', '--phdf', dest='phdf',
                            action='store_true', default=False,
                            help='Use parallel hdf5.')
        parser.add_argument('--low', dest='low_thresh',
                            type=float, default=0.1,
                            help='Lower threshold for non-integer occupancies'
                            'to include in multi-determinant exansion.')
        parser.add_argument('--high', dest='high_thresh',
                            type=float, default=0.95,
                            help='Upper threshold for non-integer occupancies'
                            'to include in multi-determinant exansion.')
        parser.add_argument('--dense', dest='dense',
                            action='store_true', default=False,
                            help='Write dense Hamiltonian.')
        parser.add_argument('-v', '--verbose', action='count', default=0,
                            help='Verbose output.')

        options = parser.parse_args(args)
    else:
        options = None
    options = comm.bcast(options, root=0)

    if not options.chk_file:
        if comm.rank == 0:
            parser.print_help()
        sys.exit()

    return options

def write_metadata(options, sha1, cwd, date_time):
    op_dict = vars(options)
    op_dict['git_hash'] = sha1
    op_dict['working_directory'] = os.getcwd()
    op_dict['date_time'] = date_time
    if options.wfn_file != options.hamil_file:
        with h5py.File(options.wfn_file, 'a') as fh5:
            try:
                fh5['metadata'] = json.dumps(op_dict)
            except RuntimeError:
                del fh5['metadata']
                fh5['metadata'] = json.dumps(op_dict)
    if not options.disable_ham:
        with h5py.File(options.hamil_file, 'a') as fh5:
            fh5['metadata'] = json.dumps(op_dict)

def main(args):
    """Generate QMCPACK input from pyscf checkpoint file.

    Parameters
    ----------
    args : list of strings
        command-line arguments.
    """
    comm = MPI.COMM_WORLD
    options = parse_args(args, comm)
    if comm.rank == 0:
        cwd = os.getcwd()
        sha1 = get_git_hash()
        date_time = time.asctime()
        print(" # Generating QMCPACK input from PYSCF checkpoint file.")
        print(" # git sha1: {}".format(sha1))
        print(" # Date/Time: {}".format(date_time))
        print(" # Working directory: {}".format(cwd))

    if options.wfn_file is None:
        if options.disable_ham:
            options.wfn_file = 'wfn.dat'
        else:
            options.wfn_file = options.hamil_file
    write_qmcpack(options.chk_file, options.hamil_file, options.thresh,
                  comm=comm,
                  ortho_ao=options.ortho_ao,
                  kpoint=options.kpoint_sym, df=options.df,
                  verbose=options.verbose, cas=options.cas,
                  qmc_input=options.qmc_input,
                  wfn_file=options.wfn_file,
                  write_hamil=(not options.disable_ham),
                  ndet_max=options.ndet_max,
                  real_chol=options.real_chol,
                  phdf=options.phdf,
                  low=options.low_thresh,
                  high=options.high_thresh,
                  dense=options.dense)
    if comm.rank == 0:
        write_metadata(options, sha1, cwd, date_time)

    if comm.rank == 0:
        print("\n # Finished.")

if __name__ == '__main__':

    main(sys.argv[1:])
