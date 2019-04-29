# AFQMCTOOLS

A collection of python scripts to aid in the running of the AFQMC code in QMCPACK.

There are two main scripts right now.

1. afqmctools/bin/pyscf_to_qmcpack.py which will generate the necessary AFQMC QMCPACK
   input from a pyscf scf calculation.
2. afqmctools/bin/fcidump_to_qmcpack.py which will generate an AFQMC QMCPACK Hamiltonian
   from a plain text FCIDUMP. The integrals are assumed to be real and 8-fold symmetric.

In addition the library afqmctools contains many useful functions to help with
generating the input for more complicated simulation workflows.  Currently most of these
routines are only useful in conjuction with pyscf, but adaptations for other software
packages should be straightforward.

See the examples in qmcpack/examples/afqmc for more details on using these scripts.

# Requirements

The tools work with the following:

* python >= 2.7 or (preferably) python > 3.6
* pyscf >= 1.6.0
* scipy >= 0.18.1
* numpy >= 1.11.2
* h5py >= 2.6.0 with parallel hdf5 support for k-point symmetric integral generation.
* pandas >= 0.18.0 (optional)
