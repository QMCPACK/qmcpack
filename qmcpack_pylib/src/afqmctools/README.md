# AFQMCTOOLS

A collection of python tools to aid in the running and analysis of the AFQMC simulations
in QMCPACK.

The afqmctools library contains useful functions to help with generating the input for
afqmc simulations.  Currently most of these routines depend on pyscf, but adaptations for
other software packages should be straightforward.

In addition to the library there are two useful scripts to automate the generation of
simulation input:

1. afqmctools/bin/pyscf_to_afqmc.py which will generate the necessary AFQMC QMCPACK
   input from a pyscf scf calculation.
2. afqmctools/bin/fcidump_to_afqmc.py which will generate an AFQMC QMCPACK Hamiltonian
   from a plain text FCIDUMP. The integrals are assumed to be real and 8-fold symmetric.

See the examples in qmcpack/examples/afqmc for more details on using these scripts or pass
-h/--help to the scripts themselves.

You will have to add afqmctools to your PYTHONPATH.

# Requirements

The tools work with the following:

* python >= 2.7 or python > 3.6
* pyscf >= 1.6.0
* scipy >= 0.18.1
* numpy >= 1.11.2
* h5py >= 2.6.0 with parallel hdf5 support for k-point symmetric integral generation
  (optional).
* mpi4py >= 2.0.0

# Tests

To run the unit tests do

```
python -m unittest discover -v
```

in the top level of the afqmctools directory. The tests should take less than three
minutes.
