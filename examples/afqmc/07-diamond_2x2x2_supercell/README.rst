Example 7: 2x2x2 Diamond supercell
----------------------------------

In this example we will show how to generate the AFQMC input from a pbc pyscf calculation
for a 2x2x2 supercell of diamond using a RHF trial wavefunction.

Again the first step is to run a pyscf calculation using the scf.py script in this
directory.

The first part of the pyscf calculation is straightforward. See pyscf/examples/pbc for
more examples on how to set up Hartree--Fock and DFT simulations.


.. code-block:: python

    import h5py
    import numpy
    import sys
    from pyscf.pbc import scf, dft, gto

    cell = gto.Cell()
    cell.verbose = 5
    alat0 = 3.6
    cell.a = (numpy.ones((3,3))-numpy.eye(3))*alat0 / 2.0
    cell.atom = (('C',0,0,0), ('C',numpy.array([0.25,0.25,0.25])*alat0))
    cell.basis = 'gth-szv'
    cell.pseudo = 'gth-pade'
    cell.mesh = [28,28,28]
    cell.build()
    nk = [2,2,2]
    kpts = cell.make_kpts(nk)

    mf = scf.KRHF(cell, kpts=kpts)
    mf.chkfile = 'scf.chk'
    mf.kernel()


In addition to a standard pyscf calculation, we add the following lines:

.. code-block:: python

    from afqmctools.utils.linalg import get_ortho_ao
    hcore = mf.get_hcore()
    fock = (hcore + mf.get_veff())
    X, nmo_per_kpt = get_ortho_ao(cell,kpts)
    with h5py.File(mf.chkfile) as fh5:
      fh5['scf/hcore'] = hcore
      fh5['scf/fock'] = fock
      fh5['scf/orthoAORot'] = X
      fh5['scf/nmo_per_kpt'] = nmo_per_kpt

essentially, storing the fock matrix, core Hamiltonian and transformation matrix to the
orthogonalised AO basis. This is currently required for running PBC AFQMC calculations.

Once the above (scf.py) script is run we will again use the `pyscf_to_afqmc.py` script
to generate the necessary AFQMC input file.

.. code-block:: bash

    mpirun -n 8 /path/to/qmcpack/utils/afqmctools/bin/pyscf_to_afqmc.py -i scf.chk -o afqmc.h5 -t 1e-5 -v -a


Note that the commands necessary to generate the integrals are identical to those for the
molecular calculations, except now we accelerate their calculation using MPI. Note that if
the number of tasks > number of kpoints then the number of MPI tasks must be divisible by
the number of kpoints.

Once this is done we will again find a Hamiltonian file and afqmc input xml file.
Inspecting these you will notice that their structure is identical to the molecular
calculations seen previously. This is because we have not exploited k-point symmetry and
are writing the integrals in a supercell basis. In the next example we will show how
exploiting k-point symmetry can be done explicitly, which leads to a faster and lower
memory algorithm for AFQMC.
