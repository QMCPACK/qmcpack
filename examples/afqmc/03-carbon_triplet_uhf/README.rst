Example 3: UHF Trial
--------------------

In this example we show how to use a unrestricted Hartree--Fock (UHF) style wavefunction
to find the ph-AFQMC (triplet) ground state energy of the carbon atom (cc-pvtz). Again we
first run the scf (scf.py) calculation followed by the integral generation script:

.. code-block:: bash

    mpirun -n 1 /path/to/qmcpack/utils/afqmctools/bin/pyscf_to_afqmc.py -i scf.chk -o afqmc.h5 -t 1e-5 -v -a

Note the new flag `-a`/`--ao` which tells the script to transform the integrals to an
orthogonalised atomic orbital basis, rather that the more typical MO basis. This is
necessary as qmcpack does not support spin dependent two electron integrals.

Running qmcpack as before should yield a mixed estimate for the energy of roughly: -37.78471 +/- 0.00014.
