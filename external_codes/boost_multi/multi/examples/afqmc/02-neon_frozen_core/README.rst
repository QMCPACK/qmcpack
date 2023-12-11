Example 2: Frozen Core
----------------------

In this example we show how to perform a frozen core calculation, which only affects the
integral generation step. We will use the the previous Neon example and freeze 2 core
electrons. The following only currently works for RHF/ROHF trial wavefunctions.

.. code-block:: bash

    mpirun -n 1 /path/to/qmcpack/utils/afqmctools/bin/pyscf_to_afqmc.py -i scf.chk -o afqmc.h5 -t 1e-5 -v -c 8,22

Again, run gen_input.py to generate the input file `afqmc.xml`.

Comparing the above to the previous example we see that we have added the `-c` or `--cas`
option followed by a comma separated list of the form N,M defining a (N,M) CAS space
containing 8 electrons in 22 spatial orbitals (freezing the lowest MO).

The rest of the qmcpack process follows as before.
