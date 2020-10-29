Example 5: CASSCF Trial
-----------------------

In this example we will show how to format a casscf trial wavefunction.

Rather than use the `pyscf_to_afqmc.py`, script we will break up the process to allow
for more flexibility and show what is going on under the hood.

The qmcpack input can be generated with the scf.py script followed by gen_input.py.

See the relevant code below for a breakdown of the steps involved.

The first step is to run a CASSCF calculation. Here we'll consider N:sub:2. This
replicates the calculations from Al-Saidi et al J. Chem. Phys. 127, 144101 (2007).
They find a CASSCF energy of -108.916484 Ha, and a ph-AFQMC energy of -109.1975(6) Ha with
a 97 determinant CASSCF trial.

.. code-block:: python

    mol = gto.M(atom=[['N', (0,0,0)], ['N', (0,0,3.0)]],
                basis='cc-pvdz',
                unit='Bohr')
    nalpha, nbeta = mol.nelec
    rhf = scf.RHF(mol)
    rhf.chkfile = 'scf.chk'
    rhf.kernel()

    M = 12
    N = 6
    nmo = rhf.mo_coeff.shape[-1]
    mc = mcscf.CASSCF(rhf, M, N)
    mc.chkfile = 'scf.chk'
    mc.kernel()

Next we unpack the wavefunction

.. code-block:: python

    nalpha = 3
    nbeta = 3
    ci, occa, occb = zip(*fci.addons.large_ci(mc.ci, M, (nalpha,nbeta),
                         tol=tol, return_strs=False))

and sort the determinants by the magnitude of their weight:

.. code-block:: python

    ixs = numpy.argsort(numpy.abs(coeff))[::-1]
    coeff = coeff[ixs]
    occa = numpy.array(occa)[ixs]
    occb = numpy.array(occb)[ixs]

Next we reinsert the frozen core as the AFQMC simulation is not run using an active space:

.. code-block:: python

    core = [i for i in range(mc.ncore)]
    occa = [numpy.array(core + [o + mc.ncore for o in oa]) for oa in occa]
    occb = [numpy.array(core + [o + mc.ncore for o in ob]) for ob in occb]

Next we need to generate the one- and two-electron integrals. Note that we need to use the
CASSCF MO coefficients to rotate the integrals.

.. code-block:: python

    scf_data = load_from_pyscf_chk_mol('scf.chk', 'mcscf')
    write_hamil_mol(scf_data, 'afqmc.h5', 1e-5, verbose=True)

Finally we can write the wavefunction to the QMCPACK format:

.. code-block:: python

    ci = numpy.array(ci, dtype=numpy.complex128)
    uhf = True # UHF always true for CI expansions.
    write_qmcpack_wfn('afqmc.h5', (ci, occa, occb), uhf, mol.nelec, nmo)

To generate the input file we again run ``gen_input.py``. Note the ``rediag`` option which
is necessary if the CI code used uses a different convention for ordering creation and
annihilation operations when defining determinant strings.
