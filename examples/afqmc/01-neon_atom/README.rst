Example 1: Neon atom
--------------------

In this example we will go through the basic steps necessary to
generate AFQMC input from a pyscf scf calculation on a simple closed
shell molecule (neon/aug-cc-pvdz).

The pyscf scf script is given below (scf.py in the current directory):

.. code-block:: python

    from pyscf import gto, scf, cc
    from pyscf.cc import ccsd_t
    import h5py

    mol = gto.Mole()
    mol.basis = 'aug-cc-pvdz'
    mol.atom = (('Ne', 0,0,0),)
    mol.verbose = 4
    mol.build()

    mf = scf.RHF(mol)
    mf.chkfile = 'scf.chk'
    ehf = mf.kernel()

    ccsd = cc.CCSD(mf)
    ecorr_ccsd = ccsd.kernel()[0]
    ecorr_ccsdt = ccsd_t.kernel(ccsd, ccsd.ao2mo())
    print("E(CCSD(T)) = {}".format(ehf+ecorr_ccsd+ecorr_ccsdt))

The most important point above is that we create a scf checkpoint file by specifying the
`mf.chkfile` mol member variable. Note we will also compute the CCSD and CCSD(T) energies
for comparison puposes since this system is trivially small.

We next run the pyscf calculation using

.. code-block:: bash

    python scf.py > scf.out

which will yield a converged restricted Hartree--Fock total energy of -128.496349730541
Ha, a CCSD value of -128.7084878405062 Ha, and a CCSD(T) value of -128.711294157 Ha.

The next step is to generate the necessary qmcpack input from this scf calculation. To
this we do (assuming afqmctools is in your PYTHONPATH):

.. code-block:: bash

    /path/to/qmcpack/utils/afqmctools/bin/pyscf_to_afqmc.py -i scf.chk -o afqmc.h5 -t 1e-5 -v

which will peform the necessary AO to MO transformation of the one and two electron
integrals and perform a modified cholesky transormation of the two electron integrals. A
full explanation of the various options available for `pyscf_to_afqmc.py` you can do

.. code-block:: bash

    pyscf_to_afqmc.py -h

In the above example, `-i` designates the input pyscf checkpoint file, `-o` speficies the
output filename to write the qmcpack hamiltonian/wavefunction to, `-t` specifies the
convergence threshold for the Cholesky decomposition, `-v` increases verbosity.
You can optionally pass the `-q/--qmcpack-input` to generate a qmcpack input
file which is based on the hamiltonian and wavefunction generated. Greater control over
input file generation can be achieved using the write_xml_input function provided with
afqmctools. Run gen_input.py after the integrals/wavefunction have been generated  to
generate the input file `afqmc.xml`.

Running the above will generate one file: `afqmc.h5`. The plain text wavefunction files
are deprecated and will be removed in later releases. The qmcpack input file `afqmc.xml`
is a *skeleton* input file, meaning that it's created from the information in `hamil.h5`
and is meant as a convenience, not as a guarantee that the convergeable parameters
(timestep, walker number, bias bound etc. are converged or appropriate).

We will next run through the relevant sections of the input file `afqmc.xml` below:

.. code-block:: xml

        <project id="qmc" series="0"/>
        <random seed="7"/>

        <AFQMCInfo name="info0">
            <parameter name="NMO">23</parameter>
            <parameter name="NAEA">5</parameter>
            <parameter name="NAEB">5</parameter>
        </AFQMCInfo>

We first specify how to name the output file. We also have fixed the random number seed so
that the results of this tutorial can be reproducible (if run on the same number of
cores).

Next comes the system description, which is mostly a sanity check, as these parameters
will be read from the hamiltonian file. They specify the number of single-particle
orbitals in the basis set (`NMO`) and the number of alpha (`NAEA`) and beta (`NAEB`)
electrons respectively.

Next we specify the Hamiltonian and wavefunction to use:

.. code-block:: xml

        <Hamiltonian name="ham0" info="info0">
          <parameter name="filetype">hdf5</parameter>
          <parameter name="filename">afqmc.h5</parameter>
        </Hamiltonian>

        <Wavefunction name="wfn0" type="NOMSD" info="info0">
          <parameter name="filetype">hdf5</parameter>
          <parameter name="filename">afqmc.h5</parameter>
        </Wavefunction>

The above should be enough for most calculations. A `NOMSD` (non-orthogonal multi-Slater
determinant) wavefunction allows for a generalised wavefunction input in the form of a
single (or multiple) matrix (matrices) of molecular orbital coefficients for the RHF
calculation we perform here.

We next set the walker options:

.. code-block:: xml

        <WalkerSet name="wset0" type="shared">
          <parameter name="walker_type">CLOSED</parameter>
        </WalkerSet>

The important point here is that as we are using a RHF trial wavefunction we must specify
that the `walker_type` is `CLOSED`. For a UHF trial wavefunction one would set this to
`COLLINEAR`.

And now the propagator options:

.. code-block:: xml

        <Propagator name="prop0" info="info0">
          <parameter name="hybrid">yes</parameter>
        </Propagator>

In the above we specify that we will be using the hybrid approach for updating the walker
weights. If you wish to use the local energy approximation you should set this flag to
false.

Finally comes the execute block which controls how the simulation is run:

.. code-block:: xml

        <execute wset="wset0" ham="ham0" wfn="wfn0" prop="prop0" info="info0">
          <parameter name="ncores">1</parameter>
          <parameter name="timestep">0.01</parameter>
          <parameter name="nWalkers">10</parameter>
          <parameter name="blocks">100</parameter>
          <parameter name="steps">10</parameter>
       </execute>

The time step (`timestep`), number of Monte Carlo samples (`blocks`*`steps`), and number
of walkers (`nWalkers`) should be adjusted as appropriate. Note that `nWalkers` sets the
number of walkers per `ncores`. For example, if we wanted to use 100 walkers we could run
the above input file on 10 cores. If the problem size is very large we may want
distribute the workload over more cores per walker, say 10. In this case we would require
100 cores to maintain the same number of walkers. Typically in this case you want to
specify fewer walkers per core anyway.

We can now run the qmcpack simulation:

.. code-block:: bash

    qmcpack afqmc.xml > qmcpack.out

Assuming the calculation finishes successfully, the very first thing you should do is
check the information in `qmcpack.out` to see confirm no warnings were raised.  The second
thing you should check is that the energy of the starting determinant matches the
Hartree--Fock energy you computed earlier from pyscf to within roughly the error threshold
you specified when generating the Cholesky decomposition. This check is not very
meaningful if using, say, DFT orbitals.  However if this energy is crazy it's a good sign
something went wrong with either the wavefunction or integral generation.  Next you should
inspect the `qmc.s000.scalar.dat` file which contains the mixed estimates for various
quantities. This can be plotted using gnuplot.  `EnergyEstim__nume_real` contains the
block averaged values for the local energy, which should be the 7th column.

Assuming everything worked correctly we need to analyse the afqmc output using:

.. code-block:: bash

    /path/to/qmcpack/nexus/bin/qmca -e num_skip -q el qmc.s000.scalar.dat

where `num_skip` is the number of blocks to skip for the equilibration stage. For a
practical calculation you may want to use more walkers and run for longer to get
meaningful statistics.

See the options for qmca for further information. Essentially we discarded the first 100
blocks as equilibaration and only computed the mixed estimate for the local energy
internally called `EnergyEstim__nume_real`, which can be specified with `-q el`. We see
that the ph-AFQMC energy agrees well with the CCSD(T) value. However, we probably did not
run the simulation for long enough to really trust the error bars.
