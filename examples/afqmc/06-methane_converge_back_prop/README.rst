Example 6: Back Propagation
---------------------------

.. note::
    matplotlib is required to generate the figure in this example.

The basic estimators printed out in the qmcpack `*.scalar.dat` files are *mixed*
estimates. Unless the operator for which the mixed estimate is computed commutes with the
Hamiltonian this result will generally be biased. To obtain pure estimates we can use
back propagation as outlined in: `Motta & Zhang, JCTC 13, 5367 (2017)`. For this example
we will look at computing the one-body energy of a methane molecule (see Fig. 2 of M&Z).

As before run scf.py and generate the integrals using `pyscf_to_afmqc.py`:

.. code-block:: bash

    mpirun -n 1 /path/to/qmcpack/utils/afqmctools/bin/pyscf_to_afqmc.py -i scf.chk -o afqmc.h5 -t 1e-5 -v

Note we are working in the MO basis. The input file is generated using `gen_input.py` and
comparing to the previous examples we can now see the estimator block:

.. code-block:: xml

      <Estimator name="back_propagation">
          <parameter name="naverages">4</parameter>
          <parameter name="block_size">2</parameter>
          <parameter name="ortho">1</parameter>
          <OneRDM />
          <parameter name="nsteps">200</parameter>
      </Estimator>

Which will tell QMCPACK to compute the back propagated one-rdm.  In the above we set
`block_size` to be 2 meaning that we average the back propagated estimates into bins of
length 2 in this case. This helps reduce the size of the hdf5 files.  We also specify the
option `nsteps`: We see that it is set to 200, meaning that we will back propagated the
bra wavefunction in the estimator by 200*.01 = 2 a.u., where the timestep has been set to
0.01 a.u. Finally `naverages` allows us to split the full path into `naverages` chunks,
so we will have averaged data at :math:`\tau_{BP}=[0.5, 1.0, 1.5, 2.0]` au.
This allows us to monitor the convergence of the estimator with back propagation time.


Running QMCPACK as before we will notice that in addition to the `qmc.s000.scalar.dat`
file we have generated a new file `qmc.s000.scalar.h5`. This file will contain the back
propagated estimates, which, for the time being, means the back propagated one-particle
reduced density matrix (1RDM), given as

.. math::

    P^{\sigma}_{ij} = \langle c_{i\sigma}^{\dagger} c_{j\sigma} \rangle

Before we analyse the output we should question why we chose a back propagation time of 2
au.  The back propagation time represents yet another parameter which must be carefully
converged.

In this example we will show how this is done.  In this directory you will find a script
`check_h1e_conv.py` which shows how to use various helper scripts provided in
`afqmctools/analysis/average.py`. The most of important of which are:

.. code-block:: python

    from afqmctools.analysis.extraction import get_metadata
    metadata = get_metadata(filename)

which returns a dict containing the RDM metadata,

.. code-block:: python

    from afqmctools.analysis.average import average_one_rdm
    rdm_av, rdm_errs = average_one_rdm(f, name='back_propagated', eqlb=3, ix=2)

which computes the average of the 1RDM, where 'i' specifies the index for the length of
back propagation time desired (e.g. :math:`i=2 \rightarrow \tau_{BP} = 1.5` au). `eqlb` is
the equilibration time, and here we skip 10 blocks of length 2 au.

.. code-block:: python

    from afqmctools.analysis.extraction import extract_observable
    dm = extract_observable(filename,
                            estimator='back_propagated'
                            name='one_rdm',
                            ix=2)

which extracts the 1RDM for all blocks and finally,

.. code-block:: python

    from afqmctools.analysis.extraction import extract_observable
    dm, weights = extract_observable(filename,
                                     estimator='back_propagated'
                                     name='one_rdm',
                                     ix=2,
                                     sample=index)

which extracts a single density matrix for block `index`.

Have a look through `check_h1e_conv.py` and run it. A plot should be produced which shows
the back propagated AFQMC one-body energy as a function of back propagation time, which
converges to a value of roughly -78.888(1). This system is sufficiently small to perform
FCI on. How does ph-AFQMC compare? Why are the error bars getting bigger with back
propagation time?

Finally, we should mention that the path restoration algorithm introduced in M&Z is also
implemented and can be turned on using the `path_restoration` parameter in the Estimator
block.

In QMCPACK path restoration restores both the cosine projection and phase along the back
propagation path. In general it was found in M&Z that path restoration always produced
better results than using the standard back propagation algorithm, and it is recommended
that it is always used. Does path restoration affect the results for methane?
