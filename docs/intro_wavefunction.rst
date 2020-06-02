.. _intro_wavefunction:

Trial wavefunction specificaion
===============================

.. _trial-intro:

Introduction
------------

This section describes the input blocks associated with the specification of the trial wavefunction in a QMCPACK calculation. These sections are contained within the ``<wavefunction>`` :math:`...`  ``</wavefunction>`` xml blocks. **Users are expected to rely on converters to generate the input blocks described in this section.** The converters and the workflows are designed such that input blocks require minimum modifications from users. Unless the workflow requires modification of wavefunction blocks (e.g., setting the cutoff in a multideterminant calculation), only expert users should directly alter them.

The trial wavefunction in QMCPACK has a general product form:

.. math::
  :label: eq1

  \Psi_T(\vec{r}) = \prod_k \Theta_k(\vec{r}) ,

where each :math:`\Theta_k(\vec{r})` is a function of the electron coordinates
(and possibly ionic coordinates and variational parameters).
For problems involving electrons, the overall trial wavefunction
must be antisymmetric with respect to electron exchange,
so at least one of the functions in the product must be
antisymmetric. Notice that, although QMCPACK allows for the
construction of arbitrary trial wavefunctions based on the
functions implemented in the code
(e.g., slater determinants, jastrow functions),
the user must make sure that a correct wavefunction is
used for the problem at hand. From here on, we assume a
standard trial wavefunction for an electronic structure problem

.. math::
  :label: eq2

  Psi_T(\vec{r}) =  \textit{A}(\vec{r}) \prod_k \textit{J}_k(\vec{r}),

where :math:`\textit{A}(\vec{r})`
is one of the antisymmetric functions: (1) slater determinant, (2) multislater determinant, or (3) pfaffian and :math:`\textit{J}_k`
is any of the Jastrow functions (described in :ref:`jastrow`).  The antisymmetric functions are built from a set of single particle orbitals ``(sposet)``. QMCPACK implements four different types of ``sposet``, described in the following section. Each ``sposet`` is designed for a different type of calculation, so their definition and generation varies accordingly.

.. _singledeterminant:

Single determinant wavefunctons
-------------------------------

Placing a single determinant for each spin is the most used ansatz for the antisymmetric part of a trial wavefunction.
The input xml block for ``slaterdeterminant`` is given in :ref:`Listing 1 <Listing 1>`. A list of options is given in
:numref:`Table2`.

.. centered:: ``slaterdeterminant`` element:


.. _Table2:
.. table::

     +-----------------+--------------------+
     | Parent elements | ``determinantset`` |
     +-----------------+--------------------+
     | Child elements  | ``determinant``    |
     +-----------------+--------------------+

.. centered:: Attribute:

+-----------------+----------+--------+---------+------------------------------+
| Name            | Datatype | Values | Default | Description                  |
+=================+==========+========+=========+==============================+
| ``delay_rank``  | Integer  | >=0    | 1       | Number of delayed updates.   |
+-----------------+----------+--------+---------+------------------------------+
| ``optimize``    | Text     | yes/no | yes     | Enable orbital optimization. |
+-----------------+----------+--------+---------+------------------------------+


.. centered:: Table 2 Options for the ``slaterdeterminant`` xml-block.

.. code-block::
      :caption: Slaterdeterminant set XML element.
      :name: Listing 1

      <slaterdeterminant delay_rank="32">
         <determinant id="updet" size="208">
           <occupation mode="ground" spindataset="0">
           </occupation>
         </determinant>
         <determinant id="downdet" size="208">
           <occupation mode="ground" spindataset="0">
           </occupation>
         </determinant>
       </slaterdeterminant>


Additional information:

- ``delay_rank`` This option enables delayed updates of the Slater matrix inverse when particle-by-particle move is used.
  By default or if ``delay_rank=0`` given in the input file, QMCPACK sets 1 for Slater matrices with a leading dimension :math:`<192` and 32 otherwise.
  ``delay_rank=1`` uses the Fahy's variant :cite:`Fahy1990` of the Sherman-Morrison rank-1 update, which is mostly using memory bandwidth-bound BLAS-2 calls.
  With ``delay_rank>1``, the delayed update algorithm :cite:`Luo2018delayedupdate,McDaniel2017` turns most of the computation to compute bound BLAS-3 calls.
  Tuning this parameter is highly recommended to gain the best performance on medium-to-large problem sizes (:math:`>200` electrons).
  We have seen up to an order of magnitude speedup on large problem sizes.
  When studying the performance of QMCPACK, a scan of this parameter is required and we recommend starting from 32.
  The best ``delay_rank`` giving the maximal speedup depends on the problem size.
  Usually the larger ``delay_rank`` corresponds to a larger problem size.
  On CPUs, ``delay_rank`` must be chosen as a multiple of SIMD vector length for good performance of BLAS libraries.
  The best ``delay_rank`` depends on the processor microarchitecture.
  GPU support is under development.

.. _singleparticle:

Single-particle orbitals
------------------------

.. _spo-spline:

Spline basis sets
~~~~~~~~~~~~~~~~~

In this section we describe the use of spline basis sets to expand the ``sposet``.
Spline basis sets are designed to work seamlessly with plane wave DFT code (e.g.,\ Quantum ESPRESSO as a trial wavefunction generator).

.. bibliography:: bibliography.bib
