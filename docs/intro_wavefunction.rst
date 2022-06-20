.. _intro_wavefunction:

Trial wavefunction specification
================================

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
is any of the Jastrow functions (described in :ref:`jastrow`).  The antisymmetric functions are built from a set of single particle orbitals (SPO) ``(sposet)``. QMCPACK implements four different types of ``sposet``, described in the following section. Each ``sposet`` is designed for a different type of calculation, so their definition and generation varies accordingly.

.. code-block::
   :caption: wavefunction XML element skeleton.
   :name: wavefunction.skeleton.xml

   <wavefunction>
     <sposet_collection ...>
       <sposet ...>
         ...
       </sposet>
     </sposet_collection>
     <determinantset>
       <slaterdeterminant ...>
         ...
       </slaterdeterminant>
       <backflow>
         ...
       </backflow>
     </determinantset>
     <jastrow ...>
     </jastrow>
   </wavefunction>

.. _singleparticle:

Single-particle orbitals
------------------------

A single particle orbital set (SPOSet) is a set of orbitals evaluated at a single electron real-space position.
A typical Slater determinant is calculated from a N-by-N matrix constructed from N orbitals at the positions of N electrons.
QMCPACK supports a range of SPOSet types:
 * :ref:`spo-spline`
 * :ref:`spo-lcao`
 * :ref:`spo-hybrid`
 * :ref:`pwbasis`


sposet_collection input style
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block::
  :caption: SPO XML element framework.
  :name: spo.collection.xml

  <!-- build a sposet collection of type bspline. /-->
  <sposet_collection type="bspline" ...>
    <sposet name="spo-up" ... /sposet>
    ...
  </sposet_collection>

The ``sposet_collection`` element forms the container for ``sposet`` and a few other tags.
The contents and attributes in a ``sposet_collection`` node and ``sposet`` node depend on the ``type`` being used.
The ``name`` of each ``sposet`` must be unique. It is used for look-up by :ref:`singledeterminant` and :ref:`multideterminants`.

``sposet_collection`` element:

.. _table1:
.. table::

  +-----------------+------------------+
  | Parent elements | ``wavefunction`` |
  +-----------------+------------------+
  | Child elements  | ``sposet``       |
  +-----------------+------------------+

attribute:

+--------------------+--------------+---------------+-------------+------------------------------------------------+
| **Name**           | **Datatype** | **Values**    | **Default** | **Description**                                |
+====================+==============+===============+=============+================================================+
| ``type``           | Text         | See below     | '' ''       | Type of ``sposet``                             |
+--------------------+--------------+---------------+-------------+------------------------------------------------+

``type`` Type of ``sposet``. Accepted values are 'spline' ('bspline' or 'einspline'), 'MolecularOrbital', 'pw', 'heg', 'composite'.

If QMCPACK printout contains `!!!!!!! Deprecated input style: creating SPO set
inside determinantset. Support for this usage will soon be removed. SPO sets
should be built outside.`, users need to update the input XML by moving all the
SPOSet construction related details out of ``determinantset``. This revised
specification keeps the basis set details separate from information about the
determinants. 

.. code-block::
  :caption: Deprecated input style.
  :name: spo.singledet.old.xml

  <determinantset type="einspline" href="pwscf.pwscf.h5" tilematrix="2 0 0 0 1 0 0 0 1" source="ion0" meshfactor="1.0" precision="double">
     <slaterdeterminant>
        <determinant id="updet" size="8">
           <occupation mode="ground" spindataset="0"/>
        </determinant>
        <determinant id="downdet" size="8">
           <occupation mode="ground" spindataset="0"/>
        </determinant>
     </slaterdeterminant>
  </determinantset>

After updating the input style.

.. code-block::
  :caption: Updated input style.
  :name: spo.singledet.xml

  <!-- all the attributes are moved from determinantset.-->
  <sposet_collection type="einspline" href="pwscf.pwscf.h5" tilematrix="2 0 0 0 1 0 0 0 1" source="ion0" meshfactor="1.0" precision="double">
    <!-- all the attributes and contents are moved from determinant.  Change 'id' tag to 'name' tag.
         Need only one sposet for unpolarized calculation.-->
    <sposet name="spo-ud" size="8">
       <occupation mode="ground" spindataset="0"/>
    </sposet>
  </sposet_collection>
  <determinantset>
     <slaterdeterminant>
        <!-- build two determinants from the same sposet named 'spo-ud'. One for each spin.-->
        <determinant sposet="spo-ud"/>
        <determinant sposet="spo-ud"/>
     </slaterdeterminant>
  </determinantset>


In the case of multi-determinants, all the attributes of ``determinantset`` need to be moved to ``sposet_collection``
and existing ``sposet`` xml nodes need to be moved under ``sposet_collection``. If there is a ``basisset`` node,
it needs to be moved under ``sposet_collection`` as well.

.. _spo-spline:

3D B-splines orbitals
~~~~~~~~~~~~~~~~~~~~~

In this section we describe the use of spline basis sets to expand the
``sposet``. Spline basis sets are designed to work seamlessly with plane wave
DFT codes (e.g.,\ Quantum ESPRESSO as a trial wavefunction generator). Codes
that utilize regular real space grids as a basis can also be seamlessly
interfaced.

In QMC algorithms, all the SPOs :math:`\{\phi(\vec{r})\}` need to be updated
every time a single electron moves. Evaluating SPOs takes a very large portion of computation time.
In principle, PW basis set can be used to express SPOs directly in QMC, as in DFT.
But it introduces an unfavorable scaling because the basis set size increases linearly as the system size.
For this reason, it is efficient to use a localized basis with compact
support and a good transferability from the plane wave basis.

In particular, 3D tricubic B-splines provide a basis in which only
64 elements are nonzero at any given point in :cite:`blips4QMC`.
The 1D cubic B-spline is given by

.. math::
  :label: eq3

  f(x) = \sum_{i'=i-1}^{i+2} b^{i'\!,3}(x)\,\,  p_{i'},

where :math:`b^{i}(x)` is the piecewise cubic polynomial basis functions
and :math:`i = \text{floor}(\Delta^{-1} x)` is the index of the first
grid point :math:`\le x`. Constructing a tensor product in each
Cartesian direction, we can represent a 3D orbital as

.. math::
 :label: eq4

 \phi_n(x,y,z) =
     \!\!\!\!\sum_{i'=i-1}^{i+2} \!\! b_x^{i'\!,3}(x)
     \!\!\!\!\sum_{j'=j-1}^{j+2} \!\! b_y^{j'\!,3}(y)
     \!\!\!\!\sum_{k'=k-1}^{k+2} \!\! b_z^{k'\!,3}(z) \,\, p_{i', j', k',n}.

This allows the rapid evaluation of each orbital in constant time unlike with a plane wave basis set where the cost increases with system size.
Furthermore, this basis is systematically improvable with a single spacing
parameter so that accuracy is not compromised compared with the plane wave basis.

The use of 3D tricubic B-splines greatly improves computational efficiency. The
gain in computation time compared to an equivalent plane wave basis set becomes
increasingly large as the system size grows. On the downside, this computational
efficiency comes at the expense of increased memory use, which is easily
overcome, however, by the large aggregate memory available per node through
OpenMP/MPI hybrid QMC.

The input xml block for the spline SPOs is given in :ref:`spline.spo.xml`. A list of options is given in
:numref:`table3`.

.. code-block::
  :caption: Spline SPO XML element
  :name: spline.spo.xml

  <sposet_collection type="bspline" source="i" href="pwscf.h5"
                  tilematrix="1 1 3 1 2 -1 -2 1 0" gpu="yes" meshfactor="0.8"
                  twist="0  0  0" precision="double">
    <sposet name="spo-up" size="208">
      <occupation mode="ground" spindataset="0"/>
    </sposet>
    <!-- spin polarized case needs two sposets /-->
    <sposet name="spo-dn" size="208">
      <occupation mode="ground" spindataset="1"/>
    </sposet>
  </sposet_collection>


``sposet_collection`` element:

.. _table3:
.. table::

  +-----------------+------------------+
  | Parent elements | ``wavefunction`` |
  +-----------------+------------------+
  | Child elements  | ``sposet``       |
  +-----------------+------------------+

attribute:

+-----------------------------+------------+--------------------------+---------+-------------------------------------------+
| Name                        | Datatype   | Values                   | Default | Description                               |
+=============================+============+==========================+=========+===========================================+
| ``type``                    | Text       | Bspline                  |         | Type of ``sposet``                        |
+-----------------------------+------------+--------------------------+---------+-------------------------------------------+
| ``href``                    | Text       |                          |         | Path to hdf5 file from pw2qmcpack.x.      |
+-----------------------------+------------+--------------------------+---------+-------------------------------------------+
| ``tilematrix``              | 9 integers |                          |         | Tiling matrix used to expand supercell.   |
+-----------------------------+------------+--------------------------+---------+-------------------------------------------+
| ``twistnum``                | Integer    |                          |         | Index of the super twist.                 |
+-----------------------------+------------+--------------------------+---------+-------------------------------------------+
| ``twist``                   | 3 floats   |                          |         | Super twist.                              |
+-----------------------------+------------+--------------------------+---------+-------------------------------------------+
| ``meshfactor``              | Float      | :math:`\le 1.0`          |         | Grid spacing ratio.                       |
+-----------------------------+------------+--------------------------+---------+-------------------------------------------+
| ``precision``               | Text       | Single/double            |         | Precision of spline coefficients          |
+-----------------------------+------------+--------------------------+---------+-------------------------------------------+
| ``gpu``                     | Text       | Yes/no                   |         | GPU switch.                               |
+-----------------------------+------------+--------------------------+---------+-------------------------------------------+
| ``gpusharing``              | Text       | Yes/no                   | No      | Share B-spline table across GPUs.         |
+-----------------------------+------------+--------------------------+---------+-------------------------------------------+
| ``Spline_Size_Limit_MB``    | Integer    |                          |         | Limit B-spline table size on GPU.         |
+-----------------------------+------------+--------------------------+---------+-------------------------------------------+
| ``check_orb_norm``          | Text       | Yes/no                   | Yes     | Check norms of orbitals from h5 file.     |
+-----------------------------+------------+--------------------------+---------+-------------------------------------------+
| ``save_coefs``              | Text       | Yes/no                   | No      | Save the spline coefficients to h5 file.  |
+-----------------------------+------------+--------------------------+---------+-------------------------------------------+
| ``source``                  | Text       | Any                      | Ion0    | Particle set with atomic positions.       |
+-----------------------------+------------+--------------------------+---------+-------------------------------------------+
| ``skip_checks``             | Text       | Yes/no                   | No      | skips checks for ion information in h5    |
+-----------------------------+------------+--------------------------+---------+-------------------------------------------+

.. centered:: Table 3 Options for the ``sposet_collection`` xml-block associated with B-spline single particle orbital sets.

Additional information:

- precision
    Only effective on CPU versions without mixed
    precision, “single" is always imposed with mixed precision. Using
    single precision not only saves memory use but also speeds up the
    B-spline evaluation. We recommend using single precision since we saw
    little chance of really compromising the accuracy of calculation.

- meshfactor
    The ratio of actual grid spacing of B-splines used in
    QMC calculation with respect to the original one calculated from h5.
    A smaller meshfactor saves memory use but reduces accuracy. The
    effects are similar to reducing plane wave cutoff in DFT
    calculations. Use with caution!

- twistnum
    We recommend not using it in the input because the ordering of orbitals
    depends on how they are being stored in the h5 file. ``twistnum`` gets
    ignored if ``twist`` exists in the input. If positive, it is the index.
    If negative, the super twist is referred by ``twist``. This input
    parameter is kept only for keeping old input files working.

- twist
    The twist angle. If neither ``twist`` nor ``twistnum`` is provided,
    Take Gamma point, (0, 0, 0).

- save_coefs
    If yes, dump the real-space B-spline coefficient
    table into an h5 file on the disk. When the orbital transformation
    from k space to B-spline requires more than the available amount of
    scratch memory on the compute nodes, users can perform this step on
    fat nodes and transfer back the h5 file for QMC calculations.

- gpusharing
    If enabled, spline data is shared across multiple
    GPUs on a given computational node. For example, on a
    two-GPU-per-node system, each GPU would have half of the orbitals.
    This enables larger overall spline tables than would normally fit in
    the memory of individual GPUs to be used, potentially up to the total
    GPU memory on a node. To obtain high performance, large electron
    counts or a high-performing CPU-GPU interconnect is required. To use
    this feature, the following needs to be done:

      -  The CUDA Multi-Process Service (MPS) needs to be used (e.g., on
         Summit use "-alloc_flags gpumps" for bsub). If MPS is not
         detected, sharing will be disabled.

      -  CUDA_VISIBLE_DEVICES needs to be properly set to control each rank’s
         visible CUDA devices (e.g., on OLCF Summit one needs to
         create a resource set containing all GPUs with the respective number
         of ranks with "jsrun –task-per-rs Ngpus -g Ngpus").

- Spline_Size_Limit_MB
    Allows distribution of the B-spline
    coefficient table between the host and GPU memory. The compute kernels
    access host memory via zero-copy. Although the performance penalty
    introduced by it is significant, it allows large calculations to go
    through.
 
- skip_checks
    When converting the wave function from convertpw4qmc instead
    of pw2qmcpack, there is missing ionic information. This flag bypasses the requirement
    that the ionic information in the eshdf.h5 file match the input xml. 

.. _spo-lcao:

Linear combination of atomic orbitals (LCAO) with Gaussian and/or Slater-type basis sets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this section we describe the use of localized basis sets to expand the ``sposet``. The general form of a single particle orbital in this case is given by:

.. math::
  :label: eq5

  \phi_i(\vec{r}) = \sum_k C_{i,k} \ \eta_k(\vec{r}),

where :math:`\{\eta_k(\vec{r})\}` is a set of M atom-centered basis
functions and :math:`C_{i,k}` is a coefficient matrix. This should be
used in calculations of finite systems employing an atom-centered basis
set and is typically generated by the *convert4qmc* converter. Examples
include calculations of molecules using Gaussian basis sets or
Slater-type basis functions. Initial support for periodic systems is
described in :ref:`LCAO`. Even though this section is called
"Gaussian basis sets" (by far the most common atom-centered basis set),
QMCPACK works with any atom-centered basis set based on either spherical
harmonic angular functions or Cartesian angular expansions. The radial
functions in the basis set can be expanded in either Gaussian functions,
Slater-type functions, or numerical radial functions.

In this section we describe the input sections of ``sposet_collection`` for the atom-centered basis set.
Here is an :ref:`example <spo.singledet.lcao.xml>` of single determinant with LCAO.
The input sections for multideterminant trial wavefunctions are described in :ref:`multideterminants`.

.. code-block::
  :caption: ``slaterdeterminant`` with an LCAO ``sposet_collection`` example
  :name: spo.singledet.lcao.xml

  <sposet_collection type="MolecularOrbital" source="ion0" cuspCorrection="no">
    <basisset name="LCAOBSet">
      <atomicBasisSet name="Gaussian-G2" angular="cartesian" elementType="H" normalized="no">
        <grid type="log" ri="1.e-6" rf="1.e2" npts="1001"/>
        <basisGroup rid="H00" n="0" l="0" type="Gaussian">
          <radfunc exponent="5.134400000000e-02" contraction="1.399098787100e-02"/>
        </basisGroup>
      </atomicBasisSet>
    </basisset>
    <sposet name="spo" basisset="LCAOBSet" size="1">
      <occupation mode="ground"/>
      <coefficient size="1" id="updetC">
        1.00000000000000e+00
      </coefficient>
    </sposet>
  </sposet_collection>
  <determinantset>
     <slaterdeterminant>
        <determinant sposet="spo" />
     </slaterdeterminant>
  </determinantset>

Here is the :ref:`basic structure <spo.lcao.xml>` for LCAO ``sposet_collection`` input block.
A list of options for ``sposet_collection`` is given in :numref:`table4`.

.. code-block::
   :caption: Basic input block for ``sposet_collection`` for LCAO.
   :name: spo.lcao.xml

   <sposet_collection type="MolecularOrbital" ...>
     <basisset name="LCAOBSet" ...>
       ...
     </basisset>
     <sposet name="spo" basisset="LCAOBSet" size="1">
       <occupation mode="ground"/>
       <coefficient size="1" id="updetC">
         1.00000000000000e+00
       </coefficient>
     </sposet>
   </sposet_collection>


The definition of the set of atom-centered basis functions is given by the ``basisset`` block and the ``sposet`` defined within ``sposet_collection``.
The ``basisset`` input block is composed from a collection of ``atomicBasisSet`` input blocks, one for each atomic species in the simulation where basis functions are centered.
The general structure for ``basisset`` and ``atomicBasisSet`` are given in :ref:`Listing 4 <Listing 4>`, and the corresponding lists of options are given in
:numref:`table5` and :numref:`table6`.

``sposet_collection`` element:

.. _table4:
.. table::

  +-----------------+---------------------------+
  | Parent elements | ``wavefunction``          |
  +-----------------+---------------------------+
  | Child elements  | ``basisset`` , ``sposet`` |
  +-----------------+---------------------------+

Attribute:

+--------------------+--------------+---------------+-------------+------------------------------------------------+
| **Name**           | **Datatype** | **Values**    | **Default** | **Description**                                |
+====================+==============+===============+=============+================================================+
| ``name/id``        | Text         | *Any*         | '' ''       | Name of determinant set                        |
+--------------------+--------------+---------------+-------------+------------------------------------------------+
| ``type``           | Text         | See below     | '' ''       | Type of ``sposet``                             |
+--------------------+--------------+---------------+-------------+------------------------------------------------+
| ``keyword``        | Text         | NMO, GTO, STO | NMO         | Type of orbital set generated                  |
+--------------------+--------------+---------------+-------------+------------------------------------------------+
| ``transform``      | Text         | Yes/no        | Yes         | Transform to numerical radial functions?       |
+--------------------+--------------+---------------+-------------+------------------------------------------------+
| ``source``         | Text         | *Any*         | Ion0        | Particle set with the position of atom centers |
+--------------------+--------------+---------------+-------------+------------------------------------------------+
| ``cuspCorrection`` | Text         | Yes/no        | No          | Apply cusp correction scheme to ``sposet``?    |
+--------------------+--------------+---------------+-------------+------------------------------------------------+

.. centered:: Table 4 Options for the ``sposet_collection`` xml-block associated with atom-centered single particle orbital sets.


- type
    Type of ``sposet``. For atom-centered based ``sposets``, use type="MolecularOrbital" or type="MO".

- keyword/key
    Type of basis set generated, which does not necessarily match the type of basis set on the input block. The three possible options are: NMO (numerical molecular orbitals), GTO (Gaussian-type orbitals), and STO (Slater-type orbitals). The default option is NMO. By default, QMCPACK will generate numerical orbitals from both GTO and STO types and use cubic or quintic spline interpolation to evaluate the radial functions. This is typically more efficient than evaluating the radial functions in the native basis (Gaussians or exponents) and allows for arbitrarily large contractions without any additional cost. To force use of the native expansion (not recommended), use GTO or STO for each type of input basis set.

- transform
    Request (or avoid) a transformation of the radial functions to NMO type. The default and recommended behavior is to transform to numerical radial functions. If ``transform`` is set to *yes*, the option ``keyword`` is ignored.

- cuspCorrection
    Enable (disable) use of the cusp correction algorithm (CASINO REFERENCE) for a ``basisset`` built with GTO functions. The algorithm is implemented as described in (CASINO REFERENCE) and works only with transform="yes" and an input GTO basis set. No further input is needed.

.. code-block::
  :caption: Basic input block for ``basisset``.
  :name: Listing 4

  <basisset name="LCAOBSet">
    <atomicBasisSet name="Gaussian-G2" angular="cartesian" elementType="C" normalized="no">
      <grid type="log" ri="1.e-6" rf="1.e2" npts="1001"/>
      <basisGroup rid="C00" n="0" l="0" type="Gaussian">
        <radfunc exponent="5.134400000000e-02" contraction="1.399098787100e-02"/>
        ...
      </basisGroup>
      ...
    </atomicBasisSet>
    <atomicBasisSet name="Gaussian-G2" angular="cartesian" type="Gaussian" elementType="C" normalized="no">
      ...
    </atomicBasisSet>
    ...
  </basisset>

``basisset`` element:

.. _table5:
.. table::

  +-----------------+----------------------+
  | Parent elements | ``sposet_collection``|
  +-----------------+----------------------+
  | Child elements  | ``atomicBasisSet``   |
  +-----------------+----------------------+

Attribute:

+-------------------+--------------+------------+-------------+----------------------------------+
| **Name**          | **Datatype** | **Values** | **Default** | **Description**                  |
+===================+==============+============+=============+==================================+
| ``name`` / ``id`` | Text         | *Any*      | " "         | Name of atom-centered basis set  |
+-------------------+--------------+------------+-------------+----------------------------------+

.. centered:: Table 5 Options for the ``basisset`` xml-block associated with atom-centered single particle orbital sets.

``AtomicBasisSet`` element:

.. _table6:
.. table::

  +-----------------+--------------------------+
  | Parent elements | ``basisset``             |
  +-----------------+--------------------------+
  | Child elements  | ``grid`` , ``basisGroup``|
  +-----------------+--------------------------+

Attribute:

+-------------------------+--------------+------------+-------------+---------------------------------------------+
| **Name**                | **Datatype** | **Values** | **Default** | **Description**                             |
+=========================+==============+============+=============+=============================================+
| ``name`` / ``id``       | Text         | *Any*      | " "         | Name of atomic basis set                    |
+-------------------------+--------------+------------+-------------+---------------------------------------------+
| ``angular``             | Text         | See below  | Default     | Type of angular functions                   |
+-------------------------+--------------+------------+-------------+---------------------------------------------+
| ``expandYlm``           | Text         | See below  | Yes         | Expand Ylm shells?                          |
+-------------------------+--------------+------------+-------------+---------------------------------------------+
| ``expM``                | Text         | See below  | Yes         | Add sign for :math:`(-1)^{m}`?              |
+-------------------------+--------------+------------+-------------+---------------------------------------------+
| ``elementType/species`` | Text         | *Any*      | e           | Atomic species where functions are centered |
+-------------------------+--------------+------------+-------------+---------------------------------------------+
| ``normalized``          | Text         | Yes/no     | Yes         | Are single particle functions normalized?   |
+-------------------------+--------------+------------+-------------+---------------------------------------------+

.. centered:: Table 6 Options for the ``atomicBasisSet`` xml-block.

- name/id
    Name of the basis set. Names should be unique.

- angular
    Type of angular functions used in the expansion. In general, two angular basis functions are allowed: "spherical" (for spherical Ylm functions) and "Cartesian" (for functions of the type :math:`x^{n}y^{m}z^{l}`).

- expandYlm
    Determines whether each basis group is expanded across the corresponding shell of m values (for spherical type) or consistent powers (for Cartesian functions). Options:

      - "No": Do not expand angular functions across corresponding angular shell.

      - "Gaussian": Expand according to Gaussian03 format. This function is compatible only with angular="spherical." For a given input (l,m), the resulting order of the angular functions becomes (1,-1,0) for l=1 and (0,1,-1,2,-2,...,l,-l) for general l.

      - "Natural": Expand angular functions according to (-l,-l+1,...,l-1,l).

      - "Gamess": Expand according to Gamess' format for Cartesian functions. Notice that this option is compatible only with angular="Cartesian." If angular="Cartesian" is used, this option is not necessary.

- expM
    Determines whether the sign of the spherical Ylm function associated with m (:math:`-1^{m}`) is included in the coefficient matrix or not.

- elementType/species
    Name of the species where basis functions are centered. Only one ``atomicBasisSet`` block is allowed per species. Additional blocks are ignored. The corresponding species must exist in the ``particleset`` given as the ``source`` option to ``determinantset``. Basis functions for all the atoms of the corresponding species are included in the basis set, based on the order of atoms in the ``particleset``.

``basicGroup`` element:

.. _table7:
.. table::

  +-----------------+-------------------+
  | Parent elements | ``AtomicBasisSet``|
  +-----------------+-------------------+
  | Child elements  | ``radfunc``       |
  +-----------------+-------------------+

Attribute:

+-------------+--------------+------------+-------------+-------------------------------+
| **Name**    | **Datatype** | **Values** | **Default** | **Description**               |
+=============+==============+============+=============+===============================+
| ``rid/id``  | Text         | *Any*      | '' ''       | Name of the basisGroup        |
+-------------+--------------+------------+-------------+-------------------------------+
| ``type``    | Text         | *Any*      | '' ''       | Type of basisGroup            |
+-------------+--------------+------------+-------------+-------------------------------+
| ``n/l/m/s`` | Integer      | *Any*      | 0           | Quantum numbers of basisGroup |
+-------------+--------------+------------+-------------+-------------------------------+

.. centered:: :numref:`table7` Options for the ``basisGroup`` xml-block.

- type
    Type of input basis radial function. Note that this refers to the type of radial function in the input xml-block, which might not match the radial function generated internally and used in the calculation (if ``transform`` is set to "yes"). Also note that different ``basisGroup`` blocks within a given ``atomicBasisSet`` can have different ``types``.

- n/l/m/s
    Quantum numbers of the basis function. Note that if
    ``expandYlm`` is set to *"yes"* in ``atomicBasisSet``, a
    full shell of basis functions with the appropriate values of
    *"m"* will be defined for the corresponding value of
    *"l."* Otherwise a single basis function will be given for the
    specific combination of *"(l,m)."*


``radfunc`` element:
  attributes for ``type`` = *"Gaussian"*:
``TBDoc``

.. _spo-hybrid:

Hybrid orbital representation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The hybrid representation of the single particle orbitals combines a localized atomic basis set around atomic cores and B-splines in the interstitial regions to reduce memory use while retaining high evaluation speed and either retaining or increasing overall accuracy. Full details are provided in :cite:`Luo2018hyb`, and **users of this feature are kindly requested to cite this paper**.
In practice, we have seen that using a meshfactor=0.5 is often possible and achieves huge memory savings.
:numref:`fig3` illustrates how the regions are assigned.

.. _fig3:
.. figure:: /figs/hybrid_new.jpg
    :width: 400
    :align: center

    Regular and hybrid orbital representation. Regular B-spline representation (left panel) contains only one region and a sufficiently fine mesh to resolve orbitals near the nucleus. The hybrid orbital representation (right panel) contains near nucleus regions (A) where spherical harmonics and radial functions are used, buffers or interpolation regions (B), and an interstitial region (C) where a coarse B-spline mesh is used.

Orbitals within region A are computed as

.. math:: \phi^A_n({\bf r})=R_{n,l,m}(r)Y_{l,m}(\hat{r})

Orbitals in region C are computed as the regular B-spline basis described in :ref:`spo-spline` above. The region B interpolates between A and C as

.. math::
  :label: eq6

  \phi^B_n({\bf r}) = S(r) \phi^A_n({\bf r}) + (1-S(r))\phi^C_n({\bf r})

.. math::
  :label: eq7

  (S(r) = \frac{1}{2}-\frac{1}{2} tanh \left[\alpha\left(\frac{r-r_{\rm A/B}}{r_{\rm B/C}-r_{\rm A/B}}-\frac{1}{2}\right)\right]

To enable hybrid orbital representation, the input XML needs to see the tag ``hybridrep="yes"`` shown in :ref:`Listing 6 <Listing 6>`.

.. code-block::
  :caption: Hybrid orbital representation input example.
  :name: Listing 6

  <sposet_collection type="bspline" source="i" href="pwscf.h5"
                tilematrix="1 1 3 1 2 -1 -2 1 0" gpu="yes" meshfactor="0.8"
                twist="0  0  0" precision="single" hybridrep="yes">
    ...
  </sposet_collection>

Second, the information describing the atomic regions is required in the particle set, shown in :ref:`Listing 7 <Listing 7>`.

.. code-block::
  :caption: particleset elements for ions with information needed by hybrid orbital representation.
  :name: Listing 7

  <group name="Ni">
    <parameter name="charge">          18 </parameter>
    <parameter name="valence">         18 </parameter>
    <parameter name="atomicnumber" >   28 </parameter>
    <parameter name="cutoff_radius" > 1.6 </parameter>
    <parameter name="inner_cutoff" >  1.3 </parameter>
    <parameter name="lmax" >            5 </parameter>
    <parameter name="spline_radius" > 1.8 </parameter>
    <parameter name="spline_npoints">  91 </parameter>
  </group>

The parameters specific to hybrid representation are listed as

``attrib`` element

Attribute:

+---------------------+--------------+------------+-------------+---------------------------------------+
| **Name**            | **Datatype** | **Values** | **Default** | **Description**                       |
+=====================+==============+============+=============+=======================================+
| ``cutoff_radius``   | Real         | >=0.0      | *None*      | Cutoff radius for B/C boundary        |
+---------------------+--------------+------------+-------------+---------------------------------------+
| ``lmax``            | Integer      | >=0        | *None*      | Largest angular channel               |
+---------------------+--------------+------------+-------------+---------------------------------------+
| ``inner_cutoff``    | Real         | >=0.0      | Dep.        | Cutoff radius for A/B boundary        |
+---------------------+--------------+------------+-------------+---------------------------------------+
| ``spline_radius``   | Real         | >0.0       | Dep.        | Radial function radius used in spine  |
+---------------------+--------------+------------+-------------+---------------------------------------+
| ``spline_npoints``  | Integer      | >0         | Dep.        | Number of spline knots                |
+---------------------+--------------+------------+-------------+---------------------------------------+

- ``cutoff_radius``  is required for every species. If a species is intended to not be covered by atomic regions, setting the value 0.0 will put default values for all the reset parameters. A good value is usually a bit larger than the core radius listed in the pseudopotential file. After a parametric scan, pick the one from the flat energy region with the smallest variance.

- ``lmax`` is required if ``cutoff_radius`` :math:`>` 0.0. This value usually needs to be at least the highest angular momentum plus 2.

- ``inner_cutoff`` is optional and set as ``cutoff_radius`` :math:`-0.3` by default, which is fine in most cases.

- ``spline_radius`` and ``spline_npoints`` are optional. By default, they are calculated based on ``cutoff_radius`` and a grid displacement 0.02 bohr.
  If users prefer inputing them, it is required that ``cutoff_radius`` <=  ``spline_radius`` :math:`-` 2 :math:`\times` ``spline_radius``/(``spline_npoints`` :math:`-` 1).

In addition, the hybrid orbital representation allows extra optimization to speed up the nonlocal pseudopotential evaluation using the batched algorithm listed in :ref:`nlpp`.

.. _pwbasis:

Plane-wave basis sets
~~~~~~~~~~~~~~~~~~~~~

.. _hegbasis:

Homogeneous electron gas
~~~~~~~~~~~~~~~~~~~~~~~~

The interacting Fermi liquid has its own special ``determinantset`` for filling up a
Fermi surface.  The shell number can be specified separately for both spin-up and spin-down.
This determines how many electrons to include of each time; only closed shells are currently
implemented.  The shells are filled according to the rules of a square box; if other lattice
vectors are used, the electrons might not fill up a complete shell.

This following example can also be used for Helium simulations by specifying the
proper pair interaction in the Hamiltonian section.

.. code-block::
  :caption: 2D Fermi liquid example: particle specification
  :name: Listing 8

  <simulationcell name="global">
    <parameter name="rs" pol="0" condition="74">6.5</parameter>
    <parameter name="bconds">p p p</parameter>
    <parameter name="LR_dim_cutoff">15</parameter>
  </simulationcell>
  <particleset name="e" random="yes">
    <group name="u" size="37">
      <parameter name="charge">-1</parameter>
      <parameter name="mass">1</parameter>
    </group>
    <group name="d" size="37">
      <parameter name="charge">-1</parameter>
      <parameter name="mass">1</parameter>
    </group>
  </particleset>

.. code-block::
  :caption: 2D Fermi liquid example (Slater Jastrow wavefunction)
  :name: Listing 9

  <wavefunction name="psi0" target="e">
    <determinantset type="electron-gas" shell="7" shell2="7" randomize="true">
  </determinantset>
  <jastrow name="J2" type="Two-Body" function="Bspline" print="no">
    <correlation speciesA="u" speciesB="u" size="8" cusp="0">
      <coefficients id="uu" type="Array" optimize="yes">
    </correlation>
    <correlation speciesA="u" speciesB="d" size="8" cusp="0">
      <coefficients id="ud" type="Array" optimize="yes">
    </correlation>
  </jastrow>
  </wavefunction>

.. _singledeterminant:

Single determinant wavefunctions
--------------------------------

Placing a single determinant for each spin is the most used ansatz for the antisymmetric part of a trial wavefunction.
The input xml block for ``slaterdeterminant`` is given in :ref:`Listing 1 <Listing 1>`. A list of options is given in
:numref:`Table2`.

``slaterdeterminant`` element:


.. _Table2:
.. table::

     +-----------------+--------------------+
     | Parent elements | ``determinantset`` |
     +-----------------+--------------------+
     | Child elements  | ``determinant``    |
     +-----------------+--------------------+

Attribute:

+-----------------------+----------+----------+---------+-------------------------------------------+
| Name                  | Datatype | Values   | Default | Description                               |
+=======================+==========+==========+=========+===========================================+
| ``delay_rank``        | Integer  | >=0      | 1       | Number of delayed updates.                |
+-----------------------+----------+----------+---------+-------------------------------------------+
| ``optimize``          | Text     | yes/no   | yes     | Enable orbital optimization.              |
+-----------------------+----------+----------+---------+-------------------------------------------+
| ``gpu``               | Text     | yes/no   | yes     | Use the GPU acceleration implementation.  |
+-----------------------+----------+----------+---------+-------------------------------------------+
| ``batch``             | Text     | yes/no   | dep.    | Select the batched walker implementation. |
+-----------------------+----------+----------+---------+-------------------------------------------+
| ``matrix_inverter``   | Text     | gpu/host | gpu     | Slater matrix inversion scheme.           |
+-----------------------+----------+----------+---------+-------------------------------------------+


.. centered:: Table 2 Options for the ``slaterdeterminant`` xml-block.

.. code-block::
   :caption: Slaterdeterminant set XML element.
   :name: Listing 1

   <sposet_collection ...>
     <sposet name="spo" size="8">
       ...
     </sposet>
   </sposet_collection>
   <determinantset>
     <slaterdeterminant delay_rank="32">
       <determinant sposet="spo"/>
       <determinant sposet="spo"/>
     </slaterdeterminant>
   </determinantset>


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

- ``gpu`` This option is only effective when GPU features are built. Use the implementation with GPU acceleration if ``yes``.

- ``batch`` The default value is ``yes`` if ``gpu=yes`` and ``no`` otherwise.

- ``matrix_inverter`` If the value is ``gpu``, the inversion happens on the GPU and additional GPU memory is needed.
  If the value is ``host``, the inversion happens on the CPU and doesn't need GPU memory.

.. _multideterminants:

Multideterminant wavefunctions
------------------------------

``multideterminant`` element:


.. _Table_msd:
.. table::

     +-----------------+--------------------+
     | Parent elements | ``determinantset`` |
     +-----------------+--------------------+
     | Child elements  | ``detlist``        |
     +-----------------+--------------------+

Attribute:

+-----------------------+----------+----------+--------------------------+-------------------------------------------+
| Name                  | Datatype | Values   | Default                  | Description                               |
+=======================+==========+==========+==========================+===========================================+
| ``optimize``          | Text     | yes/no   | yes                      | Enable optimization.                      |
+-----------------------+----------+----------+--------------------------+-------------------------------------------+
| ``spo_up``            | Text     |          |                          | The name of SPO for spin up electrons     |
+-----------------------+----------+----------+--------------------------+-------------------------------------------+
| ``spo_down``          | Text     |          |                          | The name of SPO for spin down electrons   |
+-----------------------+----------+----------+--------------------------+-------------------------------------------+
| ``algorithm``         | Text     |          | precomputed_table_method | Slater matrix inversion scheme.           |
+-----------------------+----------+----------+--------------------------+-------------------------------------------+

.. centered:: Table 3 Options for the ``multideterminant`` xml-block.

Additional information:

- ``algorithm`` algorithms used in multi-Slater determinant implementation. ``table_method`` table method of Clark et al. :cite:`Clark2011` .
  ``precomputed_table_method`` adds partial sum precomputation on top of ``table_method``.

.. code-block::
   :caption: multideterminant set XML element.
   :name: multideterminant.xml

   <sposet_collection ...>
     <sposet name="spo" size="85">
       ...
     </sposet>
   </sposet_collection>
   <determinantset>
     <multideterminant optimize="yes" spo_up="spo" spo_dn="spo">
       <detlist size="1487" type="DETS" nca="0" ncb="0" nea="2" neb="2" nstates="85" cutoff="1e-20" href="LiH.orbs.h5">
     </multideterminant>
   </determinantset>

Multiple schemes to generate a multideterminant wavefunction are
possible, from CASSF to full CI or selected CI. The QMCPACK converter can
convert MCSCF multideterminant wavefunctions from
GAMESS :cite:`schmidt93` and CIPSI :cite:`Caffarel2013` wavefunctions from
Quantum Package :cite:`QP` (QP). Full details of how to run a CIPSI
calculation and convert the wavefunction for QMCPACK are given in
:ref:`cipsi`.

The script ``utils/determinants_tools.py`` can be used to generate
useful information about the multideterminant wavefunction. This script takes, as a required argument, the path of an h5 file corresponding to the wavefunction. Used without optional arguments, it prints the number of determinants, the number of CSFs, and a histogram of the excitation degree.

::

  > determinants_tools.py ./tests/molecules/C2_pp/C2.h5
  Summary:
  excitation degree 0 count: 1
  excitation degree 1 count: 6
  excitation degree 2 count: 148
  excitation degree 3 count: 27
  excitation degree 4 count: 20

  n_det 202
  n_csf 104

If the ``--verbose`` argument is used, the script will print each determinant,
the associated CSF, and the excitation degree relative to the first determinant.

::

  > determinants_tools.py -v ./tests/molecules/C2_pp/C2.h5 | head
  1
  alpha  1111000000000000000000000000000000000000000000000000000000
  beta   1111000000000000000000000000000000000000000000000000000000
  scf    2222000000000000000000000000000000000000000000000000000000
  excitation degree  0

  2
  alpha  1011100000000000000000000000000000000000000000000000000000
  beta   1011100000000000000000000000000000000000000000000000000000
  scf    2022200000000000000000000000000000000000000000000000000000
  excitation degree  2

.. _backflow:

Backflow Wavefunctions
----------------------

One can perturb the nodal surface of a single-Slater/multi-Slater
wavefunction through use of a backflow transformation. Specifically, if
we have an antisymmetric function
:math:`D(\mathbf{x}_{0\uparrow},\cdots,\mathbf{x}_{N\uparrow}, \mathbf{x}_{0\downarrow},\cdots,\mathbf{x}_{N\downarrow})`,
and if :math:`i_\alpha` is the :math:`i`-th particle of species type
:math:`\alpha`, then the backflow transformation works by making the
coordinate transformation
:math:`\mathbf{x}_{i_\alpha} \to \mathbf{x}'_{i_\alpha}` and evaluating
:math:`D` at these new “quasiparticle" coordinates. QMCPACK currently
supports quasiparticle transformations given by

.. math::
  :label: eq24

  \mathbf{x}'_{i_\alpha}=\mathbf{x}_{i_\alpha}+\sum_{\alpha \leq \beta} \sum_{i_\alpha \neq j_\beta} \eta^{\alpha\beta}(|\mathbf{x}_{i_\alpha}-\mathbf{x}_{j_\beta}|)(\mathbf{x}_{i_\alpha}-\mathbf{x}_{j_\beta})\:.

Here, :math:`\eta^{\alpha\beta}(|\mathbf{x}_{i_\alpha}-\mathbf{x}_{j_\beta}|)`
is a radially symmetric backflow transformation between species
:math:`\alpha` and :math:`\beta`. In QMCPACK, particle :math:`i_\alpha`
is known as the “target" particle and :math:`j_\beta` is known as the
“source." The main types of transformations are so-called one-body
terms, which are between an electron and an ion
:math:`\eta^{eI}(|\mathbf{x}_{i_e}-\mathbf{x}_{j_I}|)` and two-body
terms. Two-body terms are distinguished as those between like and
opposite spin electrons:
:math:`\eta^{e(\uparrow)e(\uparrow)}(|\mathbf{x}_{i_e(\uparrow)}-\mathbf{x}_{j_e(\uparrow)}|)`
and
:math:`\eta^{e(\uparrow)e(\downarrow)}(|\mathbf{x}_{i_e(\uparrow)}-\mathbf{x}_{j_e(\downarrow)}|)`.
Henceforth, we will assume that
:math:`\eta^{e(\uparrow)e(\uparrow)}=\eta^{e(\downarrow)e(\downarrow)}`.

In the following, we explain how to describe general terms such as
:eq:`eq24` in a QMCPACK XML file. For specificity, we will
consider a particle set consisting of H and He (in that order). This
ordering will be important when we build the XML file, so you can find
this out either through your specific declaration of <particleset>, by
looking at the hdf5 file in the case of plane waves, or by looking at
the QMCPACK output file in the section labeled “Summary of QMC systems."

Input specifications
~~~~~~~~~~~~~~~~~~~~

All backflow declarations occur within a single ``<backflow> ... </backflow>`` block.  Backflow transformations occur in ``<transformation>`` blocks and have the following input parameters:

Transformation element:

  +----------+--------------+------------+-------------+----------------------------------------------------------+
  | **Name** | **Datatype** | **Values** | **Default** | **Description**                                          |
  +==========+==============+============+=============+==========================================================+
  | name     | Text         |            | (Required)  | Unique name for this Jastrow function.                   |
  +----------+--------------+------------+-------------+----------------------------------------------------------+
  | type     | Text         | "e-I"      | (Required)  | Define a one-body backflow transformation.               |
  +----------+--------------+------------+-------------+----------------------------------------------------------+
  |          | Text         | "e-e"      |             | Define a two-body backflow transformation.               |
  +----------+--------------+------------+-------------+----------------------------------------------------------+
  | function | Text         | B-spline   | (Required)  | B-spline type transformation (no other types supported). |
  +----------+--------------+------------+-------------+----------------------------------------------------------+
  | source   | Text         |            |             | "e" if two body, ion particle set if one body.           |
  +----------+--------------+------------+-------------+----------------------------------------------------------+

Just like one- and two-body jastrows, parameterization of the backflow transformations are specified within the ``<transformation>`` blocks by  ``<correlation>`` blocks.  Please refer to :ref:`onebodyjastrowspline` for more information.

Example Use Case
~~~~~~~~~~~~~~~~

Having specified the general form, we present a general example of one-body and two-body backflow transformations in a hydrogen-helium mixture.  The hydrogen and helium ions have independent backflow transformations, as do the like and unlike-spin two-body terms.  One caveat is in order:  ionic backflow transformations must be listed in the order they appear in the particle set.  If in our example, helium is listed first and hydrogen is listed second, the following example would be correct.  However, switching backflow declaration to hydrogen first then helium, will result in an error.  Outside of this, declaration of one-body blocks and two-body blocks are not sensitive to ordering.

::

  <backflow>
  <!--The One-Body term with independent e-He and e-H terms. IN THAT ORDER -->
  <transformation name="eIonB" type="e-I" function="Bspline" source="ion0">
      <correlation cusp="0.0" size="8" type="shortrange" init="no" elementType="He" rcut="3.0">
          <coefficients id="eHeC" type="Array" optimize="yes">
              0 0 0 0 0 0 0 0
          </coefficients>
      </correlation>
      <correlation cusp="0.0" size="8" type="shortrange" init="no" elementType="H" rcut="3.0">
          <coefficients id="eHC" type="Array" optimize="yes">
              0 0 0 0 0 0 0 0
          </coefficients>
      </correlation>
  </transformation>

  <!--The Two-Body Term with Like and Unlike Spins -->
  <transformation name="eeB" type="e-e" function="Bspline" >
      <correlation cusp="0.0" size="7" type="shortrange" init="no" speciesA="u" speciesB="u" rcut="1.2">
          <coefficients id="uuB1" type="Array" optimize="yes">
              0 0 0 0 0 0 0
          </coefficients>
      </correlation>
      <correlation cusp="0.0" size="7" type="shortrange" init="no" speciesA="d" speciesB="u" rcut="1.2">
          <coefficients id="udB1" type="Array" optimize="yes">
              0 0 0 0 0 0 0
          </coefficients>
      </correlation>
  </transformation>
  </backflow>

Currently, backflow works only with single-Slater determinant wavefunctions.  When a backflow transformation has been declared, it should be placed within the ``<determinantset>`` block, but outside of the ``<slaterdeterminant>`` blocks, like so:

::

  <determinantset ... >
      <!--basis set declarations go here, if there are any -->

      <backflow>
          <transformation ...>
            <!--Here is where one and two-body terms are defined -->
           </transformation>
       </backflow>

       <slaterdeterminant>
           <!--Usual determinant definitions -->
       </slaterdeterminant>
   </determinantset>

Optimization Tips
~~~~~~~~~~~~~~~~~

Backflow is notoriously difficult to optimize---it is extremely nonlinear in the variational parameters and moves the nodal surface around.  As such, it is likely that a full Jastrow+Backflow optimization with all parameters initialized to zero might not converge in a reasonable time.  If you are experiencing this problem, the following pointers are suggested (in no particular order).

Get a good starting guess for :math:`\Psi_T`:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Try optimizing the Jastrow first without backflow.

#. Freeze the Jastrow parameters, introduce only the e-e terms in the
   backflow transformation, and optimize these parameters.

#. Freeze the e-e backflow parameters, and then optimize the e-I terms.

   -  If difficulty is encountered here, try optimizing each species
      independently.


#. Unfreeze all Jastrow, e-e backflow, and e-I backflow parameters, and
   reoptimize.

Optimizing Backflow Terms
^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible that the previous prescription might grind to a halt in steps 2 or 3 with the inability to optimize the e-e or e-I backflow transformation independently, especially if it is initialized to zero.  One way to get around this is to build a good starting guess for the e-e or e-I backflow terms iteratively as follows:

#. Start off with a small number of knots initialized to zero. Set
   :math:`r_{cut}` to be small (much smaller than an interatomic distance).

#. Optimize the backflow function.

#. If this works, slowly increase :math:`r_{cut}` and/or the number of
   knots.

#. Repeat steps 2 and 3 until there is no noticeable change in energy or
   variance of :math:`\Psi_T`.

Tweaking the Optimization Run
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following modifications are worth a try in the optimization block:

-  Try setting “useDrift" to “no." This eliminates the use of
   wavefunction gradients and force biasing in the VMC algorithm. This
   could be an issue for poorly optimized wavefunctions with
   pathological gradients.

-  Try increasing “exp0" in the optimization block. Larger values of
   exp0 cause the search directions to more closely follow those
   predicted by steepest-descent than those by the linear method.

Note that the new adaptive shift optimizer has not yet been tried with
backflow wavefunctions. It should perform better than the older
optimizers, but a considered optimization process is still recommended.

.. _jastrow:

Jastrow Factors
---------------

Jastrow factors are among the simplest and most effective ways of including
dynamical correlation in the trial many body wavefunction.  The resulting many body
wavefunction is expressed as the product of an antisymmetric (in the case
of Fermions) or symmetric (for Bosons) part and a correlating Jastrow factor
like so:

.. math::
  :label: eq8

  \Psi(\vec{R}) = \mathcal{A}(\vec{R}) \exp\left[J(\vec{R})\right]

In this section we will detail the types and forms of Jastrow factor used
in QMCPACK.  Note that each type of Jastrow factor needs to be specified using
its own individual ``jastrow`` XML element.  For this reason, we have repeated the
specification of the ``jastrow`` tag in each section, with specialization for the
options available for that given type of Jastrow.

.. _onebodyjastrow:

One-body Jastrow functions
~~~~~~~~~~~~~~~~~~~~~~~~~~

The one-body Jastrow factor is a form that allows for the direct inclusion
of correlations between particles that are included in the wavefunction with
particles that are not explicitly part of it.  The most common example of
this are correlations between electrons and ions.

The Jastrow function is specified within a ``wavefunction`` element
and must contain one or more ``correlation`` elements specifying
additional parameters as well as the actual coefficients.
:ref:`1bjsplineexamples` gives examples of the typical nesting of
``jastrow``, ``correlation``, and ``coefficient`` elements.


Input Specification
^^^^^^^^^^^^^^^^^^^
Jastrow element:

    +----------+--------------+------------+--------------+----------------+
    | **name** | **datatype** | **values** | **defaults** | **description**|
    |          |              |            |              |                |
    +----------+--------------+------------+--------------+----------------+
    | name     | text         |            | (required)   | Unique name    |
    |          |              |            |              | for this       |
    |          |              |            |              | Jastrow        |
    |          |              |            |              | function       |
    +----------+--------------+------------+--------------+----------------+
    | type     | text         | One-body   | (required)   | Define a       |
    |          |              |            |              | one-body       |
    |          |              |            |              | function       |
    +----------+--------------+------------+--------------+----------------+
    | function | text         | Bspline    | (required)   | BSpline        |
    |          |              |            |              | Jastrow        |
    +----------+--------------+------------+--------------+----------------+
    |          | text         | pade2      |              | Pade form      |
    +----------+--------------+------------+--------------+----------------+
    |          | text         | …          |              | …              |
    +----------+--------------+------------+--------------+----------------+
    | source   | text         | name       | (required)   | Name of        |
    |          |              |            |              | attribute of   |
    |          |              |            |              | classical      |
    |          |              |            |              | particle set   |
    +----------+--------------+------------+--------------+----------------+
    | print    | text         | yes / no   | yes          | Jastrow        |
    |          |              |            |              | factor         |
    |          |              |            |              | printed in     |
    |          |              |            |              | external       |
    |          |              |            |              | file?          |
    +----------+--------------+------------+--------------+----------------+

    +----------+--------------+------------+--------------+--------------+
    | elements |              |            |              |              |
    +----------+--------------+------------+--------------+--------------+
    |          | Correlation  |            |              |              |
    +----------+--------------+------------+--------------+--------------+
    | Contents |              |            |              |              |
    +----------+--------------+------------+--------------+--------------+
    |          | (None)       |            |              |              |
    +----------+--------------+------------+--------------+--------------+

To be more concrete, the one-body Jastrow factors used to describe correlations
between electrons and ions take the form below:

.. math::
  :label: eq9

  J1=\sum_I^{ion0}\sum_i^e u_{ab}(|r_i-R_I|)

where I runs over all of the ions in the calculation, i runs over the
electrons and :math:`u_{ab}` describes the functional form of the
correlation between them. Many different forms of :math:`u_{ab}` are
implemented in QMCPACK. We will detail two of the most common ones
below.

.. _onebodyjastrowspline:

Spline form
...........

The one-body spline Jastrow function is the most commonly used one-body
Jastrow for solids. This form was first described and used in
:cite:`EslerKimCeperleyShulenburger2012`. Here
:math:`u_{ab}` is an interpolating 1D B-spline (tricublc spline on a
linear grid) between zero distance and :math:`r_{cut}`. In 3D periodic
systems the default cutoff distance is the Wigner Seitz cell radius. For
other periodicities, including isolated molecules, the :math:`r_{cut}`
must be specified. The cusp can be set. :math:`r_i` and :math:`R_I` are
most commonly the electron and ion positions, but any particlesets that
can provide the needed centers can be used.

Correlation element:

    +-------------+-------------+-------------+-------------+----------------+
    | **Name**    | **Datatype**| **Values**  | **Defaults**| **Description**|
    |             |             |             |             |                |
    +-------------+-------------+-------------+-------------+----------------+
    | ElementType | Text        | Name        | See below   | Classical      |
    |             |             |             |             | particle       |
    |             |             |             |             | target         |
    +-------------+-------------+-------------+-------------+----------------+
    | SpeciesA    | Text        | Name        | See below   | Classical      |
    |             |             |             |             | particle       |
    |             |             |             |             | target         |
    +-------------+-------------+-------------+-------------+----------------+
    | SpeciesB    | Text        | Name        | See below   | Quantum        |
    |             |             |             |             | species        |
    |             |             |             |             | target         |
    +-------------+-------------+-------------+-------------+----------------+
    | Size        | Integer     | :math:`> 0` | (Required)  | Number of      |
    |             |             |             |             | coefficients   |
    |             |             |             |             |                |
    +-------------+-------------+-------------+-------------+----------------+
    | Rcut        | Real        | :math:`> 0` | See below   | Distance at    |
    |             |             |             |             | which the      |
    |             |             |             |             | correlation    |
    |             |             |             |             | goes to 0      |
    +-------------+-------------+-------------+-------------+----------------+
    | Cusp        | Real        |:math:`\ge 0`| 0           | Value for      |
    |             |             |             |             | use in Kato    |
    |             |             |             |             | cusp           |
    |             |             |             |             | condition      |
    +-------------+-------------+-------------+-------------+----------------+
    | Spin        | Text        | Yes or no   | No          | Spin           |
    |             |             |             |             | dependent      |
    |             |             |             |             | Jastrow        |
    |             |             |             |             | factor         |
    +-------------+-------------+-------------+-------------+----------------+

    +----------+--------------+------------+--------------+--------------+
    | Elements |              |            |              |              |
    +----------+--------------+------------+--------------+--------------+
    |          | Coefficients |            |              |              |
    +----------+--------------+------------+--------------+--------------+
    | Contents |              |            |              |              |
    +----------+--------------+------------+--------------+--------------+
    |          | (None)       |            |              |              |
    +----------+--------------+------------+--------------+--------------+

Additional information:

- ``elementType, speciesA, speciesB, spin``
    For a spin-independent Jastrow factor (spin = “no”), elementType
    should be the name of the group of ions in the classical particleset to
    which the quantum particles should be correlated. For a spin-dependent
    Jastrow factor (spin = “yes”), set speciesA to the group name in the
    classical particleset and speciesB to the group name in the quantum
    particleset.

- ``rcut``
    The cutoff distance for the function in atomic units (bohr). For 3D
    fully periodic systems, this parameter is optional, and a default of the
    Wigner Seitz cell radius is used. Otherwise this parameter is required.

- ``cusp``
    The one-body Jastrow factor can be used to make the wavefunction
    satisfy the electron-ion cusp condition :cite:``kato``. In this
    case, the derivative of the Jastrow factor as the electron approaches
    the nucleus will be given by

.. math::
  :label: eq10

  \left(\frac{\partial J}{\partial r_{iI}}\right)_{r_{iI} = 0} = -Z .

Note that if the antisymmetric part of the wavefunction satisfies the electron-ion cusp
condition (for instance by using single-particle orbitals that respect the cusp condition)
or if a nondivergent pseudopotential is used, the Jastrow should be cuspless at the
nucleus and this value should be kept at its default of 0.

Coefficients element:

    +-----------+--------------+------------+--------------+----------------+
    | **Name**  | **Datatype** | **Values** | **Defaults** | **Description**|
    |           |              |            |              |                |
    +-----------+--------------+------------+--------------+----------------+
    | Id        | Text         |            | (Required)   | Unique         |
    |           |              |            |              | identifier     |
    +-----------+--------------+------------+--------------+----------------+
    | Type      | Text         | Array      | (Required)   |                |
    +-----------+--------------+------------+--------------+----------------+
    | Optimize  | Text         | Yes or no  | Yes          | if no,         |
    |           |              |            |              | values are     |
    |           |              |            |              | fixed in       |
    |           |              |            |              | optimizations  |
    |           |              |            |              |                |
    +-----------+--------------+------------+--------------+----------------+
    +-----------+--------------+------------+--------------+----------------+
    | Elements  |              |            |              |                |
    +-----------+--------------+------------+--------------+----------------+
    | (None)    |              |            |              |                |
    +-----------+--------------+------------+--------------+----------------+
    | Contents  |              |            |              |                |
    +-----------+--------------+------------+--------------+----------------+
    | (No name) | Real array   |            | Zeros        | Jastrow        |
    |           |              |            |              | coefficients   |
    +-----------+--------------+------------+--------------+----------------+

.. _1bjsplineexamples:

Example use cases
.................

Specify a spin-independent function with four parameters. Because rcut  is not
specified, the default cutoff of the Wigner Seitz cell radius is used; this
Jastrow must be used with a 3D periodic system such as a bulk solid. The name of
the particleset holding the ionic positions is "i."

::

  <jastrow name="J1" type="One-Body" function="Bspline" print="yes" source="i">
   <correlation elementType="C" cusp="0.0" size="4">
     <coefficients id="C" type="Array"> 0  0  0  0  </coefficients>
   </correlation>
  </jastrow>

Specify a spin-dependent function with seven up-spin and seven down-spin parameters.
The cutoff distance is set to 6 atomic units.  Note here that the particleset holding
the ions is labeled as ion0 rather than "i," as in the other example.  Also in this case,
the ion is lithium with a coulomb potential, so the cusp condition is satisfied by
setting cusp="d."

::

  <jastrow name="J1" type="One-Body" function="Bspline" source="ion0" spin="yes">
    <correlation speciesA="Li" speciesB="u" size="7" rcut="6">
      <coefficients id="eLiu" cusp="3.0" type="Array">
      0.0 0.0 0.0 0.0 0.0 0.0 0.0
      </coefficients>
    </correlation>
    <correlation speciesA="C" speciesB="d" size="7" rcut="6">
      <coefficients id="eLid" cusp="3.0" type="Array">
      0.0 0.0 0.0 0.0 0.0 0.0 0.0
      </coefficients>
    </correlation>
  </jastrow>

.. _onebodyjastrowpade:

Pade form
.........

Although the spline Jastrow factor is the most flexible and most commonly used form implemented in QMCPACK,
there are times where its flexibility can make it difficult to optimize.  As an example, a spline Jastrow
with a very large cutoff can be difficult to optimize for isolated systems such as molecules because of the small
number of samples present in the tail of the function.  In such cases, a simpler functional
form might be advantageous.  The second-order Pade Jastrow factor, given in :eq:`eq11`, is a good choice
in such cases.

.. math::
  :label: eq11

  u_{ab}(r) = \frac{a*r+c*r^2}{1+b*r}

Unlike the spline Jastrow factor, which includes a cutoff, this form has an infinite range and will be applied to every particle
pair (subject to the minimum image convention).  It also is a cuspless Jastrow factor,
so it should be used either in combination with a single particle basis set that contains the proper cusp or
with a smooth pseudopotential.

Correlation element:

      +-------------+--------------+------------+--------------+---------------------------+
      | **Name**    | **Datatype** | **Values** | **Defaults** | **Description**           |
      +-------------+--------------+------------+--------------+---------------------------+
      | ElementType | Text         | Name       | See below    | Classical particle target |
      +-------------+--------------+------------+--------------+---------------------------+
      | Elements    |              |            |              |                           |
      +-------------+--------------+------------+--------------+---------------------------+
      |             | Coefficients |            |              |                           |
      +-------------+--------------+------------+--------------+---------------------------+
      | Contents    |              |            |              |                           |
      +-------------+--------------+------------+--------------+---------------------------+
      |             | (None)       |            |              |                           |
      +-------------+--------------+------------+--------------+---------------------------+

Parameter element:

      +-----------+-------------+-------------+-------------+-----------------+
      | **Name**  |**Datatype** | **Values**  | **Defaults**| **Description** |
      |           |             |             |             |                 |
      +-----------+-------------+-------------+-------------+-----------------+
      | Id        | String      | Name        | (Required)  | Name for        |
      |           |             |             |             | variable        |
      +-----------+-------------+-------------+-------------+-----------------+
      | Name      | String      | A or B or C | (Required)  | See             |
      |           |             |             |             | :eq:`eq11`      |
      |           |             |             |             |                 |
      +-----------+-------------+-------------+-------------+-----------------+
      | Optimize  | Text        | Yes or no   | Yes         | If no,          |
      |           |             |             |             | values are      |
      |           |             |             |             | fixed in        |
      |           |             |             |             | optimizations   |
      |           |             |             |             |                 |
      +-----------+-------------+-------------+-------------+-----------------+

      +-----------+-------------+-------------+-------------+-------------+
      | Elements  |             |             |             |             |
      +-----------+-------------+-------------+-------------+-------------+
      | (None)    |             |             |             |             |
      +-----------+-------------+-------------+-------------+-------------+
      | Contents  |             |             |             |             |
      +-----------+-------------+-------------+-------------+-------------+
      | (No name) | Real        | Parameter   | (Required)  | Jastrow     |
      |           |             | value       |             | coefficients|
      |           |             |             |             |             |
      +-----------+-------------+-------------+-------------+-------------+

.. _1bjpadeexamples:

Example use case
................

Specify a spin-independent function with independent Jastrow factors for two different species (Li and H).
The name of the particleset holding the ionic positions is "i."

::

  <jastrow name="J1" function="pade2" type="One-Body" print="yes" source="i">
    <correlation elementType="Li">
      <var id="LiA" name="A">  0.34 </var>
      <var id="LiB" name="B"> 12.78 </var>
      <var id="LiC" name="C">  1.62 </var>
    </correlation>
    <correlation elementType="H"">
      <var id="HA" name="A">  0.14 </var>
      <var id="HB" name="B"> 6.88 </var>
      <var id="HC" name="C"> 0.237 </var>
    </correlation>
  </jastrow>

.. onebodyjastrowsrcusp:

Short Range Cusp Form
.....................

The idea behind this functor is to encode nuclear cusps and other details at very
short range around a nucleus in the region that the Gaussian orbitals of quantum
chemistry are not capable of describing correctly.
The functor is kept short ranged, because outside this small region, quantum chemistry
orbital expansions are already capable of taking on the correct shapes.
Unlike a pre-computed cusp correction, this optimizable functor can respond to
changes in the wave function during VMC optimization.
The functor's form is

.. math::
  :label: eq12

  u(r) = -\exp{\left(-r/R_0\right)} \left( A R_0 + \sum_{k=0}^{N-1} B_k \frac{ (r/R_0)^{k+2} }{ 1 + (r/R_0)^{k+2} } \right)

in which :math:`R_0` acts as a soft cutoff radius (:math:`u(r)` decays to zero quickly beyond roughly this distance)
and :math:`A` determines the cusp condition.

.. math::
  :label: eq13

  \lim_{r \to 0} \frac{\partial u}{\partial r} = A

The simple exponential decay is modified by the :math:`N` coefficients
:math:`B_k` that define an expansion in sigmoidal functions, thus adding
detailed structure in a short-ranged region around a nucleus while
maintaining the correct cusp condition at the nucleus. Note that
sigmoidal functions are used instead of, say, a bare polynomial
expansion, as they trend to unity past the soft cutoff radius and so
interfere less with the exponential decay that keeps the functor short
ranged. Although :math:`A`, :math:`R_0`, and the :math:`B_k`
coefficients can all be optimized as variational parameters, :math:`A`
will typically be fixed as the desired cusp condition is known.

To specify this one-body Jastrow factor, use an input section like the following.

::

  <jastrow name="J1Cusps" type="One-Body" function="shortrangecusp" source="ion0" print="yes">
    <correlation rcut="6" cusp="3" elementType="Li">
      <var id="LiCuspR0" name="R0" optimize="yes"> 0.06 </var>
      <coefficients id="LiCuspB" type="Array" optimize="yes">
        0 0 0 0 0 0 0 0 0 0
      </coefficients>
    </correlation>
    <correlation rcut="6" cusp="1" elementType="H">
      <var id="HCuspR0" name="R0" optimize="yes"> 0.2 </var>
      <coefficients id="HCuspB" type="Array" optimize="yes">
        0 0 0 0 0 0 0 0 0 0
      </coefficients>
    </correlation>
  </jastrow>

Here “rcut” is specified as the range beyond which the functor is
assumed to be zero. The value of :math:`A` can either be specified via
the “cusp” option as shown above, in which case its optimization is
disabled, or through its own “var” line as for :math:`R_0`, in which
case it can be specified as either optimizable (“yes”) or not (“no”).
The coefficients :math:`B_k` are specified via the “coefficients”
section, with the length :math:`N` of the expansion determined
automatically based on the length of the array.


Note that this one-body Jastrow form can (and probably should) be used in conjunction
with a longer ranged one-body Jastrow, such as a spline form.
Be sure to set the longer-ranged Jastrow to be cusp-free!

Two-body Jastrow functions
~~~~~~~~~~~~~~~~~~~~~~~~~~

The two-body Jastrow factor is a form that allows for the explicit inclusion
of dynamic correlation between two particles included in the wavefunction.  It
is almost always given in a spin dependent form so as to satisfy the Kato cusp
condition between electrons of different spins :cite:`kato`.

The two body Jastrow function is specified within a ``wavefunction`` element
and must contain one or more correlation elements specifying additional parameters
as well as the actual coefficients.  :ref:`2bjsplineexamples` gives
examples of the typical nesting of ``jastrow``, ``correlation`` and
``coefficient`` elements.

Input Specification
^^^^^^^^^^^^^^^^^^^

Jastrow element:

    +----------+--------------+------------+--------------+-----------------+
    | **name** | **datatype** | **values** | **defaults** | **description** |
    |          |              |            |              |                 |
    +----------+--------------+------------+--------------+-----------------+
    | name     | text         |            | (required)   | Unique name     |
    |          |              |            |              | for this        |
    |          |              |            |              | Jastrow         |
    |          |              |            |              | function        |
    +----------+--------------+------------+--------------+-----------------+
    | type     | text         | Two-body   | (required)   | Define a        |
    |          |              |            |              | one-body        |
    |          |              |            |              | function        |
    +----------+--------------+------------+--------------+-----------------+
    | function | text         | Bspline    | (required)   | BSpline         |
    |          |              |            |              | Jastrow         |
    +----------+--------------+------------+--------------+-----------------+
    | print    | text         | yes / no   | yes          | Jastrow         |
    |          |              |            |              | factor          |
    |          |              |            |              | printed in      |
    |          |              |            |              | external        |
    |          |              |            |              | file?           |
    +----------+--------------+------------+--------------+-----------------+
    +----------+--------------+------------+--------------+-----------------+
    | elements |              |            |              |                 |
    +----------+--------------+------------+--------------+-----------------+
    |          | Correlation  |            |              |                 |
    +----------+--------------+------------+--------------+-----------------+
    | Contents |              |            |              |                 |
    +----------+--------------+------------+--------------+-----------------+
    |          | (None)       |            |              |                 |
    +----------+--------------+------------+--------------+-----------------+

The two-body Jastrow factors used to describe correlations between electrons take the form

.. math::
  :label: eq14

  J2=\sum_i^{e}\sum_{j>i}^{e} u_{ab}(|r_i-r_j|)

The most commonly used form of two body Jastrow factor supported by the code is a splined
Jastrow factor, with many similarities to the one body spline Jastrow.

.. _twobodyjastrowspline:

Spline form
...........

The two-body spline Jastrow function is the most commonly used two-body
Jastrow for solids. This form was first described and used in
:cite:`EslerKimCeperleyShulenburger2012`. Here
:math:`u_{ab}` is an interpolating 1D B-spline (tricublc spline on a
linear grid) between zero distance and :math:`r_{cut}`. In 3D periodic
systems, the default cutoff distance is the Wigner Seitz cell radius.
For other periodicities, including isolated molecules, the
:math:`r_{cut}` must be specified. :math:`r_i` and :math:`r_j` are
typically electron positions. The cusp condition as :math:`r_i`
approaches :math:`r_j` is set by the relative spin of the electrons.

Correlation element:

    +----------+-------------+-------------+-------------+-----------------+
    | **Name** | **Datatype**| **Values**  | **Defaults**| **Description** |
    |          |             |             |             |                 |
    +----------+-------------+-------------+-------------+-----------------+
    | SpeciesA | Text        | U or d      | (Required)  | Quantum         |
    |          |             |             |             | species         |
    |          |             |             |             | target          |
    +----------+-------------+-------------+-------------+-----------------+
    | SpeciesB | Text        | U or d      | (Required)  | Quantum         |
    |          |             |             |             | species         |
    |          |             |             |             | target          |
    +----------+-------------+-------------+-------------+-----------------+
    | Size     | Integer     | :math:`> 0` | (Required)  | Number of       |
    |          |             |             |             | coefficients    |
    |          |             |             |             |                 |
    +----------+-------------+-------------+-------------+-----------------+
    | Rcut     | Real        | :math:`> 0` | See below   | Distance at     |
    |          |             |             |             | which the       |
    |          |             |             |             | correlation     |
    |          |             |             |             | goes to 0       |
    +----------+-------------+-------------+-------------+-----------------+
    | Spin     | Text        | Yes or no   | No          | Spin-dependent  |
    |          |             |             |             | Jastrow factor  |
    |          |             |             |             |                 |
    |          |             |             |             |                 |
    +----------+-------------+-------------+-------------+-----------------+
    +----------+-------------+-------------+-------------+-----------------+
    |Elements  |             |             |             |                 |
    +----------+-------------+-------------+-------------+-----------------+
    |          | Coefficients|             |             |                 |
    |          |             |             |             |                 |
    +----------+-------------+-------------+-------------+-----------------+
    |Contents  |             |             |             |                 |
    +----------+-------------+-------------+-------------+-----------------+
    |          | (None)      |             |             |                 |
    +----------+-------------+-------------+-------------+-----------------+

Additional information:

- ``speciesA, speciesB`` The scale function u(r) is defined for species pairs uu and ud.
  There is no need to define ud or dd since uu=dd and ud=du.  The cusp condition is computed internally
  based on the charge of the quantum particles.

Coefficients element:

    +-----------+--------------+------------+--------------+-----------------+
    | **Name**  | **Datatype** | **Values** | **Defaults** | **Description** |
    |           |              |            |              |                 |
    +-----------+--------------+------------+--------------+-----------------+
    | Id        | Text         |            | (Required)   | Unique          |
    |           |              |            |              | identifier      |
    +-----------+--------------+------------+--------------+-----------------+
    | Type      | Text         | Array      | (Required)   |                 |
    +-----------+--------------+------------+--------------+-----------------+
    | Optimize  | Text         | Yes or no  | Yes          | If no,          |
    |           |              |            |              | values are      |
    |           |              |            |              | fixed in        |
    |           |              |            |              | optimizations   |
    |           |              |            |              |                 |
    +-----------+--------------+------------+--------------+-----------------+
    +-----------+--------------+------------+--------------+-----------------+
    |Elements   |              |            |              |                 |
    +-----------+--------------+------------+--------------+-----------------+
    | (None)    |              |            |              |                 |
    +-----------+--------------+------------+--------------+-----------------+
    | Contents  |              |            |              |                 |
    +-----------+--------------+------------+--------------+-----------------+
    | (No name) | Real array   |            | Zeros        | Jastrow         |
    |           |              |            |              | coefficients    |
    +-----------+--------------+------------+--------------+-----------------+

.. _2bjsplineexamples:

Example use cases
.................

Specify a spin-dependent function with four parameters for each channel.  In this case, the cusp is set at
a radius of 4.0 bohr (rather than to the default of the Wigner Seitz cell radius).  Also, in this example,
the coefficients are set to not be optimized during an optimization step.

::

  <jastrow name="J2" type="Two-Body" function="Bspline" print="yes">
    <correlation speciesA="u" speciesB="u" size="8" rcut="4.0">
      <coefficients id="uu" type="Array" optimize="no"> 0.2309049836 0.1312646071 0.05464141356 0.01306231516</coefficients>
    </correlation>
    <correlation speciesA="u" speciesB="d" size="8" rcut="4.0">
      <coefficients id="ud" type="Array" optimize="no"> 0.4351561096 0.2377951747 0.1129144262 0.0356789236</coefficients>
    </correlation>
  </jastrow>

.. _jastrowuserform:

User defined functional form
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To aid in implementing different forms for :math:`u_{ab}(r)`, there is a
script that uses a symbolic expression to generate the appropriate code
(with spatial and parameter derivatives). The script is located in
``src/QMCWaveFunctions/Jastrow/codegen/user_jastrow.py``. The script
requires Sympy (www.sympy.org) for symbolic mathematics and code
generation.

To use the script, modify it to specify the functional form and a list
of variational parameters. Optionally, there may be fixed parameters -
ones that are specified in the input file, but are not part of the
variational optimization. Also one symbol may be specified that accepts
a cusp value in order to satisfy the cusp condition. There are several
example forms in the script. The default form is the simple Padé.

Once the functional form and parameters are specified in the script, run
the script from the ``codegen`` directory and recompile QMCPACK. The
main output of the script is the file
``src/QMCWaveFunctions/Jastrow/UserFunctor.h``. The script also prints
information to the screen, and one section is a sample XML input block
containing all the parameters.

There is a unit test in
``src/QMCWaveFunctions/test/test_user_jastrow.cpp`` to perform some
minimal testing of the Jastrow factor. The unit test will need updating
to properly test new functional forms. Most of the changes relate to the
number and name of variational parameters.

Jastrow element:

    +----------+--------------+------------+--------------+----------------+
    | **name** | **datatype** | **values** | **defaults** | **description**|
    |          |              |            |              |                |
    +----------+--------------+------------+--------------+----------------+
    | name     | text         |            | (required)   | Unique name    |
    |          |              |            |              | for this       |
    |          |              |            |              | Jastrow        |
    |          |              |            |              | function       |
    +----------+--------------+------------+--------------+----------------+
    | type     | text         | One-body   | (required)   | Define a       |
    |          |              |            |              | one-body       |
    |          |              |            |              | function       |
    +----------+--------------+------------+--------------+----------------+
    |          |              | Two-body   | (required)   | Define a       |
    |          |              |            |              | two-body       |
    |          |              |            |              | function       |
    +----------+--------------+------------+--------------+----------------+
    | function | text         | user       | (required)   | User-defined   |
    |          |              |            |              | functor        |
    +----------+--------------+------------+--------------+----------------+

    See other parameters as appropriate for one or two-body functions

    +----------+--------------+------------+--------------+--------------+
    | elements |              |            |              |              |
    +----------+--------------+------------+--------------+--------------+
    |          | Correlation  |            |              |              |
    +----------+--------------+------------+--------------+--------------+
    | Contents |              |            |              |              |
    +----------+--------------+------------+--------------+--------------+
    |          | (None)       |            |              |              |
    +----------+--------------+------------+--------------+--------------+

Long-ranged Jastrow factors
~~~~~~~~~~~~~~~~~~~~~~~~~~~

While short-ranged Jastrow factors capture the majority of the benefit
for minimizing the total energy and the energy variance, long-ranged
Jastrow factors are important to accurately reproduce the short-ranged
(long wavelength) behavior of quantities such as the static structure
factor, and are therefore essential for modern accurate finite size
corrections in periodic systems.

Below two types of long-ranged Jastrow factors are described. The first
(the k-space Jastrow) is simply an expansion of the one and/or two body
correlation functions in plane waves, with the coefficients comprising
the optimizable parameters. The second type have few variational
parameters and use the optimized breakup method of Natoli and
Ceperley :cite:`Natoli1995` (the Yukawa and Gaskell RPA
Jastrows).

Long-ranged Jastrow: k-space Jastrow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The k-space Jastrow introduces explicit long-ranged dependence commensurate with the periodic supercell.  This Jastrow is to be used in periodic boundary conditions only.

The input for the k-space Jastrow fuses both one and two-body forms into a single element and so they are discussed together here.  The one- and two-body terms in the k-Space Jastrow have the form:

.. math::
  :label: eq15

  J_1 = \sum_{G\ne 0}b_G\rho_G^I\rho_{-G}

.. math::
  :label: eq16

  J_2 = \sum_{G\ne 0}a_G\rho_G\rho_{-G}

Here :math:`\rho_G` is the Fourier transform of the instantaneous electron density:

.. math::
  :label: eq17

  \rho_G=\sum_{n\in electrons}e^{iG\cdot r_n}

and :math:`\rho_G^I` has the same form, but for the fixed ions. In both cases the coefficients are restricted to be real, though in general the coefficients for the one-body term need not be.  See :ref:`feature-kspace-jastrow` for more detail.

Input for the k-space Jastrow follows the familar nesting of ``jastrow-correlation-coefficients`` elements, with attributes unique to the k-space Jastrow at the ``correlation`` input level.

``jastrow type=kSpace`` element:

    +------------------+------------------+
    | parent elements: | ``wavefunction`` |
    +------------------+------------------+
    | child elements:  | ``correlation``  |
    +------------------+------------------+

attributes:

    +-----------------------+--------------+----------------------+-------------+--------------------------+
    | **Name**              | **Datatype** | **Values**           | **Default** | **Description**          |
    +=======================+==============+======================+=============+==========================+
    | ``type``:math:`^r`    | text         | **kSpace**           |             | must be kSpace           |
    +-----------------------+--------------+----------------------+-------------+--------------------------+
    | ``name``:math:`^r`    | text         | *anything*           | 0           | Unique name for Jastrow  |
    +-----------------------+--------------+----------------------+-------------+--------------------------+
    | ``source``:math:`^r`  | text         | ``particleset.name`` |             | Ion particleset name     |
    +-----------------------+--------------+----------------------+-------------+--------------------------+

``correlation`` element:

    +------------------+-------------------------+
    | parent elements: | ``jastrow type=kSpace`` |
    +------------------+-------------------------+
    | child elements:  | ``coefficients``        |
    +------------------+-------------------------+


attributes:


    +------------------------------+--------------+--------------------------+-------------+---------------------------+
    | **Name**                     | **Datatype** | **Values**               | **Default** | **Description**           |
    +==============================+==============+==========================+=============+===========================+
    | ``type``:math:`^r`           | text         | **One-body, Two-Body**   |             | Must be One-body/Two-body |
    +------------------------------+--------------+--------------------------+-------------+---------------------------+
    | ``kc``:math:`^r`             | real         | kc :math:`\ge` 0         | 0.0         | k-space cutoff in a.u.    |
    +------------------------------+--------------+--------------------------+-------------+---------------------------+
    | ``symmetry``:math:`^o`       | text         | crystal,isotropic,none   | crystal     | symmetry of coefficients  |
    +------------------------------+--------------+--------------------------+-------------+---------------------------+
    | ``spinDependent``:math:`^o`  | boolean      | yes,no                   | no          | *No current function*     |
    +------------------------------+--------------+--------------------------+-------------+---------------------------+

``coefficients`` element:

    +------------------+-----------------+
    | parent elements: | ``correlation`` |
    +------------------+-----------------+
    | child elements:  | *None*          |
    +------------------+-----------------+

attributes:

    +---------------------+--------------+------------+-------------+------------------+
    | **Name**            | **Datatype** | **Values** | **Default** | **Description**  |
    +=====================+==============+============+=============+==================+
    | ``id``:math:`^r`    | text         | *anything* | cG1/cG2     | Label for coeffs |
    +---------------------+--------------+------------+-------------+------------------+
    | ``type``:math:`^r`  | text         | ``Array``  | 0           | Must be Array    |
    +---------------------+--------------+------------+-------------+------------------+

    body text: The body text is a list of real values for the parameters.

Additional information:

-  It is normal to provide no coefficients as an initial guess. The
   number of coefficients will be automatically calculated according to
   the k-space cutoff + symmetry and set to zero.

-  Providing an incorrect number of parameters also results in all
   parameters being set to zero.

-  There is currently no way to turn optimization on/off for the k-space
   Jastrow. The coefficients are always optimized.

-  Spin dependence is currently not implemented for this Jastrow.

-  ``kc``: Parameters with G vectors magnitudes less than ``kc`` are
   included in the Jastrow. If ``kc`` is zero, it is the same as
   excluding the k-space term.

-  ``symmetry=crystal``: Impose crystal symmetry on coefficients
   according to the structure factor.

-  ``symmetry=isotropic``: Impose spherical symmetry on coefficients
   according to G-vector magnitude.

-  ``symmetry=none``: Impose no symmetry on the coefficients.

.. code-block::
  :caption: k-space Jastrow with one- and two-body terms.
  :name: Listing 10

  <jastrow type="kSpace" name="Jk" source="ion0">
    <correlation kc="4.0" type="One-Body" symmetry="cystal">
      <coefficients id="cG1" type="Array">
      </coefficients>
    </correlation>
    <correlation kc="4.0" type="Two-Body" symmetry="crystal">
      <coefficients id="cG2" type="Array">
      </coefficients>
   </correlation>
  </jastrow>

.. code-block::
  :caption: k-space Jastrow with one-body term only.
  :name: Listing 11

  <jastrow type="kSpace" name="Jk" source="ion0">
     <correlation kc="4.0" type="One-Body" symmetry="crystal">
        <coefficients id="cG1" type="Array">
        </coefficients>
     </correlation>
  </jastrow>

.. code-block::
  :caption: k-space Jastrow with two-body term only.
  :name: Listing 12

  <jastrow type="kSpace" name="Jk" source="ion0">
     <correlation kc="4.0" type="Two-Body" symmetry="crystal">
        <coefficients id="cG2" type="Array">
        </coefficients>
     </correlation>
  </jastrow>

.. _twobodyjastrowlr:

Long-ranged Jastrows: Gaskell RPA and Yukawa forms
..................................................

**NOTE: The Yukawa and RPA Jastrows do not work at present
and are currently being revived.  Please contact the developers if
you are interested in using them.**

The exact Jastrow correlation functions contain terms which have a
form similar to the Coulomb pair potential.  In periodic systems
the Coulomb potential is replaced by an Ewald summation of the
bare potential over all periodic image cells.  This sum is often
handled by the optimized breakup method :cite:`Natoli1995` and this
same approach is applied to the long-ranged Jastrow factors in QMCPACK.

There are two main long-ranged Jastrow factors of this type
implemented in QMCPACK: the Gaskell RPA :cite:`Gaskell1961,Gaskell1962`
form and the :cite:`Ceperley1978` form.  Both of these forms
were used by Ceperley in early studies of the electron gas :cite:`Ceperley1978`,
but they are also appropriate starting points for general solids.

The Yukawa form is defined in real space.  It's long-range form is
formally defined as

.. math::
  :label: eq18

  u_Y^{PBC}(r) = \sum_{L\ne 0}\sum_{i<j}u_Y(\left|{r_i-r_j+L}\right|)

with :math:`u_Y(r)` given by

.. math::
  :label: eq19

  u_Y(r) = \frac{a}{r}\left(1-e^{-r/b}\right)

In QMCPACK a slightly more restricted form is used:

.. math::
  :label: eq20

  u_Y(r) = \frac{r_s}{r}\left(1-e^{-r/\sqrt{r_s}}\right)

here ":math:`r_s`" is understood to be a variational parameter.

The Gaskell RPA form---which contains correct short/long range limits
and minimizes the total energy of the electron gas within the RPA---is
defined directly in k-space:

.. math::
  :label: eq21

  u_{RPA}(k) = -\frac{1}{2S_0(k)}+\frac{1}{2}\left(\frac{1}{S_0(k)^2}+\frac{4m_ev_k}{\hbar^2k^2}\right)^{1/2}

where $v_k$ is the Fourier transform of the Coulomb potential and
:math:`S_0(k)` is the static structure factor of the non-interacting
electron gas:

.. math::

  S_0(k) = \left.
      \begin{cases}
        1 &  k>2k_F \\
        \frac{3k}{4k_F}-\frac{1}{2}\left(\frac{k}{2k_F}\right)^3 & k<2k_F
      \end{cases}
      \right.

When written in atomic units, RPA Jastrow implemented in QMCPACK has the
form

.. math::
  :label: eq22

  u_{RPA}(k) = \frac{1}{2N_e}\left(-\frac{1}{S_0(k)}+\left(\frac{1}{S_0(k)^2}+\frac{12}{r_s^3k^4}\right)^{1/2}\right)

Here ":math:`r_s`" is again a variational parameter and :math:`k_F\equiv(\tfrac{9\pi}{4r_s^3})^{1/3}`.

For both the Yukawa and Gaskell RPA Jastrows, the default value for :math:`r_s` is :math:`r_s=(\tfrac{3\Omega}{4\pi N_e})^{1/3}`.

``jastrow type=Two-Body function=rpa/yukawa`` element:

    +------------------+-----------------+
    | parent elements: | ``wavefunction``|
    +------------------+-----------------+
    | child elements:  | ``correlation`` |
    +------------------+-----------------+

  attributes:

    +--------------------------+--------------+----------------+-------------+--------------------------+
    | **Name**                 | **Datatype** | **Values**     | **Default** | **Description**          |
    +==========================+==============+================+=============+==========================+
    | ``type``:math:`^r`       | text         | **Two-body**   |             | Must be two-body         |
    +--------------------------+--------------+----------------+-------------+--------------------------+
    | ``function``:math:`^r`   | text         | **rpa/yukawa** |             | Must be rpa or yukawa    |
    +--------------------------+--------------+----------------+-------------+--------------------------+
    | ``name``:math:`^r`       | text         | *anything*     | RPA_Jee     | Unique name for Jastrow  |
    +--------------------------+--------------+----------------+-------------+--------------------------+
    | ``longrange``:math:`^o`  | boolean      | yes/no         | yes         | Use long-range part      |
    +--------------------------+--------------+----------------+-------------+--------------------------+
    | ``shortrange``:math:`^o` | boolean      | yes/no         | yes         | Use short-range part     |
    +--------------------------+--------------+----------------+-------------+--------------------------+

  parameters:

    +-------------------+--------------+----------------+----------------------------------------------------------------------+-------------------------+
    | **Name**          | **Datatype** | **Values**     | **Default**                                                          | **Description**         |
    +===================+==============+================+======================================================================+=========================+
    | ``rs``:math:`^o`  | rs           | :math:`r_s>0`  | :math:`\tfrac{3\Omega}{4\pi N_e}`                                    | Avg. elec-elec distance |
    +-------------------+--------------+----------------+----------------------------------------------------------------------+-------------------------+
    | ``kc``:math:`^o`  | kc           | :math:`k_c>0`  | :math:`2\left(\tfrac{9\pi}{4}\right)^{1/3}\tfrac{4\pi N_e}{3\Omega}` | k-space cutoff          |
    +-------------------+--------------+----------------+----------------------------------------------------------------------+-------------------------+

.. code-block::
  :caption: Two body RPA Jastrow with long- and short-ranged parts.
  :name: Listing 13

  <jastrow name=''Jee'' type=''Two-Body'' function=''rpa''>
  </jastrow>

Three-body Jastrow functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Explicit three-body correlations can be included in the wavefunction via
the three-body Jastrow factor. The three-body electron-electron-ion
correlation function (:math:`u_{\sigma\sigma'I}`) currently used in is
identical to the one proposed in :cite:`Drummond2004`:

.. math::
  :label: eq23

   \begin{aligned}
   u_{\sigma\sigma'I}(r_{\sigma I},r_{\sigma'I},r_{\sigma\sigma'}) &= \sum_{\ell=0}^{M_{eI}}\sum_{m=0}^{M_{eI}}\sum_{n=0}^{M_{ee}}\gamma_{\ell mn} r_{\sigma I}^\ell r_{\sigma'I}^m r_{\sigma\sigma'}^n \\
      &\times \left(r_{\sigma I}-\frac{r_c}{2}\right)^3 \Theta\left(r_{\sigma I}-\frac{r_c}{2}\right) \nonumber \\
      &\times \left(r_{\sigma' I}-\frac{r_c}{2}\right)^3 \Theta\left(r_{\sigma' I}-\frac{r_c}{2}\right) \nonumber\end{aligned}

Here :math:`M_{eI}` and :math:`M_{ee}` are the maximum polynomial orders
of the electron-ion and electron-electron distances, respectively,
:math:`\{\gamma_{\ell mn}\}` are the optimizable parameters (modulo
constraints), :math:`r_c` is a cutoff radius, and :math:`r_{ab}` are the
distances between electrons or ions :math:`a` and :math:`b`. i.e. The
correlation function is only a function of the interparticle distances
and not a more complex function of the particle positions,
:math:`\mathbf{r}`. As indicated by the :math:`\Theta` functions,
correlations are set to zero beyond a distance of :math:`r_c/2` in
either of the electron-ion distances and the largest meaningful
electron-electron distance is :math:`r_c`. This is the highest-order
Jastrow correlation function currently implemented.

Today, solid state applications of QMCPACK usually utilize one and
two-body B-spline Jastrow functions, with calculations on heavier
elements often also using the three-body term described above.

Example use case
^^^^^^^^^^^^^^^^

Here is an example of H2O molecule. After optimizing one and two body Jastrow factors, add the following block in the wavefunction.
The coefficients will be filled zero automatically if not given.

::

  <jastrow name="J3" type="eeI" function="polynomial" source="ion0" print="yes">
    <correlation ispecies="O" especies="u" isize="3" esize="3" rcut="10">
      <coefficients id="uuO" type="Array" optimize="yes"> </coefficients>
    </correlation>
    <correlation ispecies="O" especies1="u" especies2="d" isize="3" esize="3" rcut="10">
      <coefficients id="udO" type="Array" optimize="yes"> </coefficients>
    </correlation>
    <correlation ispecies="H" especies="u" isize="3" esize="3" rcut="10">
      <coefficients id="uuH" type="Array" optimize="yes"> </coefficients>
    </correlation>
    <correlation ispecies="H" especies1="u" especies2="d" isize="3" esize="3" rcut="10">
      <coefficients id="udH" type="Array" optimize="yes"> </coefficients>
    </correlation>
  </jastrow>

.. _ionwf:

Gaussian Product Wavefunction
-----------------------------

The Gaussian Product wavefunction implements :eq:`eq27`

.. math::
  :label: eq27

  \Psi(\vec{R}) = \prod_{i=1}^N \exp\left[ -\frac{(\vec{R}_i-\vec{R}_i^o)^2}{2\sigma_i^2} \right]

where :math:`\vec{R}_i` is the position of the :math:`i^{\text{th}}`
quantum particle and :math:`\vec{R}_i^o` is its center. :math:`\sigma_i`
is the width of the Gaussian orbital around center :math:`i`.

This variational wavefunction enhances single-particle density at chosen
spatial locations with adjustable strengths. It is useful whenever such
localization is physically relevant yet not captured by other parts of
the trial wavefunction. For example, in an electron-ion simulation of a
solid, the ions are localized around their crystal lattice sites. This
single-particle localization is not captured by the ion-ion Jastrow.
Therefore, the addition of this localization term will improve the
wavefunction. The simplest use case of this wavefunction is perhaps the
quantum harmonic oscillator (please see the “tests/models/sho” folder
for examples).

.. centered:: Input specification


Gaussian Product Wavefunction (ionwf):

  +----------+--------------+------------+-------------+-----------------------------------+
  | **Name** | **Datatype** | **Values** | **Default** | **Description**                   |
  +==========+==============+============+=============+===================================+
  | Name     | Text         | ionwf      | (Required)  | Unique name for this wavefunction |
  +----------+--------------+------------+-------------+-----------------------------------+
  | Width    | Floats       | 1.0 -1     | (Required)  | Widths of Gaussian orbitals       |
  +----------+--------------+------------+-------------+-----------------------------------+
  | Source   | Text         | ion0       | (Required)  | Name of classical particle set    |
  +----------+--------------+------------+-------------+-----------------------------------+

Additional information:

-  ``width`` There must be one width provided for each quantum particle.
   If a negative width is given, then its corresponding Gaussian orbital
   is removed. Negative width is useful if one wants to use Gaussian
   wavefunction for a subset of the quantum particles.

-  ``source`` The Gaussian centers must be specified in the form of a
   classical particle set. This classical particle set is likely the ion
   positions “ion0,” hence the name “ionwf.” However, arbitrary centers
   can be defined using a different particle set. Please refer to the
   examples in “tests/models/sho.”

Example Use Case
~~~~~~~~~~~~~~~~

::

  <qmcsystem>
    <simulationcell>
      <parameter name="bconds">
            n n n
      </parameter>
    </simulationcell>
    <particleset name="e">
      <group name="u" size="1">
        <parameter name="mass">5.0</parameter>
        <attrib name="position" datatype="posArray" condition="0">
          0.0001 -0.0001 0.0002
        </attrib>
      </group>
    </particleset>
    <particleset name="ion0" size="1">
      <group name="H">
        <attrib name="position" datatype="posArray" condition="0">
          0 0 0
        </attrib>
      </group>
    </particleset>
    <wavefunction target="e" id="psi0">
      <ionwf name="iwf" source="ion0" width="0.8165"/>
    </wavefunction>
    <hamiltonian name="h0" type="generic" target="e">
      <extpot type="HarmonicExt" mass="5.0" energy="0.3"/>
      <estimator type="latticedeviation" name="latdev"
        target="e"    tgroup="u"
        source="ion0" sgroup="H"/>
    </hamiltonian>
  </qmcsystem>

.. bibliography:: /bibs/intro_wavefunction.bib
