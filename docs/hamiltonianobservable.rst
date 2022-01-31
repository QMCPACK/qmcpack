.. _hamiltobs:

Hamiltonian and Observables
===========================

QMCPACK is capable of the simultaneous measurement of the Hamiltonian and many other quantum operators.  The Hamiltonian attains a special status among the available operators (also referred to as observables) because it ultimately generates all available information regarding the quantum system.  This is evident from an algorithmic standpoint as well since the Hamiltonian (embodied in the projector) generates the imaginary time dynamics of the walkers in DMC and reptation Monte Carlo (RMC).

This section covers how the Hamiltonian can be specified, component by component, by the user in the XML format native to \qmcpack. It also covers the input structure of statistical estimators corresponding to quantum observables such as the density, static structure factor, and forces.

The Hamiltonian
---------------

The many-body Hamiltonian in Hartree units is given by

.. math::
  :label: eq28

  \hat{H} = -\sum_i\frac{1}{2m_i}\nabla_i^2 + \sum_iv^{ext}(r_i) + \sum_{i<j}v^{qq}(r_i,r_j)   + \sum_{i\ell}v^{qc}(r_i,r_\ell)   + \sum_{\ell<m}v^{cc}(r_\ell,r_m)\:.

Here, the sums indexed by :math:`i/j` are over quantum particles, while
:math:`\ell/m` are reserved for classical particles. Often the quantum
particles are electrons, and the classical particles are ions, though is
not limited in this way. The mass of each quantum particle is denoted
:math:`m_i`, :math:`v^{qq}/v^{qc}/v^{cc}` are pair potentials between
quantum-quantum/quantum-classical/classical-classical particles, and
:math:`v^{ext}` denotes a purely external potential.

QMCPACK is designed modularly so that any potential can be supported with
minimal additions to the code base. Potentials currently supported
include Coulomb interactions in open and periodic boundary conditions,
the MPC potential, nonlocal pseudopotentials, helium pair potentials,
and various model potentials such as hard sphere, Gaussian, and modified
Poschl-Teller.

Reference information and examples for the ``<hamiltonian/>`` XML
element are provided subsequently. Detailed descriptions of the input
for individual potentials is given in the sections that follow.

``hamiltonian`` element:

  +------------------+----------------------------------------------------+
  | parent elements: | ``simulation, qmcsystem``                          |
  +------------------+----------------------------------------------------+
  | child elements:  | ``pairpot extpot estimator constant`` (deprecated) |
  +------------------+----------------------------------------------------+

attributes:

  +------------------------+--------------+----------------------+-------------+------------------------------------------+
  | **Name**               | **Datatype** | **Values**           | **Default** | **Description**                          |
  +========================+==============+======================+=============+==========================================+
  | ``name/id``:math:`^o`  | text         | *anything*           | h0          | Unique id for this Hamiltonian instance  |
  +------------------------+--------------+----------------------+-------------+------------------------------------------+
  | ``type``:math:`^o`     | text         |                      | generic     | *No current function*                    |
  +------------------------+--------------+----------------------+-------------+------------------------------------------+
  | ``role``:math:`^o`     | text         | primary/extra        | extra       | Designate as Hamiltonian or not          |
  +------------------------+--------------+----------------------+-------------+------------------------------------------+
  | ``source``:math:`^o`   | text         | ``particleset.name`` | i           | Identify classical ``particleset``       |
  +------------------------+--------------+----------------------+-------------+------------------------------------------+
  | ``target``:math:`^o`   | text         | ``particleset.name`` | e           | Identify quantum ``particleset``         |
  +------------------------+--------------+----------------------+-------------+------------------------------------------+
  | ``default``:math:`^o`  | boolean      | yes/no               | yes         | Include kinetic energy term implicitly   |
  +------------------------+--------------+----------------------+-------------+------------------------------------------+

Additional information:

-  **target:** Must be set to the name of the quantum ``particleset``.
   The default value is typically sufficient. In normal usage, no other
   attributes are provided.

.. code-block::
  :caption: All electron Hamiltonian XML element.
  :name: Listing 14

  <hamiltonian target="e">
    <pairpot name="ElecElec" type="coulomb" source="e" target="e"/>
    <pairpot name="ElecIon"  type="coulomb" source="i" target="e"/>
    <pairpot name="IonIon"   type="coulomb" source="i" target="i"/>
  </hamiltonian>


.. code-block::
  :caption: Pseudopotential Hamiltonian XML element.
  :name: Listing 15

  <hamiltonian target="e">
    <pairpot name="ElecElec"  type="coulomb" source="e" target="e"/>
    <pairpot name="PseudoPot" type="pseudo"  source="i" wavefunction="psi0" format="xml">
      <pseudo elementType="Li" href="Li.xml"/>
      <pseudo elementType="H" href="H.xml"/>
    </pairpot>
    <pairpot name="IonIon"    type="coulomb" source="i" target="i"/>
  </hamiltonian>

Pair potentials
---------------

Many pair potentials are supported.  Though only the most commonly used pair potentials are covered in detail in this section, all currently available potentials are listed subsequently.  If a potential you desire is not listed, or is not present at all, feel free to contact the developers.

``pairpot`` factory element:

  +------------------+--------------------+
  | parent elements: | ``hamiltonian``    |
  +------------------+--------------------+
  | child elements:  | ``type`` attribute |
  +------------------+--------------------+

  +------------------+---------+-----------------------------------------------+
  | **type options** | coulomb | Coulomb/Ewald potential                       |
  +------------------+---------+-----------------------------------------------+
  |                  | pseudo  | Semilocal pseudopotential                     |
  +------------------+---------+-----------------------------------------------+
  |                  | mpc     | Model periodic Coulomb interaction/correction |
  +------------------+---------+-----------------------------------------------+
  |                  | cpp     | Core polarization potential                   |
  +------------------+---------+-----------------------------------------------+
  |                  | skpot   | *Unknown*                                     |
  +------------------+---------+-----------------------------------------------+

shared attributes:

  +-----------------------+--------------+----------------------+------------------------+---------------------------------+
  | **Name**              | **Datatype** | **Values**           | **Default**            | **Description**                 |
  +=======================+==============+======================+========================+=================================+
  | ``type``:math:`^r`    | text         | *See above*          | 0                      | Select pairpot type             |
  +-----------------------+--------------+----------------------+------------------------+---------------------------------+
  | ``name``:math:`^r`    | text         | *Anything*           | any                    | Unique name for this pairpot    |
  +-----------------------+--------------+----------------------+------------------------+---------------------------------+
  | ``source``:math:`^r`  | text         | ``particleset.name`` | ``hamiltonian.target`` | Identify interacting particles  |
  +-----------------------+--------------+----------------------+------------------------+---------------------------------+
  | ``target``:math:`^r`  | text         | ``particleset.name`` | ``hamiltonian.target`` | Identify interacting particles  |
  +-----------------------+--------------+----------------------+------------------------+---------------------------------+
  | ``units``:math:`^o`   | text         |                      | hartree                | *No current function*           |
  +-----------------------+--------------+----------------------+------------------------+---------------------------------+

Additional information:

-  **type:** Used to select the desired pair potential. Must be selected
   from the list of type options.

-  **name:** A unique name used to identify this pair potential. Block
   averaged output data will appear under this name in ``scalar.dat``
   and/or ``stat.h5`` files.

-  **source/target:** These specify the particles involved in a pair
   interaction. If an interaction is between classical (e.g., ions) and
   quantum (e.g., electrons), ``source``/``target`` should be the name
   of the classical/quantum ``particleset``.

-  Only ``Coulomb, pseudo``, and ``mpc`` are described in detail in the
   following subsections. The older or less-used types (``cpp, skpot``)
   are not covered.

-  Available only if ``QMC_CUDA`` is not defined: ``skpot``.

-  Available only if ``OHMMS_DIM==3``: ``mpc, vhxc, pseudo``.

-  Available only if ``OHMMS_DIM==3`` and ``QMC_CUDA`` is not defined:
   ``cpp``.

Coulomb potentials
~~~~~~~~~~~~~~~~~~

The bare Coulomb potential is used in open boundary conditions:

.. math::
  :label: eq29

  V_c^{open} = \sum_{i<j}\frac{q_iq_j}{\left|{r_i-r_j}\right|}\:.

When periodic boundary conditions are selected, Ewald summation is used automatically:

.. math::
  :label: eq30

  V_c^{pbc} = \sum_{i<j}\frac{q_iq_j}{\left|{r_i-r_j}\right|} + \frac{1}{2}\sum_{L\ne0}\sum_{i,j}\frac{q_iq_j}{\left|{r_i-r_j+L}\right|}\:.

The sum indexed by :math:`L` is over all nonzero simulation cell lattice vectors.  In practice, the Ewald sum is broken into short- and long-range parts in a manner optimized for efficiency (see :cite:`Natoli1995`) for details.

For information on how to set the boundary conditions, consult :ref:`simulationcell`.

``pairpot type=coulomb`` element:

  +------------------+-----------------+
  | parent elements: | ``hamiltonian`` |
  +------------------+-----------------+
  | child elements:  | *None*          |
  +------------------+-----------------+

attributes:

  +-------------------------+--------------+----------------------+------------------------+---------------------------------+
  | **Name**                | **Datatype** | **Values**           | **Default**            | **Description**                 |
  +=========================+==============+======================+========================+=================================+
  | ``type``:math:`^r`      | text         | **coulomb**          |                        | Must be coulomb                 |
  +-------------------------+--------------+----------------------+------------------------+---------------------------------+
  | ``name/id``:math:`^r`   | text         | *anything*           | ElecElec               | Unique name for interaction     |
  +-------------------------+--------------+----------------------+------------------------+---------------------------------+
  | ``source``:math:`^r`    | text         | ``particleset.name`` | ``hamiltonian.target`` | Identify interacting particles  |
  +-------------------------+--------------+----------------------+------------------------+---------------------------------+
  | ``target``:math:`^r`    | text         | ``particleset.name`` | ``hamiltonian.target`` | Identify interacting particles  |
  +-------------------------+--------------+----------------------+------------------------+---------------------------------+
  | ``pbc``:math:`^o`       | boolean      | yes/no               | yes                    | Use Ewald summation             |
  +-------------------------+--------------+----------------------+------------------------+---------------------------------+
  | ``physical``:math:`^o`  | boolean      | yes/no               | yes                    | Hamiltonian(yes)/Observable(no) |
  +-------------------------+--------------+----------------------+------------------------+---------------------------------+
  | ``forces``              | boolean      | yes/no               | no                     | *Deprecated*                    |
  +-------------------------+--------------+----------------------+------------------------+---------------------------------+

Additional information:

-  **type/source/target:** See description for the previous generic
   ``pairpot`` factory element.

-  **name:** Traditional user-specified names for electron-electron,
   electron-ion, and ion-ion terms are ``ElecElec``, ``ElecIon``, and
   ``IonIon``, respectively. Although any choice can be used, the data
   analysis tools expect to find columns in ``*.scalar.dat`` with these
   names.

-  **pbc**: Ewald summation will not be performed if
   ``simulationcell.bconds== n n n``, regardless of the value of
   ``pbc``. Similarly, the ``pbc`` attribute can only be used to turn
   off Ewald summation if ``simulationcell.bconds!= n n n``. The default
   value is recommended.

-  **physical**: If ``physical==yes``, this pair potential is included
   in the Hamiltonian and will factor into the ``LocalEnergy`` reported
   by QMCPACK and also in the DMC branching weight. If ``physical==no``,
   then the pair potential is treated as a passive observable but not as
   part of the Hamiltonian itself. As such it does not contribute to the
   outputted ``LocalEnergy``. Regardless of the value of ``physical``
   output data will appear in ``scalar.dat`` in a column headed by
   ``name``.

.. code-block::
  :caption: QMCPXML element for Coulomb interaction between electrons.
  :name: Listing 16

  <pairpot name="ElecElec" type="coulomb" source="e" target="e"/>

.. code-block::
  :caption: QMCPXML element for Coulomb interaction between electrons and ions (all-electron only).
  :name: Listing 17

  <pairpot name="ElecIon"  type="coulomb" source="i" target="e"/>

.. code-block::
  :caption: QMCPXML element for Coulomb interaction between ions.
  :name: Listing 18

  <pairpot name="IonIon"   type="coulomb" source="i" target="i"/>

.. _nlpp:

Pseudopotentials
~~~~~~~~~~~~~~~~

QMCPACK supports pseudopotentials in semilocal form, which is local in the
radial coordinate and nonlocal in angular coordinates. When all angular
momentum channels above a certain threshold (:math:`\ell_{max}`) are
well approximated by the same potential
(:math:`V_{\bar{\ell}}\equiv V_{loc}`), the pseudopotential separates
into a fully local channel and an angularly nonlocal component:

.. math::
  :label: eq31

  V^{PP} = \sum_{ij}\Big(V_{\bar{\ell}}(\left|{r_i-\tilde{r}_j}\right|) + \sum_{\ell\ne\bar{\ell}}^{\ell_{max}}\sum_{m=-\ell}^\ell |{Y_{\ell m}}\rangle{\big[V_\ell(\left|{r_i-\tilde{r}_j}\right|) - V_{\bar{\ell}}(\left|{r_i-\tilde{r}_j}\right|) \big]}\langle{Y_{\ell m}}| \Big)\:.

Here the electron/ion index is :math:`i/j`, and only one type of ion is
shown for simplicity.

Evaluation of the localized pseudopotential energy
:math:`\Psi_T^{-1}V^{PP}\Psi_T` requires additional angular integrals.
These integrals are evaluated on a randomly shifted angular grid. The
size of this grid is determined by :math:`\ell_{max}`. See
:cite:`Mitas1991` for further detail.

uses the FSAtom pseudopotential file format associated with the “Free
Software Project for Atomic-scale Simulations” initiated in 2002. See
http://www.tddft.org/fsatom/manifest.php for more information. The
FSAtom format uses XML for structured data. Files in this format do not
use a specific identifying file extension; instead they are simply
suffixed with “``.xml``.” The tabular data format of CASINO is also
supported.

In addition to the semilocal pseudopotential above, spin-orbit 
interactions can also be included through the use of spin-orbit
pseudopotentials. The spin-orbit contribution can be written as

.. math::
  :label: eqn32

  V^{\rm SO} = \sum_{ij} \left(\sum_{\ell = 1}^{\ell_{max}-1} \frac{2}{2\ell+1} V^{\rm SO}_\ell \left( \left|r_i - \tilde{r}_j \right| \right) \sum_{m,m'=-\ell}^{\ell} | Y_{\ell m} \rangle  \langle Y_{\ell m} | \vec{\ell} \cdot \vec{s} | Y_{\ell m'}\rangle\langle Y_{\ell m'}|\right)\:.

Here, :math:`\vec{s}` is the spin operator. For each atom with a spin-orbit contribution,
the radial functions :math:`V_{\ell}^{\rm SO}` can be included in the pseudopotential 
“``.xml``” file.

``pairpot type=pseudo`` element:

  +------------------+-----------------+
  | parent elements: | ``hamiltonian`` |
  +------------------+-----------------+
  | child elements:  | ``pseudo``      |
  +------------------+-----------------+

attributes:

  +-----------------------------+--------------+-----------------------+------------------------+--------------------------------------------------+
  | **Name**                    | **Datatype** | **Values**            | **Default**            | **Description**                                  |
  +=============================+==============+=======================+========================+==================================================+
  | ``type``:math:`^r`          | text         | **pseudo**            |                        | Must be pseudo                                   |
  +-----------------------------+--------------+-----------------------+------------------------+--------------------------------------------------+
  | ``name/id``:math:`^r`       | text         | *anything*            | PseudoPot              | *No current function*                            |
  +-----------------------------+--------------+-----------------------+------------------------+--------------------------------------------------+
  | ``source``:math:`^r`        | text         | ``particleset.name``  | i                      | Ion ``particleset`` name                         |
  +-----------------------------+--------------+-----------------------+------------------------+--------------------------------------------------+
  | ``target``:math:`^r`        | text         | ``particleset.name``  | ``hamiltonian.target`` | Electron ``particleset`` name                    |
  +-----------------------------+--------------+-----------------------+------------------------+--------------------------------------------------+
  | ``pbc``:math:`^o`           | boolean      | yes/no                | yes*                   | Use Ewald summation                              |
  +-----------------------------+--------------+-----------------------+------------------------+--------------------------------------------------+
  | ``forces``                  | boolean      | yes/no                | no                     | *Deprecated*                                     |
  +-----------------------------+--------------+-----------------------+------------------------+--------------------------------------------------+
  | ``wavefunction``:math:`^r`  | text         | ``wavefunction.name`` | invalid                | Identify wavefunction                            |
  +-----------------------------+--------------+-----------------------+------------------------+--------------------------------------------------+
  | ``format``:math:`^r`        | text         | xml/table             | table                  | Select file format                               |
  +-----------------------------+--------------+-----------------------+------------------------+--------------------------------------------------+
  | ``algorithm``:math:`^o`     | text         | batched/non-batched   | batched                | Choose NLPP algorithm                            |
  +-----------------------------+--------------+-----------------------+------------------------+--------------------------------------------------+
  | ``DLA``:math:`^o`           | text         | yes/no                | no                     | Use determinant localization approximation       |
  +-----------------------------+--------------+-----------------------+------------------------+--------------------------------------------------+
  | ``physicalSO``:math:`^o`    | boolean      | yes/no                | yes                    | Include the SO contribution in the local energy  |
  +-----------------------------+--------------+-----------------------+------------------------+--------------------------------------------------+

Additional information:

-  **type/source/target** See description for the generic ``pairpot``
   factory element.

-  **name:** Ignored. Instead, default names will be present in
   ``*scalar.dat`` output files when pseudopotentials are used. The
   field ``LocalECP`` refers to the local part of the pseudopotential.
   If nonlocal channels are present, a ``NonLocalECP`` field will be
   added that contains the nonlocal energy summed over all angular
   momentum channels.

-  **pbc:** Ewald summation will not be performed if
   ``simulationcell.bconds== n n n``, regardless of the value of
   ``pbc``. Similarly, the ``pbc`` attribute can only be used to turn
   off Ewald summation if ``simulationcell.bconds!= n n n``.

-  **format:** If ``format``\ ==table, QMCPACK looks for ``*.psf`` files
   containing pseudopotential data in a tabular format. The files must
   be named after the ionic species provided in ``particleset`` (e.g.,
   ``Li.psf`` and ``H.psf``). If ``format``\ ==xml, additional
   ``pseudo`` child XML elements must be provided (see the following).
   These elements specify individual file names and formats (both the
   FSAtom XML and CASINO tabular data formats are supported).

-  **algorithm** The ``non-batched`` algorithm evaluates the ratios of
   wavefunction components together for each quadrature point and then
   one point after another. The ``batched`` algorithm evaluates the ratios
   of quadrature points together for each wavefunction component and
   then one component after another. Internally, it uses
   ``VirtualParticleSet`` for quadrature points. Hybrid orbital
   representation has an extra optimization enabled when using the
   batched algorithm. When OpenMP offload build is enabled, the default
   value is ``batched``. Otherwise, ``non-batched`` is the default.

-  **DLA** Determinant localization approximation
   (DLA) :cite:`Zen2019DLA` uses only the fermionic part of
   the wavefunction when calculating NLPP.

-  **physicalSO** If the spin-orbit components are included in the 
   ``.xml`` file, this flag allows control over whether the SO contribution
   is included in the local energy. 

.. code-block::
  :caption: QMCPXML element for pseudopotential electron-ion interaction (psf files).
  :name: Listing 19

    <pairpot name="PseudoPot" type="pseudo"  source="i" wavefunction="psi0" format="psf"/>

.. code-block::
  :caption: QMCPXML element for pseudopotential electron-ion interaction (xml files). If SOC terms present in xml, they are added to local energy
  :name: Listing 20

    <pairpot name="PseudoPot" type="pseudo"  source="i" wavefunction="psi0" format="xml">
      <pseudo elementType="Li" href="Li.xml"/>
      <pseudo elementType="H" href="H.xml"/>
    </pairpot>

.. code-block::
  :caption: QMCPXML element for pseudopotential to accumulate the spin-orbit energy, but do not include in local energy
  :name: Listing 21
  
    <pairpot name="PseudoPot" type="pseudo" source="i" wavefunction="psi0" format="xml" physicalSO="no">
      <pseudo elementType="Pb" href="Pb.xml"/>
    </pairpot>
Details of ``<pseudo/>`` input elements are shown in the following. It
is possible to include (or construct) a full pseudopotential directly in
the input file without providing an external file via ``href``. The full
XML format for pseudopotentials is not yet covered.

``pseudo`` element:

  +------------------+-----------------------------+
  | parent elements: | ``pairpot type=pseudo``     |
  +------------------+-----------------------------+
  | child elements:  | ``header local grid``       |
  +------------------+-----------------------------+

attributes:

  +-----------------------------------+--------------+-----------------+-------------+---------------------------+
  | **Name**                          | **Datatype** | **Values**      | **Default** | **Description**           |
  +===================================+==============+=================+=============+===========================+
  | ``elementType/symbol``:math:`^r`  | text         | ``groupe.name`` | none        | Identify ionic species    |
  +-----------------------------------+--------------+-----------------+-------------+---------------------------+
  | ``href``:math:`^r`                | text         | *filepath*      | none        | Pseudopotential file path |
  +-----------------------------------+--------------+-----------------+-------------+---------------------------+
  | ``format``:math:`^r`              | text         | xml/casino      | xml         | Specify file format       |
  +-----------------------------------+--------------+-----------------+-------------+---------------------------+
  | ``cutoff``:math:`^o`              | real         |                 |             | Nonlocal cutoff radius    |
  +-----------------------------------+--------------+-----------------+-------------+---------------------------+
  | ``lmax``:math:`^o`                | integer      |                 |             | Largest angular momentum  |
  +-----------------------------------+--------------+-----------------+-------------+---------------------------+
  | ``nrule``:math:`^o`               | integer      |                 |             | Integration grid order    |
  +-----------------------------------+--------------+-----------------+-------------+---------------------------+

.. code-block::
  :caption: QMCPXML element for pseudopotential of single ionic species.
  :name: Listing 21b

    <pseudo elementType="Li" href="Li.xml"/>

MPC Interaction/correction
~~~~~~~~~~~~~~~~~~~~~~~~~~

The MPC interaction is an alternative to direct Ewald summation. The MPC
corrects the exchange correlation hole to more closely match its
thermodynamic limit. Because of this, the MPC exhibits smaller
finite-size errors than the bare Ewald interaction, though a few
alternative and competitive finite-size correction schemes now exist.
The MPC is itself often used just as a finite-size correction in
post-processing (set ``physical=false`` in the input).

``pairpot type=mpc`` element:

  +------------------+-----------------+
  | parent elements: | ``hamiltonian`` |
  +------------------+-----------------+
  | child elements:  | *None*          |
  +------------------+-----------------+

attributes:

  +-------------------------+--------------+----------------------+------------------------+---------------------------------+
  | **Name**                | **Datatype** | **Values**           | **Default**            | **Description**                 |
  +=========================+==============+======================+========================+=================================+
  | ``type``:math:`^r`      | text         | **mpc**              |                        | Must be MPC                     |
  +-------------------------+--------------+----------------------+------------------------+---------------------------------+
  | ``name/id``:math:`^r`   | text         | *anything*           | MPC                    | Unique name for interaction     |
  +-------------------------+--------------+----------------------+------------------------+---------------------------------+
  | ``source``:math:`^r`    | text         | ``particleset.name`` | ``hamiltonian.target`` | Identify interacting particles  |
  +-------------------------+--------------+----------------------+------------------------+---------------------------------+
  | ``target``:math:`^r`    | text         | ``particleset.name`` | ``hamiltonian.target`` | Identify interacting particles  |
  +-------------------------+--------------+----------------------+------------------------+---------------------------------+
  | ``physical``:math:`^o`  | boolean      | yes/no               | no                     | Hamiltonian(yes)/observable(no) |
  +-------------------------+--------------+----------------------+------------------------+---------------------------------+
  | ``cutoff``              | real         | :math:`>0`           | 30.0                   | Kinetic energy cutoff           |
  +-------------------------+--------------+----------------------+------------------------+---------------------------------+

Remarks:

-  ``physical``: Typically set to ``no``, meaning the standard Ewald
   interaction will be used during sampling and MPC will be measured as
   an observable for finite-size post-correction. If ``physical`` is
   ``yes``, the MPC interaction will be used during sampling. In this
   case an electron-electron Coulomb ``pairpot`` element should not be
   supplied.

-  **Developer note:** Currently the ``name`` attribute for the MPC
   interaction is ignored. The name is always reset to ``MPC``.

.. code-block::
  :caption: MPC for finite-size postcorrection.
  :name: Listing 22

    <pairpot type="MPC" name="MPC" source="e" target="e" ecut="60.0" physical="no"/>

General estimators
------------------

A broad range of estimators for physical observables are available in QMCPACK.
The following sections contain input details for the total number
density (``density``), number density resolved by particle spin
(``spindensity``), spherically averaged pair correlation function
(``gofr``), static structure factor (``sk``), static structure factor
(``skall``), energy density (``energydensity``), one body reduced
density matrix (``dm1b``), :math:`S(k)` based kinetic energy correction
(``chiesa``), forward walking (``ForwardWalking``), and force
(``Force``) estimators. Other estimators are not yet covered.

When an ``<estimator/>`` element appears in ``<hamiltonian/>``, it is
evaluated for all applicable chained QMC runs (e.g.,
VMC\ :math:`\rightarrow`\ DMC\ :math:`\rightarrow`\ DMC). Estimators are
generally not accumulated during wavefunction optimization sections. If
an ``<estimator/>`` element is instead provided in a particular
``<qmc/>`` element, that estimator is only evaluated for that specific
section (e.g., during VMC only).

``estimator`` factory element:

  +------------------+----------------------+
  | parent elements: | ``hamiltonian, qmc`` |
  +------------------+----------------------+
  | type selector:   | ``type`` attribute   |
  +------------------+----------------------+

  +------------------+------------------+-----------------------------------------------------------+
  | **type options** | density          | Density on a grid                                         |
  +------------------+------------------+-----------------------------------------------------------+
  |                  | spindensity      | Spin density on a grid                                    |
  +------------------+------------------+-----------------------------------------------------------+
  |                  | gofr             | Pair correlation function (quantum species)               |
  +------------------+------------------+-----------------------------------------------------------+
  |                  | sk               | Static structure factor                                   |
  +------------------+------------------+-----------------------------------------------------------+
  |                  | SkAll            | Static structure factor needed for finite size correction |
  +------------------+------------------+-----------------------------------------------------------+
  |                  | structurefactor  | Species resolved structure factor                         |
  +------------------+------------------+-----------------------------------------------------------+
  |                  | species kinetic  | Species resolved kinetic energy                           |
  +------------------+------------------+-----------------------------------------------------------+
  |                  | latticedeviation | Spatial deviation between two particlesets                |
  +------------------+------------------+-----------------------------------------------------------+
  |                  | momentum         | Momentum distribution                                     |
  +------------------+------------------+-----------------------------------------------------------+
  |                  | energydensity    | Energy density on uniform or Voronoi grid                 |
  +------------------+------------------+-----------------------------------------------------------+
  |                  | dm1b             | One body density matrix in arbitrary basis                |
  +------------------+------------------+-----------------------------------------------------------+
  |                  | chiesa           | Chiesa-Ceperley-Martin-Holzmann kinetic energy correction |
  +------------------+------------------+-----------------------------------------------------------+
  |                  | Force            | Family of "force" estimators (see :ref:`force-est`)       |
  +------------------+------------------+-----------------------------------------------------------+
  |                  | ForwardWalking   | Forward walking values for existing estimators            |
  +------------------+------------------+-----------------------------------------------------------+
  |                  | orbitalimages    | Create image files for orbitals, then exit                |
  +------------------+------------------+-----------------------------------------------------------+
  |                  | flux             | Checks sampling of kinetic energy                         |
  +------------------+------------------+-----------------------------------------------------------+
  |                  | localmoment      | Atomic spin polarization within cutoff radius             |
  +------------------+------------------+-----------------------------------------------------------+
  |                  | Pressure         | *No current function*                                     |
  +------------------+------------------+-----------------------------------------------------------+

shared attributes:

  +---------------------+--------------+-------------+-------------+--------------------------------+
  | **Name**            | **Datatype** | **Values**  | **Default** | **Description**                |
  +=====================+==============+=============+=============+================================+
  | ``type``:math:`^r`  | text         | *See above* | 0           | Select estimator type          |
  +---------------------+--------------+-------------+-------------+--------------------------------+
  | ``name``:math:`^r`  | text         | *anything*  | any         | Unique name for this estimator |
  +---------------------+--------------+-------------+-------------+--------------------------------+

Chiesa-Ceperley-Martin-Holzmann kinetic energy correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This estimator calculates a finite-size correction to the kinetic energy following the formalism laid out in :cite:`Chiesa2006`.  The total energy can be corrected for finite-size effects by using this estimator in conjunction with the MPC correction.

``estimator type=chiesa`` element:

  +------------------+----------------------+
  | parent elements: | ``hamiltonian, qmc`` |
  +------------------+----------------------+
  | child elements:  | *None*               |
  +------------------+----------------------+

attributes:

  +-----------------------+--------------+------------------------+-------------+----------------------------+
  | **Name**              | **Datatype** | **Values**             | **Default** | **Description**            |
  +=======================+==============+========================+=============+============================+
  | ``type``:math:`^r`    | text         | **chiesa**             |             | Must be chiesa             |
  +-----------------------+--------------+------------------------+-------------+----------------------------+
  | ``name``:math:`^o`    | text         | *anything*             | KEcorr      | Always reset to KEcorr     |
  +-----------------------+--------------+------------------------+-------------+----------------------------+
  | ``source``:math:`^o`  | text         | ``particleset.name``   | e           | Identify quantum particles |
  +-----------------------+--------------+------------------------+-------------+----------------------------+
  | ``psi``:math:`^o`     | text         | ``wavefunction.name``  | psi0        | Identify wavefunction      |
  +-----------------------+--------------+------------------------+-------------+----------------------------+

.. code-block::
  :caption: "Chiesa" kinetic energy finite-size postcorrection.
  :name: Listing 23

     <estimator name="KEcorr" type="chiesa" source="e" psi="psi0"/>

Density estimator
~~~~~~~~~~~~~~~~~

The particle number density operator is given by

.. math::
  :label: eq32

  \hat{n}_r = \sum_i\delta(r-r_i)\:.

The ``density`` estimator accumulates the number density on a uniform
histogram grid over the simulation cell. The value obtained for a grid
cell :math:`c` with volume :math:`\Omega_c` is then the average number
of particles in that cell:

.. math::
  :label: eq33

  n_c = \int dR \left|{\Psi}\right|^2 \int_{\Omega_c}dr \sum_i\delta(r-r_i)\:.

``estimator type=density`` element:

  +------------------+----------------------+
  | parent elements: | ``hamiltonian, qmc`` |
  +------------------+----------------------+
  | child elements:  | *None*               |
  +------------------+----------------------+

attributes:

  +--------------------------+---------------+-------------------------------------------+------------------------------------------+------------------------------------------+
  | **Name**                 | **Datatype**  | **Values**                                | **Default**                              | **Description**                          |
  +==========================+===============+===========================================+==========================================+==========================================+
  | ``type``:math:`^r`       | text          | **density**                               |                                          | Must be density                          |
  +--------------------------+---------------+-------------------------------------------+------------------------------------------+------------------------------------------+
  | ``name``:math:`^r`       | text          | *anything*                                | any                                      | Unique name for estimator                |
  +--------------------------+---------------+-------------------------------------------+------------------------------------------+------------------------------------------+
  | ``delta``:math:`^o`      | real array(3) | :math:`0\le v_i \le 1`                    | 0.1 0.1 0.1                              | Grid cell spacing, unit coords           |
  +--------------------------+---------------+-------------------------------------------+------------------------------------------+------------------------------------------+
  | ``x_min``:math:`^o`      | real          | :math:`>0`                                | 0                                        | Grid starting point in x (Bohr)          |
  +--------------------------+---------------+-------------------------------------------+------------------------------------------+------------------------------------------+
  | ``x_max``:math:`^o`      | real          | :math:`>0`                                | :math:`|` ``lattice[0]`` :math:`|`       | Grid ending point in x (Bohr)            |
  +--------------------------+---------------+-------------------------------------------+------------------------------------------+------------------------------------------+
  | ``y_min``:math:`^o`      | real          | :math:`>0`                                | 0                                        | Grid starting point in y (Bohr)          |
  +--------------------------+---------------+-------------------------------------------+------------------------------------------+------------------------------------------+
  | ``y_max``:math:`^o`      | real          | :math:`>0`                                | :math:`|` ``lattice[1]`` :math:`|`       | Grid ending point in y (Bohr)            |
  +--------------------------+---------------+-------------------------------------------+------------------------------------------+------------------------------------------+
  | ``z_min``:math:`^o`      | real          | :math:`>0`                                | 0                                        | Grid starting point in z (Bohr)          |
  +--------------------------+---------------+-------------------------------------------+------------------------------------------+------------------------------------------+
  | ``z_max``:math:`^o`      | real          | :math:`>0`                                | :math:`|` ``lattice[2]`` :math:`|`       | Grid ending point in z (Bohr)            |
  +--------------------------+---------------+-------------------------------------------+------------------------------------------+------------------------------------------+
  | ``potential``:math:`^o`  | boolean       | yes/no                                    | no                                       | Accumulate local potential, *deprecated* |
  +--------------------------+---------------+-------------------------------------------+------------------------------------------+------------------------------------------+
  | ``debug``:math:`^o`      | boolean       | yes/no                                    | no                                       | *No current function*                    |
  +--------------------------+---------------+-------------------------------------------+------------------------------------------+------------------------------------------+

Additional information:

-  ``name``: The name provided will be used as a label in the
   ``stat.h5`` file for the blocked output data. Postprocessing tools
   expect ``name="Density."``

-  ``delta``: This sets the histogram grid size used to accumulate the
   density:
   ``delta="0.1 0.1 0.05"``\ :math:`\rightarrow 10\times 10\times 20`
   grid,
   ``delta="0.01 0.01 0.01"``\ :math:`\rightarrow 100\times 100\times 100`
   grid. The density grid is written to a ``stat.h5`` file at the end of
   each MC block. If you request many :math:`blocks` in a ``<qmc/>``
   element, or select a large grid, the resulting ``stat.h5`` file could
   be many gigabytes in size.

-  ``*_min/*_max``: Can be used to select a subset of the simulation
   cell for the density histogram grid. For example if a (cubic)
   simulation cell is 20 Bohr on a side, setting ``*_min=5.0`` and
   ``*_max=15.0`` will result in a density histogram grid spanning a
   :math:`10\times 10\times 10` Bohr cube about the center of the box.
   Use of ``x_min, x_max, y_min, y_max, z_min, z_max`` is only
   appropriate for orthorhombic simulation cells with open boundary
   conditions.

-  When open boundary conditions are used, a ``<simulationcell/>``
   element must be explicitly provided as the first subelement of
   ``<qmcsystem/>`` for the density estimator to work. In this case the
   molecule should be centered around the middle of the simulation cell
   (:math:`L/2`) and not the origin (:math:`0` since the space within
   the cell, and hence the density grid, is defined from :math:`0` to
   :math:`L`).

.. code-block::
  :caption: QMCPXML,caption=Density estimator (uniform grid).
  :name: Listing 24

     <estimator name="Density" type="density" delta="0.05 0.05 0.05"/>

Spin density estimator
~~~~~~~~~~~~~~~~~~~~~~

The spin density is similar to the total density described previously.  In this case, the sum over particles is performed independently for each spin component.

``estimator type=spindensity`` element:

  +------------------+----------------------+
  | parent elements: | ``hamiltonian, qmc`` |
  +------------------+----------------------+
  | child elements:  | *None*               |
  +------------------+----------------------+

attributes:

  +-----------------------+--------------+-----------------+-------------+-------------------------------+
  | **Name**              | **Datatype** | **Values**      | **Default** | **Description**               |
  +=======================+==============+=================+=============+===============================+
  | ``type``:math:`^r`    | text         | **spindensity** |             | Must be spindensity           |
  +-----------------------+--------------+-----------------+-------------+-------------------------------+
  | ``name``:math:`^r`    | text         | *anything*      | any         | Unique name for estimator     |
  +-----------------------+--------------+-----------------+-------------+-------------------------------+
  | ``report``:math:`^o`  | boolean      | yes/no          | no          | Write setup details to stdout |
  +-----------------------+--------------+-----------------+-------------+-------------------------------+

parameters:

  +----------------------------+------------------+----------------------+-------------+----------------------------------+
  | **Name**                   | **Datatype**     | **Values**           | **Default** | **Description**                  |
  +============================+==================+======================+=============+==================================+
  | ``grid``:math:`^o`         | integer array(3) | :math:`v_i>`         |             | Grid cell count                  |
  +----------------------------+------------------+----------------------+-------------+----------------------------------+
  | ``dr``:math:`^o`           | real array(3)    | :math:`v_i>`         |             | Grid cell spacing (Bohr)         |
  +----------------------------+------------------+----------------------+-------------+----------------------------------+
  | ``cell``:math:`^o`         | real array(3,3)  | *anything*           |             | Volume grid exists in            |
  +----------------------------+------------------+----------------------+-------------+----------------------------------+
  | ``corner``:math:`^o`       | real array(3)    | *anything*           |             | Volume corner location           |
  +----------------------------+------------------+----------------------+-------------+----------------------------------+
  | ``center``:math:`^o`       | real array (3)   | *anything*           |             | Volume center/origin location    |
  +----------------------------+------------------+----------------------+-------------+----------------------------------+
  | ``voronoi``:math:`^o`      | text             | ``particleset.name`` |             | *Under development*              |
  +----------------------------+------------------+----------------------+-------------+----------------------------------+
  | ``test_moves``:math:`^o`   | integer          | :math:`>=0`          | 0           | Test estimator with random moves |
  +----------------------------+------------------+----------------------+-------------+----------------------------------+

Additional information:

-  ``name``: The name provided will be used as a label in the
   ``stat.h5`` file for the blocked output data. Postprocessing tools
   expect ``name="SpinDensity."``

-  ``grid``: The grid sets the dimension of the histogram grid. Input
   like ``<parameter name="grid"> 40 40 40 </parameter>`` requests a
   :math:`40 \times 40\times 40` grid. The shape of individual grid
   cells is commensurate with the supercell shape.

-  ``dr``: The ``dr`` sets the real-space dimensions of grid cell edges
   (Bohr units). Input like
   ``<parameter name="dr"> 0.5 0.5 0.5 </parameter>`` in a supercell
   with axes of length 10 Bohr each (but of arbitrary shape) will
   produce a :math:`20\times 20\times 20` grid. The inputted ``dr``
   values are rounded to produce an integer number of grid cells along
   each supercell axis. Either ``grid`` or ``dr`` must be provided, but
   not both.

-  ``cell``: When ``cell`` is provided, a user-defined grid volume is
   used instead of the global supercell. This must be provided if open
   boundary conditions are used. Additionally, if ``cell`` is provided,
   the user must specify where the volume is located in space in
   addition to its size/shape (``cell``) using either the ``corner`` or
   ``center`` parameters.

-  ``corner``: The grid volume is defined as
   :math:`corner+\sum_{d=1}^3u_dcell_d` with :math:`0<u_d<1` (“cell”
   refers to either the supercell or user-provided cell).

-  ``center``: The grid volume is defined as
   :math:`center+\sum_{d=1}^3u_dcell_d` with :math:`-1/2<u_d<1/2`
   (“cell” refers to either the supercell or user-provided cell).
   ``corner/center`` can be used to shift the grid even if ``cell`` is
   not specified. Simultaneous use of ``corner`` and ``center`` will
   cause QMCPACK to abort.

.. code-block::
  :caption: Spin density estimator (uniform grid).
  :name: Listing 25

  <estimator type="spindensity" name="SpinDensity" report="yes">
    <parameter name="grid"> 40 40 40 </parameter>
  </estimator>

.. code-block::
  :caption: Spin density estimator (uniform grid centered about origin).
  :name: Listing 26

  <estimator type="spindensity" name="SpinDensity" report="yes">
    <parameter name="grid">
      20 20 20
    </parameter>
    <parameter name="center">
      0.0 0.0 0.0
    </parameter>
    <parameter name="cell">
      10.0  0.0  0.0
       0.0 10.0  0.0
       0.0  0.0 10.0
    </parameter>
  </estimator>

Pair correlation function, :math:`g(r)`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The functional form of the species-resolved radial pair correlation function operator is

.. math::
  :label: eq34

  g_{ss'}(r) = \frac{V}{4\pi r^2N_sN_{s'}}\sum_{i_s=1}^{N_s}\sum_{j_{s'}=1}^{N_{s'}}\delta(r-|r_{i_s}-r_{j_{s'}}|)\:,

where :math:`N_s` is the number of particles of species :math:`s` and
:math:`V` is the supercell volume. If :math:`s=s'`, then the sum is
restricted so that :math:`i_s\ne j_s`.

In QMCPACK, an estimate of :math:`g_{ss'}(r)` is obtained as a radial
histogram with a set of :math:`N_b` uniform bins of width
:math:`\delta r`. This can be expressed analytically as

.. math::
  :label: eq35

  \tilde{g}_{ss'}(r) = \frac{V}{4\pi r^2N_sN_{s'}}\sum_{i=1}^{N_s}\sum_{j=1}^{N_{s'}}\frac{1}{\delta r}\int_{r-\delta r/2}^{r+\delta r/2}dr'\delta(r'-|r_{si}-r_{s'j}|)\:,

where the radial coordinate :math:`r` is restricted to reside at the bin
centers, :math:`\delta r/2, 3 \delta r/2, 5 \delta r/2, \ldots`.

``estimator type=gofr`` element:

  +------------------+----------------------+
  | parent elements: | ``hamiltonian, qmc`` |
  +------------------+----------------------+
  | child elements:  | *None*               |
  +------------------+----------------------+

attributes:

  +-------------------------------+--------------+----------------------+------------------------+-------------------------+
  | **Name**                      | **Datatype** | **Values**           | **Default**            | **Description**         |
  +===============================+==============+======================+========================+=========================+
  | ``type``:math:`^r`            | text         | **gofr**             |                        | Must be gofr            |
  +-------------------------------+--------------+----------------------+------------------------+-------------------------+
  | ``name``:math:`^o`            | text         | *anything*           | any                    | *No current function*   |
  +-------------------------------+--------------+----------------------+------------------------+-------------------------+
  | ``num_bin``:math:`^r`         | integer      | :math:`>1`           | 20                     | # of histogram bins     |
  +-------------------------------+--------------+----------------------+------------------------+-------------------------+
  | ``rmax``:math:`^o`            | real         | :math:`>0`           | 10                     | Histogram extent (Bohr) |
  +-------------------------------+--------------+----------------------+------------------------+-------------------------+
  | ``dr``:math:`^o`              | real         | :math:`0`            | 0.5                    | *No current function*   |
  +-------------------------------+--------------+----------------------+------------------------+-------------------------+
  | ``debug``:math:`^o`           | boolean      | yes/no               | no                     | *No current function*   |
  +-------------------------------+--------------+----------------------+------------------------+-------------------------+
  | ``target``:math:`^o`          | text         | ``particleset.name`` | ``hamiltonian.target`` | Quantum particles       |
  +-------------------------------+--------------+----------------------+------------------------+-------------------------+
  | ``source/sources``:math:`^o`  | text array   | ``particleset.name`` | ``hamiltonian.target`` | Classical particles     |
  +-------------------------------+--------------+----------------------+------------------------+-------------------------+

Additional information:

-  ``num_bin:`` This is the number of bins in each species pair radial
   histogram.

-  ``rmax:`` This is the maximum pair distance included in the
   histogram. The uniform bin width is
   :math:`\delta r=\texttt{rmax/num\_bin}`. If periodic boundary
   conditions are used for any dimension of the simulation cell, then
   the default value of ``rmax`` is the simulation cell radius instead
   of 10 Bohr. For open boundary conditions, the volume (:math:`V`) used
   is 1.0 Bohr\ :math:`^3`.

-  ``source/sources:`` If unspecified, only pair correlations between
   each species of quantum particle will be measured. For each classical
   particleset specified by ``source/sources``, additional pair
   correlations between each quantum and classical species will be
   measured. Typically there is only one classical particleset (e.g.,
   ``source="ion0"``), but there can be several in principle (e.g.,
   ``sources="ion0 ion1 ion2"``).

-  ``target:`` The default value is the preferred usage (i.e.,
   ``target`` does not need to be provided).

-  Data is output to the ``stat.h5`` for each QMC subrun. Individual
   histograms are named according to the quantum particleset and index
   of the pair. For example, if the quantum particleset is named “e" and
   there are two species (up and down electrons, say), then there will
   be three sets of histogram data in each ``stat.h5`` file named
   ``gofr_e_0_0``, ``gofr_e_0_1``, and ``gofr_e_1_1`` for up-up,
   up-down, and down-down correlations, respectively.

.. code-block::
  :caption: Pair correlation function estimator element.
  :name: Listing 27

  <estimator type="gofr" name="gofr" num_bin="200" rmax="3.0" />

.. code-block::
  :caption: Pair correlation function estimator element with additional electron-ion correlations.
  :name: Listing 28

  <estimator type="gofr" name="gofr" num_bin="200" rmax="3.0" source="ion0" />

Static structure factor, :math:`S(k)`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let
:math:`\rho^e_{\mathbf{k}}=\sum_j e^{i \mathbf{k}\cdot\mathbf{r}_j^e}`
be the Fourier space electron density, with :math:`\mathbf{r}^e_j` being
the coordinate of the j-th electron. :math:`\mathbf{k}` is a wavevector
commensurate with the simulation cell. QMCPACK allows the user to
accumulate the static electron structure factor :math:`S(\mathbf{k})` at
all commensurate :math:`\mathbf{k}` such that
:math:`|\mathbf{k}| \leq (LR\_DIM\_CUTOFF) r_c`. :math:`N^e` is the
number of electrons, ``LR_DIM_CUTOFF`` is the optimized breakup
parameter, and :math:`r_c` is the Wigner-Seitz radius. It is defined as
follows:

.. math::
  :label: eq36

  S(\mathbf{k}) = \frac{1}{N^e}\langle \rho^e_{-\mathbf{k}} \rho^e_{\mathbf{k}} \rangle\:.

``estimator type=sk`` element:

  +------------------+----------------------+
  | parent elements: | ``hamiltonian, qmc`` |
  +------------------+----------------------+
  | child elements:  | *None*               |
  +------------------+----------------------+

attributes:

  +---------------------+--------------+------------+-------------+-----------------------------------------------------+
  | **Name**            | **Datatype** | **Values** | **Default** | **Description**                                     |
  +=====================+==============+============+=============+=====================================================+
  | ``type``:math:`^r`  | text         | sk         |             | Must sk                                             |
  +---------------------+--------------+------------+-------------+-----------------------------------------------------+
  | ``name``:math:`^r`  | text         | *anything* | any         | Unique name for estimator                           |
  +---------------------+--------------+------------+-------------+-----------------------------------------------------+
  | ``hdf5``:math:`^o`  | boolean      | yes/no     | no          |  Output to ``stat.h5`` (yes) or ``scalar.dat`` (no) |
  +---------------------+--------------+------------+-------------+-----------------------------------------------------+

Additional information:

-  ``name:`` This is the unique name for estimator instance. A data
   structure of the same name will appear in ``stat.h5`` output files.

-  ``hdf5:`` If ``hdf5==yes``, output data for :math:`S(k)` is directed
   to the ``stat.h5`` file (recommended usage). If ``hdf5==no``, the
   data is instead routed to the ``scalar.dat`` file, resulting in many
   columns of data with headings prefixed by ``name`` and postfixed by
   the k-point index (e.g., ``sk_0 sk_1 …sk_1037 …``).

-  This estimator only works in periodic boundary conditions. Its
   presence in the input file is ignored otherwise.

-  This is not a species-resolved structure factor. Additionally, for
   :math:`\mathbf{k}` vectors commensurate with the unit cell,
   :math:`S(\mathbf{k})` will include contributions from the static
   electronic density, thus meaning it will not accurately measure the
   electron-electron density response.

.. code-block::
  :caption: Static structure factor estimator element.
  :name: Listing 29

    <estimator type="sk" name="sk" hdf5="yes"/>

Static structure factor, ``SkAll``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to compute the finite size correction to the potential energy,
records of :math:`\rho(\mathbf{k})` is required. What sets ``SkAll``
apart from ``sk`` is that ``SkAll`` records :math:`\rho(\mathbf{k})` in
addition to :math:`s(\mathbf{k})`.

``estimator type=SkAll`` element:

  +------------------+----------------------+
  | parent elements: | ``hamiltonian, qmc`` |
  +------------------+----------------------+
  | child elements:  | *None*               |
  +------------------+----------------------+

attributes:

  +----------------------------+--------------+---------------------------+-------------+-------------------------------------------------------------------------------------------------+
  | **Name**                   | **Datatype** | **Values**                | **Default** | **Description**                                                                                 |
  +============================+==============+===========================+=============+=================================================================================================+
  | ``type``:math:`^r`         | text         | sk                        |             | Must be sk                                                                                      |
  +----------------------------+--------------+---------------------------+-------------+-------------------------------------------------------------------------------------------------+
  | ``name``:math:`^r`         | text         | *anything*                | any         | Unique name for estimator                                                                       |
  +----------------------------+--------------+---------------------------+-------------+-------------------------------------------------------------------------------------------------+
  | ``source``:math:`^r`       | text         | Ion ParticleSet name      | None        | `-`                                                                                             |
  +----------------------------+--------------+---------------------------+-------------+-------------------------------------------------------------------------------------------------+
  | ``target``:math:`^r`       | text         | Electron ParticleSet name | None        | `-`                                                                                             |
  +----------------------------+--------------+---------------------------+-------------+-------------------------------------------------------------------------------------------------+
  | ``hdf5``:math:`^o`         | boolean      | yes/no                    | no          | Output to ``stat.h5`` (yes) or ``scalar.dat`` (no)                                              |
  +----------------------------+--------------+---------------------------+-------------+-------------------------------------------------------------------------------------------------+
  | ``writeionion``:math:`^o`  | boolean      | yes/no                    | no          | Writes file rhok_IonIon.dat containing :math:`s(\mathbf{k})` for the ions                       |
  +----------------------------+--------------+---------------------------+-------------+-------------------------------------------------------------------------------------------------+

Additional information:

-  ``name:`` This is the unique name for estimator instance. A data
   structure of the same name will appear in ``stat.h5`` output files.

-  ``hdf5:`` If ``hdf5==yes``, output data is directed to the
   ``stat.h5`` file (recommended usage). If ``hdf5==no``, the data is
   instead routed to the ``scalar.dat`` file, resulting in many columns
   of data with headings prefixed by ``rhok`` and postfixed by the
   k-point index.

-  This estimator only works in periodic boundary conditions. Its
   presence in the input file is ignored otherwise.

-  This is not a species-resolved structure factor. Additionally, for
   :math:`\mathbf{k}` vectors commensurate with the unit cell,
   :math:`S(\mathbf{k})` will include contributions from the static
   electronic density, thus meaning it wil not accurately measure the
   electron-electron density response.

.. code-block::
  :caption: SkAll estimator element.
  :name: Listing 30

    <estimator type="skall" name="SkAll" source="ion0" target="e" hdf5="yes"/>

Species kinetic energy
~~~~~~~~~~~~~~~~~~~~~~

Record species-resolved kinetic energy instead of the total kinetic
energy in the ``Kinetic`` column of scalar.dat. ``SpeciesKineticEnergy``
is arguably the simplest estimator in QMCPACK. The implementation of
this estimator is detailed in
``manual/estimator/estimator_implementation.pdf``.

``estimator type=specieskinetic`` element:

  +------------------+----------------------+
  | parent elements: | ``hamiltonian, qmc`` |
  +------------------+----------------------+
  | child elements:  | *None*               |
  +------------------+----------------------+

attributes:

  +---------------------+--------------+----------------+-------------+-----------------------------+
  | **Name**            | **Datatype** | **Values**     | **Default** | **Description**             |
  +=====================+==============+================+=============+=============================+
  | ``type``:math:`^r`  | text         | specieskinetic |             | Must be specieskinetic      |
  +---------------------+--------------+----------------+-------------+-----------------------------+
  | ``name``:math:`^r`  | text         | *anything*     | any         | Unique name for estimator   |
  +---------------------+--------------+----------------+-------------+-----------------------------+
  | ``hdf5``:math:`^o`  | boolean      | yes/no         | no          | Output to ``stat.h5`` (yes) |
  +---------------------+--------------+----------------+-------------+-----------------------------+

.. code-block::
  :caption: Species kinetic energy estimator element.
  :name: Listing 31

    <estimator type="specieskinetic" name="skinetic" hdf5="no"/>

Lattice deviation estimator
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Record deviation of a group of particles in one particle set (target) from a group of particles in another particle set (source).

``estimator type=latticedeviation`` element:

  +------------------+----------------------+
  | parent elements: | ``hamiltonian, qmc`` |
  +------------------+----------------------+
  | child elements:  | *None*               |
  +------------------+----------------------+

attributes:

  +-------------------------+--------------+------------------+-------------+------------------------------+
  | **Name**                | **Datatype** | **Values**       | **Default** | **Description**              |
  +=========================+==============+==================+=============+==============================+
  | ``type``:math:`^r`      | text         | latticedeviation |             | Must be latticedeviation     |
  +-------------------------+--------------+------------------+-------------+------------------------------+
  | ``name``:math:`^r`      | text         | *anything*       | any         | Unique name for estimator    |
  +-------------------------+--------------+------------------+-------------+------------------------------+
  | ``hdf5``:math:`^o`      | boolean      | yes/no           | no          | Output to ``stat.h5`` (yes)  |
  +-------------------------+--------------+------------------+-------------+------------------------------+
  | ``per_xyz``:math:`^o`   | boolean      | yes/no           | no          | Directionally resolved (yes) |
  +-------------------------+--------------+------------------+-------------+------------------------------+
  | ``source``:math:`^r`    | text         | e/ion0/...       | no          | source particleset           |
  +-------------------------+--------------+------------------+-------------+------------------------------+
  | ``sgroup``:math:`^r`    | text         | u/d/...          | no          | source particle group        |
  +-------------------------+--------------+------------------+-------------+------------------------------+
  | ``target``:math:`^r`    | text         | e/ion0/...       | no          | target particleset           |
  +-------------------------+--------------+------------------+-------------+------------------------------+
  | ``tgroup``:math:`^r`    | text         | u/d/...          | no          | target particle group        |
  +-------------------------+--------------+------------------+-------------+------------------------------+

Additional information:

-  ``source``: The “reference” particleset to measure distances from;
   actual reference points are determined together with ``sgroup``.

-  ``sgroup``: The “reference” particle group to measure distances from.

-  ``source``: The “target” particleset to measure distances to.

-  ``sgroup``: The “target” particle group to measure distances to. For
   example, in :ref:`Listing 32 <Listing 32>` the distance from the up
   electron (“u”) to the origin of the coordinate system is recorded.

-  ``per_xyz``: Used to record direction-resolved distance. In
   :ref:`Listing 32 <Listing 32>`, the x,y,z coordinates of the up electron
   will be recorded separately if ``per_xyz=yes``.

-  ``hdf5``: Used to record particle-resolved distances in the h5 file
   if ``gdf5=yes``.

.. code-block::
  :caption: Lattice deviation estimator element.
  :name: Listing 32

  <particleset name="e" random="yes">
    <group name="u" size="1" mass="1.0">
       <parameter name="charge"              >    -1                    </parameter>
       <parameter name="mass"                >    1.0                   </parameter>
    </group>
    <group name="d" size="1" mass="1.0">
       <parameter name="charge"              >    -1                    </parameter>
       <parameter name="mass"                >    1.0                   </parameter>
    </group>
  </particleset>

  <particleset name="wf_center">
    <group name="origin" size="1">
      <attrib name="position" datatype="posArray" condition="0">
               0.00000000        0.00000000        0.00000000
      </attrib>
    </group>
  </particleset>

  <estimator type="latticedeviation" name="latdev" hdf5="yes" per_xyz="yes"
    source="wf_center" sgroup="origin" target="e" tgroup="u"/>

Energy density estimator
~~~~~~~~~~~~~~~~~~~~~~~~

An energy density operator, :math:`\hat{\mathcal{E}}_r`, satisfies

.. math::
  :label: eq37

  \int dr \hat{\mathcal{E}}_r = \hat{H},

where the integral is over all space and :math:`\hat{H}` is the
Hamiltonian. In QMCPACK, the energy density is split into kinetic and potential
components

.. math::
  :label: eq38

  \hat{\mathcal{E}}_r = \hat{\mathcal{T}}_r + \hat{\mathcal{V}}_r\:,

with each component given by

.. math::
  :label: eq39

  \begin{aligned}
      \hat{\mathcal{T}}_r &=  \frac{1}{2}\sum_i\delta(r-r_i)\hat{p}_i^2 \\
      \hat{\mathcal{V}}_r &=  \sum_{i<j}\frac{\delta(r-r_i)+\delta(r-r_j)}{2}\hat{v}^{ee}(r_i,r_j)
                 + \sum_{i\ell}\frac{\delta(r-r_i)+\delta(r-\tilde{r}_\ell)}{2}\hat{v}^{eI}(r_i,\tilde{r}_\ell) \nonumber\\
       &\qquad   + \sum_{\ell< m}\frac{\delta(r-\tilde{r}_\ell)+\delta(r-\tilde{r}_m)}{2}\hat{v}^{II}(\tilde{r}_\ell,\tilde{r}_m)\:.\nonumber\end{aligned}

Here, :math:`r_i` and :math:`\tilde{r}_\ell` represent electron and ion
positions, respectively; :math:`\hat{p}_i` is a single electron momentum
operator; and :math:`\hat{v}^{ee}(r_i,r_j)`,
:math:`\hat{v}^{eI}(r_i,\tilde{r}_\ell)`, and
:math:`\hat{v}^{II}(\tilde{r}_\ell,\tilde{r}_m)` are the
electron-electron, electron-ion, and ion-ion pair potential operators
(including nonlocal pseudopotentials, if present). This form of the
energy density is size consistent; that is, the partially integrated
energy density operators of well-separated atoms gives the isolated
Hamiltonians of the respective atoms. For periodic systems with
twist-averaged boundary conditions, the energy density is formally
correct only for either a set of supercell k-points that correspond to
real-valued wavefunctions or a k-point set that has inversion symmetry
around a k-point having a real-valued wavefunction. For more information
about the energy density, see :cite:`Krogel2013`.

In QMCPACK, the energy density can be accumulated on piecewise uniform 3D grids in generalized Cartesian, cylindrical, or spherical coordinates.  The energy density integrated within Voronoi volumes centered on ion positions is also available.  The total particle number density is also accumulated on the same grids by the energy density estimator for convenience so that related quantities, such as the regional energy per particle, can be computed easily.

``estimator type=EnergyDensity`` element:

  +------------------+---------------------------------+
  | parent elements: | ``hamiltonian, qmc``            |
  +------------------+---------------------------------+
  | child elements:  | ``reference_points, spacegrid`` |
  +------------------+---------------------------------+

attributes:

  +----------------------------+--------------+----------------------+-------------+----------------------------------------------+
  | **Name**                   | **Datatype** | **Values**           | **Default** | **Description**                              |
  +============================+==============+======================+=============+==============================================+
  | ``type``:math:`^r`         | text         | **EnergyDensity**    |             | Must be EnergyDensity                        |
  +----------------------------+--------------+----------------------+-------------+----------------------------------------------+
  | ``name``:math:`^r`         | text         | *anything*           |             | Unique name for estimator                    |
  +----------------------------+--------------+----------------------+-------------+----------------------------------------------+
  | ``dynamic``:math:`^r`      | text         | ``particleset.name`` |             | Identify electrons                           |
  +----------------------------+--------------+----------------------+-------------+----------------------------------------------+
  | ``static``:math:`^o`       | text         | ``particleset.name`` |             | Identify ions                                |
  +----------------------------+--------------+----------------------+-------------+----------------------------------------------+
  | ``ion_points``:math:`^o`   | text         | yes/no               |  no         | Separate ion energy density onto point field |
  +----------------------------+--------------+----------------------+-------------+----------------------------------------------+

Additional information:

-  ``name:`` Must be unique. A dataset with blocked statistical data for
   the energy density will appear in the ``stat.h5`` files labeled as
   ``name``.
- **Important:** in order for the estimator to work, a traces XML input element (<traces array="yes" write="no"/>) must appear following the <qmcsystem/> element and prior to any <qmc/> element.

.. code-block::
  :caption: Energy density estimator accumulated on a :math:`20 \times  10 \times 10` grid over the simulation cell.
  :name: Listing 33

  <estimator type="EnergyDensity" name="EDcell" dynamic="e" static="ion0">
     <spacegrid coord="cartesian">
       <origin p1="zero"/>
       <axis p1="a1" scale=".5" label="x" grid="-1 (.05) 1"/>
       <axis p1="a2" scale=".5" label="y" grid="-1 (.1) 1"/>
       <axis p1="a3" scale=".5" label="z" grid="-1 (.1) 1"/>
     </spacegrid>
  </estimator>

.. code-block::
  :caption: Energy density estimator accumulated within spheres of radius 6.9 Bohr centered on the first and second atoms in the ion0 particleset.
  :name: Listing 34

  <estimator type="EnergyDensity" name="EDatom" dynamic="e" static="ion0">
    <reference_points coord="cartesian">
      r1 1 0 0
      r2 0 1 0
      r3 0 0 1
    </reference_points>
    <spacegrid coord="spherical">
      <origin p1="ion01"/>
      <axis p1="r1" scale="6.9" label="r"     grid="0 1"/>
      <axis p1="r2" scale="6.9" label="phi"   grid="0 1"/>
      <axis p1="r3" scale="6.9" label="theta" grid="0 1"/>
    </spacegrid>
    <spacegrid coord="spherical">
      <origin p1="ion02"/>
      <axis p1="r1" scale="6.9" label="r"     grid="0 1"/>
      <axis p1="r2" scale="6.9" label="phi"   grid="0 1"/>
      <axis p1="r3" scale="6.9" label="theta" grid="0 1"/>
    </spacegrid>
  </estimator>

.. code-block::
  :caption: Energy density estimator accumulated within Voronoi polyhedra centered on the ions.
  :name: Listing 35

  <estimator type="EnergyDensity" name="EDvoronoi" dynamic="e" static="ion0">
    <spacegrid coord="voronoi"/>
  </estimator>

The ``<reference_points/>`` element provides a set of points for later
use in specifying the origin and coordinate axes needed to construct a
spatial histogramming grid. Several reference points on the surface of
the simulation cell (see :numref:`table8`), as well as the
positions of the ions (see the ``energydensity.static`` attribute), are
made available by default. The reference points can be used, for
example, to construct a cylindrical grid along a bond with the origin on
the bond center.

  ``reference_points`` element:

    +------------------+---------------------------------+
    | parent elements: | ``estimator type=EnergyDensity``|
    +------------------+---------------------------------+
    | child elements:  | *None*                          |
    +------------------+---------------------------------+

  attributes:

    +----------------------+--------------+----------------+-------------+---------------------------+
    | **Name**             | **Datatype** | **Values**     | **Default** | **Description**           |
    +======================+==============+================+=============+===========================+
    | ``coord``:math:`^r`  | text         | Cartesian/cell |             | Specify coordinate system |
    +----------------------+--------------+----------------+-------------+---------------------------+

  body text: The body text is a line formatted list of points with labels

Additional information:

-  ``coord:`` If ``coord=cartesian``, labeled points are in Cartesian
   (x,y,z) format in units of Bohr. If ``coord=cell``, then labeled
   points are in units of the simulation cell axes.

-  ``body text:`` The list of points provided in the body text are line
   formatted, with four entries per line (*label* *coor1* *coor2*
   *coor3*). A set of points referenced to the simulation cell is
   available by default (see :numref:`table8`). If
   ``energydensity.static`` is provided, the location of each individual
   ion is also available (e.g., if ``energydensity.static=ion0``, then
   the location of the first atom is available with label ion01, the
   second with ion02, etc.). All points can be used by label when
   constructing spatial histogramming grids (see the following
   ``spacegrid`` element) used to collect energy densities.

.. _table8:
.. table::

     ========= ======================== =================
     ``label`` ``point``                ``description``
     ========= ======================== =================
     ``zero``  0 0 0                    Cell center
     ``a1``    :math:`a_1`              Cell axis 1
     ``a2``    :math:`a_2`              Cell axis 2
     ``a3``    :math:`a_3`              Cell axis 3
     ``f1p``   :math:`a_1`/2            Cell face 1+
     ``f1m``   -:math:`a_1`/2           Cell face 1-
     ``f2p``   :math:`a_2`/2            Cell face 2+
     ``f2m``   -:math:`a_2`/2           Cell face 2-
     ``f3p``   :math:`a_3`/2            Cell face 3+
     ``f3m``   -:math:`a_3`/2           Cell face 3-
     ``cppp``  :math:`(a_1+a_2+a_3)/2`  Cell corner +,+,+
     ``cppm``  :math:`(a_1+a_2-a_3)/2`  Cell corner +,+,-
     ``cpmp``  :math:`(a_1-a_2+a_3)/2`  Cell corner +,-,+
     ``cmpp``  :math:`(-a_1+a_2+a_3)/2` Cell corner -,+,+
     ``cpmm``  :math:`(a_1-a_2-a_3)/2`  Cell corner +,-,-
     ``cmpm``  :math:`(-a_1+a_2-a_3)/2` Cell corner -,+,-
     ``cmmp``  :math:`(-a_1-a_2+a_3)/2` Cell corner -,-,+
     ``cmmm``  :math:`(-a_1-a_2-a_3)/2` Cell corner -,-,-
     ========= ======================== =================

.. centered:: Table 8 Reference points available by default. Vectors :math:`a_1`, :math:`a_2`, and :math:`a_3` refer to the simulation cell axes. The representation of the cell is centered around ``zero``.

The ``<spacegrid/>`` element is used to specify a spatial histogramming
grid for the energy density. Grids are constructed based on a set of,
potentially nonorthogonal, user-provided coordinate axes. The axes are
based on information available from ``reference_points``. Voronoi grids
are based only on nearest neighbor distances between electrons and ions.
Any number of space grids can be provided to a single energy density
estimator.

``spacegrid`` element:

  +------------------+---------------------------------+
  | parent elements: | ``estimator type=EnergyDensity``|
  +------------------+---------------------------------+
  | child elements:  | ``origin, axis``                |
  +------------------+---------------------------------+

attributes:

  +----------------------+--------------+--------------+-------------+---------------------------+
  | **Name**             | **Datatype** | **Values**   | **Default** | **Description**           |
  +======================+==============+==============+=============+===========================+
  | ``coord``:math:`^r`  | text         | Cartesian    |             | Specify coordinate system |
  +----------------------+--------------+--------------+-------------+---------------------------+
  |                      |              | cylindrical  |             |                           |
  +----------------------+--------------+--------------+-------------+---------------------------+
  |                      |              | spherical    |             |                           |
  +----------------------+--------------+--------------+-------------+---------------------------+
  |                      |              | Voronoi      |             |                           |
  +----------------------+--------------+--------------+-------------+---------------------------+

The ``<origin/>`` element gives the location of the origin for a
non-Voronoi grid.

Additional information:

-  ``p1/p2/fraction:`` The location of the origin is set to
   ``p1+fraction*(p2-p1)``. If only ``p1`` is provided, the origin is at
   ``p1``.

``origin`` element:

  +------------------+---------------------------------+
  | parent elements: | ``spacegrid``                   |
  +------------------+---------------------------------+
  | child elements:  | *None*                          |
  +------------------+---------------------------------+

attributes:

  +-------------------------+--------------+----------------------------+-------------+------------------------+
  | **Name**                | **Datatype** | **Values**                 | **Default** | **Description**        |
  +=========================+==============+============================+=============+========================+
  | ``p1``:math:`^r`        | text         | ``reference_point.label``  |             | Select end point       |
  +-------------------------+--------------+----------------------------+-------------+------------------------+
  | ``p2``:math:`^o`        | text         | ``reference_point.label``  |             | Select end point       |
  +-------------------------+--------------+----------------------------+-------------+------------------------+
  | ``fraction``:math:`^o`  | real         |                            | 0           | Interpolation fraction |
  +-------------------------+--------------+----------------------------+-------------+------------------------+

The ``<axis/>`` element represents a coordinate axis used to construct the, possibly curved, coordinate system for the histogramming grid.  Three ``<axis/>`` elements must be provided to a non-Voronoi ``<spacegrid/>`` element.

``axis`` element:

  +------------------+---------------------------------+
  | parent elements: | ``spacegrid``                   |
  +------------------+---------------------------------+
  | child elements:  | *None*                          |
  +------------------+---------------------------------+

attributes:

  +----------------------+--------------+----------------------------+-------------+------------------------+
  | **Name**             | **Datatype** | **Values**                 | **Default** | **Description**        |
  +======================+==============+============================+=============+========================+
  | ``label``:math:`^r`  | text         | *See below*                |             | Axis/dimension label   |
  +----------------------+--------------+----------------------------+-------------+------------------------+
  | ``grid``:math:`^r`   | text         |                            | "0 1"       | Grid ranges/intervals  |
  +----------------------+--------------+----------------------------+-------------+------------------------+
  | ``p1``:math:`^r`     | text         | ``reference_point.label``  |             | Select end point       |
  +----------------------+--------------+----------------------------+-------------+------------------------+
  | ``p2``:math:`^o`     | text         | ``reference_point.label``  |             | Select end point       |
  +----------------------+--------------+----------------------------+-------------+------------------------+
  | ``scale``:math:`^o`  | real         |                            |             | Interpolation fraction |
  +----------------------+--------------+----------------------------+-------------+------------------------+

Additional information:

-  ``label:`` The allowed set of axis labels depends on the coordinate
   system (i.e., ``spacegrid.coord``). Labels are ``x/y/z`` for
   ``coord=cartesian``, ``r/phi/z`` for ``coord=cylindrical``,
   ``r/phi/theta`` for ``coord=spherical``.

-  ``p1/p2/scale:`` The axis vector is set to ``p1+scale*(p2-p1)``. If
   only ``p1`` is provided, the axis vector is ``p1``.

-  ``grid:`` The grid specifies the histogram grid along the direction
   specified by ``label``. The allowed grid points fall in the range
   [-1,1] for ``label=x/y/z`` or [0,1] for ``r/phi/theta``. A grid of 10
   evenly spaced points between 0 and 1 can be requested equivalently by
   ``grid="0 (0.1) 1"`` or ``grid="0 (10) 1."`` Piecewise uniform grids
   covering portions of the range are supported, e.g.,
   ``grid="-0.7 (10) 0.0 (20) 0.5."``

-  Note that ``grid`` specifies the histogram grid along the (curved)
   coordinate given by ``label``. The axis specified by ``p1/p2/scale``
   does not correspond one-to-one with ``label`` unless ``label=x/y/z``,
   but the full set of axes provided defines the (sheared) space on top
   of which the curved (e.g., spherical) coordinate system is built.

One body density matrix
~~~~~~~~~~~~~~~~~~~~~~~

The N-body density matrix in DMC is
:math:`\hat{\rho}_N=\left|{\Psi_{T}}\rangle{}\langle{\Psi_{FN}}\right|` (for VMC,
substitute :math:`\Psi_T` for :math:`\Psi_{FN}`). The one body reduced
density matrix (1RDM) is obtained by tracing out all particle
coordinates but one:

.. math::
  :label: eq40

  \hat{n}_1 = \sum_nTr_{R_n}\left|{\Psi_{T}}\rangle{}\langle{\Psi_{FN}}\right|

In this formula, the sum is over all electron indices and
:math:`Tr_{R_n}(*)\equiv\int dR_n\langle{R_n}\left|{*}\right|{R_n}\rangle` with
:math:`R_n=[r_1,...,r_{n-1},r_{n+1},...,r_N]`. When the sum is
restricted over spin-up or spin-down electrons, one obtains a density
matrix for each spin species. The 1RDM computed by is partitioned in
this way.

In real space, the matrix elements of the 1RDM are

.. math::
  :label: eq41

  \begin{aligned}
     n_1(r,r') &= \langle{r}\left|{\hat{n}_1}\right|{r'}\rangle = \sum_n\int dR_n \Psi_T(r,R_n)\Psi_{FN}^*(r',R_n)\:. \end{aligned}

A more efficient and compact representation of the 1RDM is obtained by
expanding in the SPOs obtained from a Hartree-Fock or DFT calculation,
:math:`\{\phi_i\}`:

.. math::
  :label: eq42

  n_1(i,j) &= \langle{\phi_i}\left|{\hat{n}_1}\right|{\phi_j}\rangle \nonumber \\
           &= \int dR \Psi_{FN}^*(R)\Psi_{T}(R) \sum_n\int dr'_n \frac{\Psi_T(r_n',R_n)}{\Psi_T(r_n,R_n)}\phi_i(r_n')^* \phi_j(r_n)\:.

The integration over :math:`r'` in :eq:`eq42` is inefficient when one is also interested in obtaining matrices involving energetic quantities, such as the energy density matrix of :cite:`Krogel2014` or the related (and more well known) generalized Fock matrix.  For this reason, an approximation is introduced as follows:

.. math::
  :label: eq43

  \begin{aligned}
      n_1(i,j) \approx \int dR \Psi_{FN}(R)^*\Psi_T(R)  \sum_n \int dr_n' \frac{\Psi_T(r_n',R_n)^*}{\Psi_T(r_n,R_n)^*}\phi_i(r_n)^* \phi_j(r_n')\:. \end{aligned}

For VMC, FN-DMC, FP-DMC, and RN-DMC this formula represents an exact
sampling of the 1RDM corresponding to :math:`\hat{\rho}_N^\dagger` (see
appendix A of :cite:`Krogel2014` for more detail).

``estimtor type=dm1b`` element:

  +------------------+----------------------+
  | parent elements: | ``hamiltonian, qmc`` |
  +------------------+----------------------+
  | child elements:  | *None*               |
  +------------------+----------------------+

attributes:

  +---------------------+--------------+------------+-------------+---------------------------+
  | **Name**            | **Datatype** | **Values** | **Default** | **Description**           |
  +=====================+==============+============+=============+===========================+
  | ``type``:math:`^r`  | text         | **dm1b**   |             | Must be dm1b              |
  +---------------------+--------------+------------+-------------+---------------------------+
  | ``name``:math:`^r`  | text         | *anything* |             | Unique name for estimator |
  +---------------------+--------------+------------+-------------+---------------------------+

parameters:

  +-----------------------------------+---------------+-------------------------------+---------------+---------------------------+
  | **Name**                          | **Datatype**  | **Values**                    | **Default**   | **Description**           |
  +===================================+===============+===============================+===============+===========================+
  | ``basis``:math:`^r`               | text array    | sposet.name(s)                |               | Orbital basis             |
  +-----------------------------------+---------------+-------------------------------+---------------+---------------------------+
  | ``integrator``:math:`^o`          | text          | uniform_grid uniform density  | uniform_grid  | Integration method        |
  +-----------------------------------+---------------+-------------------------------+---------------+---------------------------+
  | ``evaluator``:math:`^o`           | text          | loop/matrix                   | loop          | Evaluation method         |
  +-----------------------------------+---------------+-------------------------------+---------------+---------------------------+
  | ``scale``:math:`^o`               | real          | :math:`0<scale<1`             | 1.0           | Scale integration cell    |
  +-----------------------------------+---------------+-------------------------------+---------------+---------------------------+
  | ``center``:math:`^o`              | real array(3) | *any point*                   |               | Center of cell            |
  +-----------------------------------+---------------+-------------------------------+---------------+---------------------------+
  | ``points``:math:`^o`              | integer       | :math:`>0`                    | 10            | Grid points in each dim   |
  +-----------------------------------+---------------+-------------------------------+---------------+---------------------------+
  | ``samples``:math:`^o`             | integer       | :math:`>0`                    | 10            | MC samples                |
  +-----------------------------------+---------------+-------------------------------+---------------+---------------------------+
  | ``warmup``:math:`^o`              | integer       | :math:`>0`                    | 30            | MC warmup                 |
  +-----------------------------------+---------------+-------------------------------+---------------+---------------------------+
  | ``timestep``:math:`^o`            | real          | :math:`>0`                    | 0.5           | MC time step              |
  +-----------------------------------+---------------+-------------------------------+---------------+---------------------------+
  | ``use_drift``:math:`^o`           | boolean       | yes/no                        | no            | Use drift in VMC          |
  +-----------------------------------+---------------+-------------------------------+---------------+---------------------------+
  | ``check_overlap``:math:`^o`       | boolean       | yes/no                        | no            | Print overlap matrix      |
  +-----------------------------------+---------------+-------------------------------+---------------+---------------------------+
  | ``check_derivatives``:math:`^o`   | boolean       | yes/no                        | no            | Check density derivatives |
  +-----------------------------------+---------------+-------------------------------+---------------+---------------------------+
  | ``acceptance_ratio``:math:`^o`    | boolean       | yes/no                        | no            | Print accept ratio        |
  +-----------------------------------+---------------+-------------------------------+---------------+---------------------------+
  | ``rstats``:math:`^o`              | boolean       | yes/no                        | no            | Print spatial stats       |
  +-----------------------------------+---------------+-------------------------------+---------------+---------------------------+
  | ``normalized``:math:`^o`          | boolean       | yes/no                        | yes           | ``basis`` comes norm'ed   |
  +-----------------------------------+---------------+-------------------------------+---------------+---------------------------+
  | ``volume_normed``:math:`^o`       | boolean       | yes/no                        | yes           | ``basis`` norm is volume  |
  +-----------------------------------+---------------+-------------------------------+---------------+---------------------------+
  | ``energy_matrix``:math:`^o`       | boolean       | yes/no                        | no            | Energy density matrix     |
  +-----------------------------------+---------------+-------------------------------+---------------+---------------------------+

Additional information:

-  ``name:`` Density matrix results appear in ``stat.h5`` files labeled
   according to ``name``.

-  ``basis:`` List ``sposet.name``\ ’s. The total set of orbitals
   contained in all ``sposet``\ ’s comprises the basis (subspace) onto
   which the one body density matrix is projected. This set of orbitals
   generally includes many virtual orbitals that are not occupied in a
   single reference Slater determinant.

-  ``integrator:`` Select the method used to perform the additional
   single particle integration. Options are ``uniform_grid`` (uniform
   grid of points over the cell), ``uniform`` (uniform random sampling
   over the cell), and ``density`` (Metropolis sampling of approximate
   density, :math:`\sum_{b\in \texttt{basis}}\left|{\phi_b}\right|^2`, is not
   well tested, please check results carefully!). Depending on the
   integrator selected, different subsets of the other input parameters
   are active.

-  ``evaluator:`` Select for-loop or matrix multiply implementations.
   Matrix is preferred for speed. Both implementations should give the
   same results, but please check as this has not been exhaustively
   tested.

-  ``scale:`` Resize the simulation cell by scale for use as an
   integration volume (active for ``integrator=uniform/uniform_grid``).

-  ``center:`` Translate the integration volume to center at this point
   (active for ``integrator=uniform/ uniform_grid``). If ``center`` is
   not provided, the scaled simulation cell is used as is.

-  ``points:`` Number of grid points in each dimension for
   ``integrator=uniform_grid``. For example, ``points=10`` results in a
   uniform :math:`10 \times 10 \times 10` grid over the cell.

-  ``samples:`` Sets the number of MC samples collected for each step
   (active for ``integrator=uniform/ density``).

-  ``warmup:`` Number of warmup Metropolis steps at the start of the run
   before data collection (active for ``integrator=density``).

-  ``timestep:`` Drift-diffusion time step used in Metropolis sampling
   (active for ``integrator=density``).

-  ``use_drift:`` Enable drift in Metropolis sampling (active for
   ``integrator=density``).

-  ``check_overlap:`` Print the overlap matrix (computed via simple
   Riemann sums) to the log, then abort. Note that subsequent analysis
   based on the 1RDM is simplest if the input orbitals are orthogonal.

-  ``check_derivatives:`` Print analytic and numerical derivatives of
   the approximate (sampled) density for several sample points, then
   abort.

-  ``acceptance_ratio:`` Print the acceptance ratio of the density
   sampling to the log for each step.

-  ``rstats:`` Print statistical information about the spatial motion of
   the sampled points to the log for each step.

-  ``normalized:`` Declare whether the inputted orbitals are normalized
   or not. If ``normalized=no``, direct Riemann integration over a
   :math:`200 \times 200 \times 200` grid will be used to compute the
   normalizations before use.

-  ``volume_normed:`` Declare whether the inputted orbitals are
   normalized to the cell volume (default) or not (a norm of 1.0 is
   assumed in this case). Currently, B-spline orbitals coming from QE
   and HEG planewave orbitals native to QMCPACK are known to be volume
   normalized.

-  ``energy_matrix:`` Accumulate the one body reduced energy density
   matrix, and write it to ``stat.h5``. This matrix is not covered in
   any detail here; the interested reader is referred to
   :cite:`Krogel2014`.

.. code-block::
  :caption: One body density matrix with uniform grid integration.
  :name: Listing 36

  <estimator type="dm1b" name="DensityMatrices">
    <parameter name="basis"        >  spo_u spo_uv  </parameter>
    <parameter name="evaluator"    >  matrix        </parameter>
    <parameter name="integrator"   >  uniform_grid  </parameter>
    <parameter name="points"       >  4             </parameter>
    <parameter name="scale"        >  1.0           </parameter>
    <parameter name="center"       >  0 0 0         </parameter>
  </estimator>

.. code-block::
  :caption: One body density matrix with uniform sampling.
  :name: Listing 37

  <estimator type="dm1b" name="DensityMatrices">
    <parameter name="basis"        >  spo_u spo_uv  </parameter>
    <parameter name="evaluator"    >  matrix        </parameter>
    <parameter name="integrator"   >  uniform       </parameter>
    <parameter name="samples"      >  64            </parameter>
    <parameter name="scale"        >  1.0           </parameter>
    <parameter name="center"       >  0 0 0         </parameter>
  </estimator>

.. code-block::
  :caption: One body density matrix with density sampling.
  :name: Listing 38

  <estimator type="dm1b" name="DensityMatrices">
    <parameter name="basis"        >  spo_u spo_uv  </parameter>
    <parameter name="evaluator"    >  matrix        </parameter>
    <parameter name="integrator"   >  density       </parameter>
    <parameter name="samples"      >  64            </parameter>
    <parameter name="timestep"     >  0.5           </parameter>
    <parameter name="use_drift"    >  no            </parameter>
  </estimator>

.. code-block::
  :caption: Example ``sposet`` initialization for density matrix use.  Occupied and virtual orbital sets are created separately, then joined (``basis="spo_u spo_uv"``).
  :name: Listing 39

  <sposet_builder type="bspline" href="../dft/pwscf_output/pwscf.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" meshfactor="1.0" gpu="no" precision="single">
    <sposet type="bspline" name="spo_u"  group="0" size="4"/>
    <sposet type="bspline" name="spo_d"  group="0" size="2"/>
    <sposet type="bspline" name="spo_uv" group="0" index_min="4" index_max="10"/>
  </sposet_builder>

.. code-block::
  :caption: Example ``sposet`` initialization for density matrix use. Density matrix orbital basis created separately (``basis="dm_basis"``).
  :name: Listing 40

  <sposet_builder type="bspline" href="../dft/pwscf_output/pwscf.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" meshfactor="1.0" gpu="no" precision="single">
    <sposet type="bspline" name="spo_u"  group="0" size="4"/>
    <sposet type="bspline" name="spo_d"  group="0" size="2"/>
    <sposet type="bspline" name="dm_basis" size="50" spindataset="0"/>
  </sposet_builder>

.. _forward-walking:

Forward-Walking Estimators
--------------------------

Forward walking is a method for sampling the pure fixed-node
distribution :math:`\langle \Phi_0 | \Phi_0\rangle`. Specifically, one
multiplies each walker’s DMC mixed estimate for the observable
:math:`\mathcal{O}`,
:math:`\frac{\mathcal{O}(\mathbf{R})\Psi_T(\mathbf{R})}{\Psi_T(\mathbf{R})}`,
by the weighting factor
:math:`\frac{\Phi_0(\mathbf{R})}{\Psi_T(\mathbf{R})}`. As it turns out,
this weighting factor for any walker :math:`\mathbf{R}` is proportional
to the total number of descendants the walker will have after a
sufficiently long projection time :math:`\beta`.

To forward walk on an observable, declare a generic forward-walking
estimator within a ``<hamiltonian>`` block, and then specify the
observables to forward walk on and the forward-walking parameters. Here
is a summary.

``estimator type=ForwardWalking`` element:

  +------------------+----------------------+
  | parent elements: | ``hamiltonian, qmc`` |
  +------------------+----------------------+
  | child elements:  | ``Observable``       |
  +------------------+----------------------+

  attributes:

    +---------------------+--------------+--------------------+-------------+---------------------------+
    | **Name**            | **Datatype** | **Values**         | **Default** | **Description**           |
    +=====================+==============+====================+=============+===========================+
    | ``type``:math:`^r`  | text         | **ForwardWalking** |             | Must be "ForwardWalking"  |
    +---------------------+--------------+--------------------+-------------+---------------------------+
    | ``name``:math:`^r`  | text         | *anything*         | any         | Unique name for estimator |
    +---------------------+--------------+--------------------+-------------+---------------------------+

``Observable`` element:

  +------------------+---------------------------------+
  | parent elements: | ``estimator, hamiltonian, qmc`` |
  +------------------+---------------------------------+
  | child elements:  | *None*                          |
  +------------------+---------------------------------+

    +--------------------------+--------------+---------------+-------------+---------------------------------------------------------------------------------+
    | **Name**                 | **Datatype** | **Values**    | **Default** | **Description**                                                                 |
    +==========================+==============+===============+=============+=================================================================================+
    | ``name``:math:`^r`       | text         | *anything*    | any         | Registered name of existing estimator on which to forward walk                  |
    +--------------------------+--------------+---------------+-------------+---------------------------------------------------------------------------------+
    | ``max``:math:`^r`        | integer      | :math:`>0`    |             | Maximum projection time in steps (``max``:math:`=\beta/\tau`)                   |
    +--------------------------+--------------+---------------+-------------+---------------------------------------------------------------------------------+
    | ``frequency``:math:`^r`  | text         | :math:`\geq 1`|             | Dump data only for every ``frequency``-th to ``scalar.dat`` file                |
    +--------------------------+--------------+---------------+-------------+---------------------------------------------------------------------------------+

Additional information:

-  **Cost**: Because histories of observables up to ``max`` time steps
   have to be stored, the memory cost of storing the nonforward-walked
   observables variables should be multiplied by :math:`\texttt{max}`.
   Although this is not an issue for items such as potential energy, it
   could be prohibitive for observables such as density, forces, etc.

-  **Naming Convention**: Forward-walked observables are automatically
   named ``FWE_name_i``, where ``i`` is the forward-walked expectation
   value at time step ``i``, and ``name`` is whatever name appears in
   the ``<Observable>`` block. This is also how it will appear in the
   ``scalar.dat`` file.

In the following example case, QMCPACK forward walks on the potential
energy for 300 time steps and dumps the forward-walked value at every
time step.

.. code-block::
  :caption: Forward-walking estimator element.
  :name: Listing 41

  <estimator name="fw" type="ForwardWalking">
      <Observable name="LocalPotential" max="300" frequency="1"/>
       <!--- Additional Observable blocks go here -->
   </estimator>

.. _force-est:

"Force" estimators
------------------

QMCPACK supports force estimation by use of the Chiesa-Ceperly-Zhang
(CCZ) estimator. Currently, open and periodic boundary conditions are
supported but for all-electron calculations only.

Without loss of generality, the CCZ estimator for the z-component of the
force on an ion centered at the origin is given by the following
expression:

.. math::
  :label: eq44

  F_z = -Z \sum_{i=1}^{N_e}\frac{z_i}{r_i^3}[\theta(r_i-\mathcal{R}) + \theta(\mathcal{R}-r_i)\sum_{\ell=1}^{M}c_\ell r_i^\ell]\:.

Z is the ionic charge, :math:`M` is the degree of the smoothing
polynomial, :math:`\mathcal{R}` is a real-space cutoff of the sphere
within which the bare-force estimator is smoothed, and :math:`c_\ell`
are predetermined coefficients. These coefficients are chosen to
minimize the weighted mean square error between the bare force estimate
and the s-wave filtered estimator. Specifically,

.. math::
  :label: eq45

  \chi^2 = \int_0^\mathcal{R}dr\,r^m\,[f_z(r) - \tilde{f}_z(r)]^2\:.

Here, :math:`m` is the weighting exponent, :math:`f_z(r)` is the
unfiltered radial force density for the z force component, and
:math:`\tilde{f}_z(r)` is the smoothed polynomial function for the same
force density. The reader is invited to refer to the original paper for
a more thorough explanation of the methodology, but with the notation in
hand, QMCPACK takes the following parameters.

``estimator type=Force`` element:

  +------------------+----------------------+
  | parent elements: | ``hamiltonian, qmc`` |
  +------------------+----------------------+
  | child elements:  | ``parameter``        |
  +------------------+----------------------+

  attributes:

    +--------------------------+--------------+-----------------+-------------+--------------------------------------------------------------+
    | **Name**                 | **Datatype** | **Values**      | **Default** | **Description**                                              |
    +==========================+==============+=================+=============+==============================================================+
    | ``mode``:math:`^o`       | text         | *See above*     | bare        | Select estimator type                                        |
    +--------------------------+--------------+-----------------+-------------+--------------------------------------------------------------+
    | ``lrmethod``:math:`^o`   | text         | ewald or srcoul | ewald       | Select long-range potential breakup method                   |
    +--------------------------+--------------+-----------------+-------------+--------------------------------------------------------------+
    | ``type``:math:`^r`       | text         | Force           |             | Must be "Force"                                              |
    +--------------------------+--------------+-----------------+-------------+--------------------------------------------------------------+
    | ``name``:math:`^o`       | text         | *Anything*      | ForceBase   | Unique name for this estimator                               |
    +--------------------------+--------------+-----------------+-------------+--------------------------------------------------------------+
    | ``pbc``:math:`^o`        | boolean      | yes/no          | yes         | Using periodic BCs or not                                    |
    +--------------------------+--------------+-----------------+-------------+--------------------------------------------------------------+
    | ``addionion``:math:`^o`  | boolean      | yes/no          | no          | Add the ion-ion force contribution to output force estimate  |
    +--------------------------+--------------+-----------------+-------------+--------------------------------------------------------------+

  parameters:

    +--------------------------+--------------+------------+-------------+----------------------------------------------------------+
    | **Name**                 | **Datatype** | **Values** | **Default** | **Description**                                          |
    +==========================+==============+============+=============+==========================================================+
    | ``rcut``:math:`^o`       | real         | :math:`>0` | 1.0         | Real-space cutoff :math:`\mathcal{R}` in bohr            |
    +--------------------------+--------------+------------+-------------+----------------------------------------------------------+
    | ``nbasis``:math:`^o`     | integer      | :math:`>0` | 2           | Degree of smoothing polynomial :math:`M`                 |
    +--------------------------+--------------+------------+-------------+----------------------------------------------------------+
    | ``weightexp``:math:`^o`  | integer      | :math:`>0` | 2           | :math:`\chi^2` weighting exponent :math`m`               |
    +--------------------------+--------------+------------+-------------+----------------------------------------------------------+

Additional information:

-  **Naming Convention**: The unique identifier ``name`` is appended
   with ``name_X_Y`` in the ``scalar.dat`` file, where ``X`` is the ion
   ID number and ``Y`` is the component ID (an integer with x=0, y=1,
   z=2). All force components for all ions are computed and dumped to
   the ``scalar.dat`` file.

-  **Long-range breakup**: With periodic boundary conditions, it is
   important to converge the lattice sum when calculating Coulomb
   contribution to the forces. As a quick test, increase the
   ``LR_dim_cutoff`` parameter until ion-ion forces are converged. The
   Ewald method converges more slowly than optimized method, but the
   optimized method can break down in edge cases, eg. too large
   ``LR_dim_cutoff``.

-  **Miscellaneous**: Usually, the default choice of ``weightexp`` is
   sufficient. Different combinations of ``rcut`` and ``nbasis`` should
   be tested though to minimize variance and bias. There is, of course,
   a tradeoff, with larger ``nbasis`` and smaller ``rcut`` leading to
   smaller biases and larger variances.

The following is an example use case.

::

  <simulationcell>
    ...
    <parameter name="LR_handler">  opt_breakup_original  </parameter>
    <parameter name="LR_dim_cutoff">  20  </parameter>
  </simulationcell>
  <hamiltonian>
    <estimator name="F" type="Force" mode="cep" addionion="yes">
      <parameter name="rcut">0.1</parameter>
      <parameter name="nbasis">4</parameter>
      <parameter name="weightexp">2</parameter>
    </estimator>
  </hamiltonian>

.. _stress-est:

Stress estimators
------------------

QMCPACK takes the following parameters.

  +------------------+----------------------+
  | parent elements: | ``hamiltonian``      |
  +------------------+----------------------+

  attributes:

    +--------------------------+--------------+-----------------+-------------+--------------------------------------------------------------+
    | **Name**                 | **Datatype** | **Values**      | **Default** | **Description**                                              |
    +==========================+==============+=================+=============+==============================================================+
    | ``mode``:math:`^r`       | text         | stress          | bare        | Must be "stress"                                             |
    +--------------------------+--------------+-----------------+-------------+--------------------------------------------------------------+
    | ``type``:math:`^r`       | text         | Force           |             | Must be "Force"                                              |
    +--------------------------+--------------+-----------------+-------------+--------------------------------------------------------------+
    | ``source``:math:`^r`     | text         | ion0            |             | Name of ion particleset                                      |
    +--------------------------+--------------+-----------------+-------------+--------------------------------------------------------------+
    | ``name``:math:`^o`       | text         | *Anything*      | ForceBase   | Unique name for this estimator                               |
    +--------------------------+--------------+-----------------+-------------+--------------------------------------------------------------+
    | ``addionion``:math:`^o`  | boolean      | yes/no          | no          | Add the ion-ion stress contribution to output                |
    +--------------------------+--------------+-----------------+-------------+--------------------------------------------------------------+

Additional information:

-  **Naming Convention**: The unique identifier ``name`` is appended
   with ``name_X_Y`` in the ``scalar.dat`` file, where ``X`` and ``Y``
   are the component IDs (an integer with x=0, y=1, z=2).

-  **Long-range breakup**: With periodic boundary conditions, it is
   important to converge the lattice sum when calculating Coulomb
   contribution to the forces. As a quick test, increase the
   ``LR_dim_cutoff`` parameter until ion-ion stresses are converged.
   Check using QE "Ewald contribution", for example. The stress
   estimator is implemented only with the Ewald method.

The following is an example use case.

::

  <simulationcell>
    ...
    <parameter name="LR_handler">  ewald  </parameter>
    <parameter name="LR_dim_cutoff">  45  </parameter>
  </simulationcell>
  <hamiltonian>
    <estimator name="S" type="Force" mode="stress" source="ion0"/>
  </hamiltonian>

.. bibliography:: /bibs/hamiltonianobservable.bib
