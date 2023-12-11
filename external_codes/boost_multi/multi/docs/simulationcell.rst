.. _simulationcell:

Specifying the system to be simulated
=====================================

Specifying the Simulation Cell
------------------------------

The :code:`simulationcell` block specifies the geometry of the cell, how the boundary
conditions should be handled, and how ewald summation should be broken
up.

``simulationcell`` Element:

  +------------------+-------------------------------------------------------------------------------------------------------+
  | Parent elements: | :code:`qmcsystem`                                                                                     |
  +------------------+-------------------------------------------------------------------------------------------------------+
  | Child elements:  | None                                                                                                  |
  +------------------+-------------------------------------------------------------------------------------------------------+

Attribute:

  +---------------------+--------------+---------------------------+-------------------+----------------------------------------------------+
  | **parameter name**  | **datatype** | **values**                | **default**       | **description**                                    |
  +=====================+==============+===========================+===================+====================================================+
  | ``lattice``         | 9 floats     | any float                 | Must be specified | Specification of lattice vectors.                  |
  +---------------------+--------------+---------------------------+-------------------+----------------------------------------------------+
  | ``bconds``          | string       | "p" or "n"                | "n n n "          | Boundary conditions for each axis.                 |
  +---------------------+--------------+---------------------------+-------------------+----------------------------------------------------+
  | ``vacuum``          | float        | :math:`\geq 1.0`          | 1.0               | Vacuum scale.                                      |
  +---------------------+--------------+---------------------------+-------------------+----------------------------------------------------+
  | ``LR_handler``      | string       | string                    | "opt_breakup"     | Ewald breakup method.                              |
  +---------------------+--------------+---------------------------+-------------------+----------------------------------------------------+
  | ``LR_dim_cutoff``   | float        | float                     | 15                | Ewald breakup distance.                            |
  +---------------------+--------------+---------------------------+-------------------+----------------------------------------------------+
  | ``LR_tol``          | float        | float                     | 3e-4              | Tolerance in Ha for Ewald ion-ion energy per atom. |
  +---------------------+--------------+---------------------------+-------------------+----------------------------------------------------+


An example of a block is given below:

::

   <simulationcell>
       <parameter name="lattice">
         3.8       0.0       0.0
         0.0       3.8       0.0
         0.0       0.0       3.8
       </parameter>
       <parameter name="bconds">
          p p p
       </parameter>
       <parameter name="LR_dim_cutoff"> 20 </parameter>
     </simulationcell>

Here, a cubic cell 3.8 bohr on a side will be used. This simulation will
use periodic boundary conditions, and the maximum :math:`k` vector will
be :math:`20/r_{wigner-seitz}` of the cell.

Lattice
~~~~~~~

The cell is specified using 3 lattice vectors.

Boundary conditions
~~~~~~~~~~~~~~~~~~~

QMCPACK offers the capability to use a mixture of open and periodic
boundary conditions. The parameter expects a single string of three
characters separated by spaces, *e.g.* “p p p” for purely periodic
boundary conditions. These characters control the behavior of the
:math:`x`, :math:`y`, and :math:`z`, axes, respectively. Non periodic
directions must be placed after the periodic ones. The only supported
combinations are:

**“p p p”** Periodic boundary conditions. Corresponds to a 3D crystal.

**“p p n”** Slab geometry. Corresponds to a 2D crystal.

**“p n n”** Wire geometry. Corresponds to a 1D crystal.

**“n n n”**
Open boundary conditions. Corresponds to an isolated molecule in a vacuum.

Vacuum
~~~~~~

The vacuum option allows adding a vacuum region in slab or wire boundary
conditions (``bconds= p p n`` or ``bconds= p n n``, respectively). The main use is to save memory with
spline or plane-wave basis trial wavefunctions, because no basis
functions are required inside the vacuum region. For example, a large
vacuum region can be added above and below a graphene sheet without
having to generate the trial wavefunction in such a large box or to have
as many splines as would otherwise be required. Note that the trial
wavefunction must still be generated in a large enough box to
sufficiently reduce periodic interactions in the underlying electronic
structure calculation.

With the vacuum option, the box used for Ewald summation increases along
the axis labeled by a factor of ``vacuum``. Note that all the particles remain in
the original box without altering their positions. i.e. Bond lengths are
not changed by this option. The default value is 1, no change to the
specified axes.

An example of a ``simulationcell`` block using is given below. The size of the box along
the z-axis increases from 12 to 18 by the vacuum scale of 1.5.

::

   <simulationcell>
       <parameter name="lattice">
         3.8       0.0       0.0
         0.0       3.8       0.0
         0.0       0.0      12.0
       </parameter>
       <parameter name="bconds">
          p p n
       </parameter>
       <parameter name="vacuum"> 1.5 </parameter>
       <parameter name="LR_dim_cutoff"> 20 </parameter>
       <parameter name="LR_handler"> ewald </parameter>
     </simulationcell>

LR_handler
~~~~~~~~~~

When using periodic boundary conditions direct calculation of the
Coulomb energy is conditionally convergent. As a result, QMCPACK uses an
optimized short-range/long-range breakup technique to compute the Coulomb
interaction in a rapidly convergent lattice sum. :cite:`Natoli1995`

In this summation, the energy is broken into short- and long-ranged
terms. The short-ranged term is computed directly in real space, while
the long-ranged term is computed in reciprocal space.

.. math:: v(r) = 1/r = v^{sr}(r) + v^{lr}(r)

`LR_handler` determines the functional form of :math:`v^{sr}` and :math:`v^{lr}`.
For example, the Ewald forms are

.. math:: v^{sr}(r) = \text{erfc}(\alpha r)/r

.. math:: v^{lr}(r) = \text{erf}(\alpha r)/r

Implemented choices for 3D systems are: ``ewald``, ``opt_breakup``, and ``opt_breakup_original``.
The choice for a 2D system is ``ewald_strict2d``.
The choice for a quasi-2D (e.g. slab) system is ``ewald_quasi2d``.

LR_dim_cutoff
~~~~~~~~~~~~~

QMCPACK chooses the short-range part to terminate at the image radius of
the simulation cell. This way only one real-space cell needs to be considered
using the minimum image convention.
`LR_dim_cutoff` controls the number of terms to include in the long-range sum.
The real-space cutoff :math:`r_{c}` and reciprocal-space cutoff :math:`k_{c}` are related by

.. math:: \mathrm{LR\_dim\_cutoff} = r_{c} \times k_{c}

where :math:`r_{c}` is the Wigner-Seitz (simulation cell image) radius,
and :math:`k_{c}` is the
length of the maximum :math:`k`-vector used in the long-ranged term.
Larger values of increase the accuracy of the evaluation.
A value of 15 tends to be conservative for the ``opt_breakup`` handler in 3D.

.. _particleset:

Specifying the particle set
---------------------------

The :code:`particleset` blocks specify the particles in the QMC simulations: their types,
attributes (mass, charge, valence), and positions.

Input specification
~~~~~~~~~~~~~~~~~~~

``particleset`` element:

  +-----------------+-----------------------+
  | Parent elements | ``simulation``        |
  +-----------------+-----------------------+
  | Child elements  | ``group``, ``attrib`` |
  +-----------------+-----------------------+

Attribute:

  +----------------------------------------+----------+----------------------+---------+-------------------------------+
  | Name                                   | Datatype | Values               | Default | Description                   |
  +========================================+==========+======================+=========+===============================+
  | ``name/id``                            | Text     | *Any*                | e       | Name of particle set          |
  +----------------------------------------+----------+----------------------+---------+-------------------------------+
  | ``size``:math:`^o`                     | Integer  | *Any*                | 0       | Number of particles in set    |
  +----------------------------------------+----------+----------------------+---------+-------------------------------+
  | ``random``\ :math:`^o`                 | Text     | Yes/no               | No      | Randomize starting positions  |
  +----------------------------------------+----------+----------------------+---------+-------------------------------+
  | ``randomsrc``/``randomsrc``:math:`^o`  | Text     | ``particleset.name`` | *None*  | Particle set to randomize     |
  +----------------------------------------+----------+----------------------+---------+-------------------------------+
  | ``spinor``:math:`^o`                   | Text     | Yes/no               | No      | particleset treated as spinor |
  +----------------------------------------+----------+----------------------+---------+-------------------------------+

Detailed attribute description
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Required particleset attributes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  | ``name/id``
   | Unique name for the particle set. Default is “e" for electrons. “i"
     or “ion0" is typically used for ions. For special cases where an
     empty particle set is needed, the special name “empty" can be used
     to bypass the zero-size error check.

Optional particleset attributes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  | ``size``
   | Number of particles in set.

``Group`` element:

  +-----------------+---------------------------+
  | Parent elements | ``particleset``           |
  +-----------------+---------------------------+
  | Child elements  | ``parameter``, ``attrib`` |
  +-----------------+---------------------------+

  Attribute:

  +---------------------+----------+--------+---------+-----------------------------+
  | Name                | Datatype | Values | Default | Description                 |
  +=====================+==========+========+=========+=============================+
  | ``name``            | Text     | *Any*  | e       | Name of particle set        |
  +---------------------+----------+--------+---------+-----------------------------+
  | ``size``:math:`^o`  | Integer  | *Any*  | 0       | Number of particles in set  |
  +---------------------+----------+--------+---------+-----------------------------+
  | ``mass``:math:`^o`  | Real     | *Any*  | 1       | Mass of particles in set    |
  +---------------------+----------+--------+---------+-----------------------------+
  | ``unit``:math:`^o`  | Text     | au/amu | au      | Units for mass of particles |
  +---------------------+----------+--------+---------+-----------------------------+

  Parameters:

  +------------------+----------+--------+---------+------------------------------------+
  | Name             | Datatype | Values | Default | Description                        |
  +==================+==========+========+=========+====================================+
  | ``charge``       | Real     | *Any*  | 0       | Charge of particles in set         |
  +------------------+----------+--------+---------+------------------------------------+
  | ``valence``      | Real     | *Any*  | 0       | Valence charge of particles in set |
  +------------------+----------+--------+---------+------------------------------------+
  | ``atomicnumber`` | Integer  | *Any*  | 0       | Atomic number of particles in set  |
  +------------------+----------+--------+---------+------------------------------------+

``attrib`` element:

  +---------------------+------------------------------------+
  | Parent elements     | ``particleset``, ``group``         |
  +---------------------+------------------------------------+

  Attribute:

  +--------------------+--------------+--------------------------------------------+-------------+------------------------+
  | **Name**           | **Datatype** | **Values**                                 | **Default** | **Description**        |
  +====================+==============+============================================+=============+========================+
  | ``name``           | String       | *Any*                                      | *None*      | Name of attrib         |
  +--------------------+--------------+--------------------------------------------+-------------+------------------------+
  | ``datatype``       | String       | IntArray, realArray, posArray, stringArray | *None*      | Type of data in attrib |
  +--------------------+--------------+--------------------------------------------+-------------+------------------------+
  | ``size``:math:`^o` | String       | *Any*                                      | *None*      | Size of data in attrib |
  +--------------------+--------------+--------------------------------------------+-------------+------------------------+

-  | ``random``
   | Randomize starting positions of particles. Each component of each
     particle’s position is randomized independently in the range of the
     simulation cell in that component’s direction.

-  | ``randomsrc``/``random_source``
   | Specify source particle set around which to randomize the initial
     positions of this particle set.

-  | ``spinor``
   | Sets an internal flag that the particleset (usually for electrons) is
     a spinor object. This is used in the wavefunction builders and QMC drivers
     to determiane if spin sampling will be used

Required name attributes
^^^^^^^^^^^^^^^^^^^^^^^^

-  | ``name``/``id``
   | Unique name for the particle set group. Typically, element symbols
     are used for ions and “u" or “d" for spin-up and spin-down electron
     groups, respectively.

Optional group attributes
^^^^^^^^^^^^^^^^^^^^^^^^^

-  | ``mass``
   | Mass of particles in set.

-  | ``unit``
   | Units for mass of particles in set (au[:math:`m_e` = 1] or
     amu[:math:`\frac{1}{12}m_{\rm ^{12}C}` = 1]).

Example use cases
~~~~~~~~~~~~~~~~~

.. _listing1:

.. centered:: Particleset elements for ions and electrons randomizing electron start positions.

::

     <particleset name="i" size="2">
       <group name="Li">
         <parameter name="charge">3.000000</parameter>
         <parameter name="valence">3.000000</parameter>
         <parameter name="atomicnumber">3.000000</parameter>
       </group>
       <group name="H">
         <parameter name="charge">1.000000</parameter>
         <parameter name="valence">1.000000</parameter>
         <parameter name="atomicnumber">1.000000</parameter>
       </group>
       <attrib name="position" datatype="posArray" condition="1">
       0.0   0.0   0.0
       0.5   0.5   0.5
       </attrib>
       <attrib name="ionid" datatype="stringArray">
          Li H
       </attrib>
     </particleset>
     <particleset name="e" random="yes" randomsrc="i">
       <group name="u" size="2">
         <parameter name="charge">-1</parameter>
       </group>
       <group name="d" size="2">
         <parameter name="charge">-1</parameter>
       </group>
     </particleset>

.. centered:: Particleset elements for ions and electrons specifying electron start positions.

::

     <particleset name="e">
       <group name="u" size="4">
         <parameter name="charge">-1</parameter>
         <attrib name="position" datatype="posArray">
       2.9151687332e-01 -6.5123272502e-01 -1.2188463918e-01
       5.8423636048e-01  4.2730406357e-01 -4.5964306231e-03
       3.5228575807e-01 -3.5027014639e-01  5.2644808295e-01
          -5.1686250912e-01 -1.6648002292e+00  6.5837023441e-01
         </attrib>
       </group>
       <group name="d" size="4">
         <parameter name="charge">-1</parameter>
         <attrib name="position" datatype="posArray">
       3.1443445436e-01  6.5068682609e-01 -4.0983449009e-02
          -3.8686061749e-01 -9.3744432997e-02 -6.0456005388e-01
       2.4978241724e-02 -3.2862514649e-02 -7.2266047173e-01
          -4.0352404772e-01  1.1927734805e+00  5.5610824921e-01
         </attrib>
       </group>
     </particleset>
     <particleset name="ion0" size="3">
       <group name="O">
         <parameter name="charge">6</parameter>
         <parameter name="valence">4</parameter>
         <parameter name="atomicnumber">8</parameter>
       </group>
       <group name="H">
         <parameter name="charge">1</parameter>
         <parameter name="valence">1</parameter>
         <parameter name="atomicnumber">1</parameter>
       </group>
       <attrib name="position" datatype="posArray">
         0.0000000000e+00  0.0000000000e+00  0.0000000000e+00
         0.0000000000e+00 -1.4308249289e+00  1.1078707576e+00
         0.0000000000e+00  1.4308249289e+00  1.1078707576e+00
       </attrib>
       <attrib name="ionid" datatype="stringArray">
         O H H
       </attrib>
     </particleset>

.. centered:: Particleset elements for ions specifying positions by ion type.

::

     <particleset name="ion0">
       <group name="O" size="1">
         <parameter name="charge">6</parameter>
         <parameter name="valence">4</parameter>
         <parameter name="atomicnumber">8</parameter>
         <attrib name="position" datatype="posArray">
           0.0000000000e+00  0.0000000000e+00  0.0000000000e+00
         </attrib>
       </group>
       <group name="H" size="2">
         <parameter name="charge">1</parameter>
         <parameter name="valence">1</parameter>
         <parameter name="atomicnumber">1</parameter>
         <attrib name="position" datatype="posArray">
           0.0000000000e+00 -1.4308249289e+00  1.1078707576e+00
           0.0000000000e+00  1.4308249289e+00  1.1078707576e+00
         </attrib>
       </group>
     </particleset>

.. bibliography:: /bibs/simulationcell.bib
