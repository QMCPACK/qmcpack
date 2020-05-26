.. _chap:simulationcell:

Specifying the Simulation Cell
==============================

The :code:`simulationcell` block specifies the geometry of the cell, how the boundary
conditions should be handled, and how ewald summation should be broken
up.

``simulationcell`` Element:

+------------------+-------------------------------------------------------------------------------------------------------+
| Parent elements: | :code:`qmcsystem`                                                                                     |
+------------------+-------------------------------------------------------------------------------------------------------+
| Child elements:  | None                                                                                                  |
+------------------+-------------------------------------------------------------------------------------------------------+

Attributes:

+---------------------+--------------+---------------------------+-------------------+----------------------------------------------------+
| **parameter name**  | **datatype** | **values**                | **default**       | **description**                                    |
+=====================+==============+===========================+===================+====================================================+
| ``lattice``         | 9 floats     | any float                 | Must be specified | Specification of lattice vectors.                  |
+---------------------+--------------+---------------------------+-------------------+----------------------------------------------------+
| ``bconds``          | string       | "p" or "n"                | "n n n "          | Boundary conditions for each axis.                 |
+---------------------+--------------+---------------------------+-------------------+----------------------------------------------------+
| ``vacuum``          | float        | :math:`\geq 1.0`          | 1.0               | Vacuum scale.                                      |
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
-------

The cell is specified using 3 lattice vectors.

Boundary conditions
-------------------

QMCPACK offers the capability to use a mixture of open and periodic
boundary conditions. The parameter expects a single string of three
characters separated by spaces, *e.g.* “p p p” for purely periodic
boundary conditions. These characters control the behavior of the
:math:`x`, :math:`y`, and :math:`z`, axes, respectively. Non periodic
directions must be placed after the periodic ones. Examples of valid
include:

**“p p p”** Periodic boundary conditions. Corresponds to a 3D crystal.

**“p p n”** Slab geometry. Corresponds to a 2D crystal.

**“p n n”** Wire geometry. Corresponds to a 1D crystal.

**“n n n”**
Open boundary conditions. Corresponds to an isolated molecule in a vacuum.

Vacuum
------

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
the axis labeled by a factor of``vacuum``. Note that all the particles remain in
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
     </simulationcell>

LR_dim_cutoff
-------------

When using periodic boundary conditions direct calculation of the
Coulomb energy is not well behaved. As a result, QMCPACK uses an
optimized Ewald summation technique to compute the Coulomb
interaction. :cite:`Natoli1995`

In the Ewald summation, the energy is broken into short- and long-ranged
terms. The short-ranged term is computed directly in real space, while
the long-ranged term is computed in reciprocal space. controls where the
short-ranged term ends and the long-ranged term begins. The real-space
cutoff, reciprocal-space cutoff, and are related via:

.. math:: \mathrm{LR\_dim\_cutoff} = r_{c} \times k_{c}

where :math:`r_{c}` is the Wigner-Seitz radius, and :math:`k_{c}` is the
length of the maximum :math:`k`-vector used in the long-ranged term.
Larger values of increase the accuracy of the evaluation. A value of 15
tends to be conservative.

.. bibliography:: bibliography.bib
