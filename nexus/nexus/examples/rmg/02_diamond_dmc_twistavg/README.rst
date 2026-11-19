Nexus RMG+QMCPACK Example 2: Diamond cubic cell (twist averaged DMC)
=====================================================================

This example runs a complete RMG to QMCPACK workflow for an 8-atom cell of
diamond with a 2x2x2 twist grid:

1. an RMG SCF run in the primitive cell on a 4x4x4 k-point mesh,
2. an RMG NSCF run on the k-points of the twists, which writes the orbitals
   for QMCPACK (``write_qmcpack_restart = True``),
3. a QMCPACK optimization of one- and two-body Jastrow factors at the first
   twist, and
4. a twist-averaged QMCPACK DMC run.

The script is ``diamond_rmg_dmc_twistavg.py``.  The NSCF run depends on the
``charge_density`` of the SCF run, and the QMCPACK runs depend on the
``orbitals`` of the NSCF run.  No converter step is needed because RMG writes
the orbital file (``Waves/wave.out.h5``) itself:

.. parsed-literal::

    nscf = generate_rmg(
        ...
        **calculation_mode      = 'NSCF',**
        **write_qmcpack_restart = True,**
        **dependencies          = (scf,'charge_density'),**
        ...
        )

    qmc = generate_qmcpack(
        ...
        **dependencies = [(nscf,'orbitals'),**
                        **(opt,'jastrow')],**
        )

The NSCF run restarts from the restart files of the SCF run, so both runs
use the same RMG job, wavefunction grid and ``kpoint_distribution``.  Do not
set ``kpoint_mesh`` in the NSCF run: Nexus then writes the k-points of the
twists into the RMG input, and every twist is present in the orbital file.

The script expects ``rmg-cpu`` (RMG built with HDF5 and QMCPACK support) and
``qmcpack`` on your ``PATH``, and uses the pseudopotentials in
``../../qmcpack/pseudopotentials``.
