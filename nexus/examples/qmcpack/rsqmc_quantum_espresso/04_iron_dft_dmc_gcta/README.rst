Nexus QE+QMCPACK Example 4: Grand Canonical Twist Averaging in BCC Iron
===========================================================================
In this example, we consider twist averaging in a more complicated material: ferromagnetic body-centered cubic iron (BCC Fe).
It presents a challenge due to two reasons.
One is that it is a metal (conductor) and thus has a varying number of electron occupations at different twists due to band crossings at the Fermi level.
The second reason is that it is a ferromagnet, and we would like to closely reproduce the reference SCF cell magnetization in QMC twist averaging.

The pseudopotentials for this example are not included here but can be obtained from `pseudopotentiallibrary <https://pseudopotentiallibrary.org/>`_.
In the Nexus script ``iron_ldaU_dmc_gcta.py``, we carry out the standard sequence of SCF, NSCF, conversion, and Jastrow optimizations.
As usual, the Jastrow optimizations are carried out at the Gamma twist and with a spin close to the SCF value (5.66 Bohr mag/cell).

The main difference in this example is the ``gcta`` keyword in ``generate_qmcpack()`` in the final QMC run:

.. code-block:: python

    qmc = generate_qmcpack(
        ...
        gcta = 'safl',
        ...
        )

This keyword activates the grand canonical twist averaging (GCTA) with varying twist occupations.
The ``safl`` argument stands for "spin-adapted Fermi level", which guarantees charge neutrality in the system and closely reproduces the reference magnetization.
This is achieved by appropriately determining the up and down Fermi levels to meet the above criteria.
A summary of GCTA preprocessing is written in the ``gcta_report.txt`` file located in the same path as the QMC run.
The top portion of this file shows the summary of the GCTA occupations:

.. code-block::

    SUMMARY FOR GCTA OCCUPATIONS:
    ==================================================
    GCTA Flavor:                    safl
    Fermi Level Up [eV]            17.0858966293881203
    Fermi Level Dn [eV]            17.4080907712283910
    Net Charge:                     0
    Net Charge / Prim Cell:         0.0000000000000000
    Net Magnetization / Prim Cell:  5.6574074074074074
    SCF Magnetization (Reference):  5.6565051775613684
    
    
     TWISTNUM  NELEC_UP  NELEC_DN   CHARGE     SPIN
    ==================================================
        0         19        11         2         8
        1         19        11         2         8
        2         19        11         2         8
        3         20        10         2        10
        4         19        11         2         8
        5         19        11         2         8
        6         19        11         2         8
        7         19        12         1         7
        8         19        12         1         7
        9         20        13        -1         7
        ...

The user is encouraged to inspect this file.
As noted, we see that the net charge is zero, and the net magnetization is very close to the SCF magnetization.
The following section in ``gcta_report.txt`` shows how the number of up and down electrons in the inputs has been modified for each twist number.
Note that now, each twist has a varying charge and spin. 

The use of ``gcta = 'safl'`` in the final QMC run results in a significant drop in energy.
Specifically, without this keyword, the twist averaged DMC energy is ``-123.9042(3) Ha/atom``, while with ``gcta = 'safl'`` one obtains ``-123.9313(3) Ha/atom``.

An alternative way to run GCTA is to use ``gcta = 'afl'``, which uses an "adapted Fermi level" that guarantees charge neutrality but does not target the SCF magnetization.
This can be useful in noncollinear calculations where a single total magnetization is not defined.
Detailed description and performance of ``'safl'`` and ``'afl'`` can be found in `this article <https://pubs.acs.org/doi/10.1021/acs.jctc.4c00058>`_.
Also see the Nexus documentation for the details of the ``gcta`` implementation.
