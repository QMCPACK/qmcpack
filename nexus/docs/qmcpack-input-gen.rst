.. _qmcpack-input-gen:

Generating QMCPACK input files
==============================

Nexus provides many ways to generate QMCPACK input.  Documentation in this
section is intentionally brief at present and will be expanded in the future.
Currently covered features are:

* ``nrule``

The ``nrule`` input selects a spherical quadrature rule from 1 through 8 for
nonlocal pseudopotentials.  A single integer applies the same rule to every
atomic species:

.. code-block:: python

  qmc = generate_qmcpack(...,pseudos=['Fe.ccECP.xml'],nrule=7)

.. code-block:: xml

  <hamiltonian name="h0" type="generic" target="e">
     <pairpot type="coulomb" name="ElecElec" source="e" target="e"/>
     <pairpot type="coulomb" name="IonIon" source="ion0" target="ion0"/>
     <pairpot type="pseudo" name="PseudoPot" source="ion0" wavefunction="psi0" format="xml">
        <pseudo elementType="Fe" href="Fe.ccECP.xml" nrule="7"/>
     </pairpot>
  </hamiltonian>

A dictionary can instead select a different rule for each atomic species:

.. code-block:: python

  qmc = generate_qmcpack(
      ...,
      pseudos = ['V.ccECP.xml','O.ccECP.xml'],
      nrule   = {'V':6,'O':4},
      )

.. code-block:: xml

  <hamiltonian name="h0" type="generic" target="e">
     <pairpot type="coulomb" name="ElecElec" source="e" target="e"/>
     <pairpot type="coulomb" name="IonIon" source="ion0" target="ion0"/>
     <pairpot type="pseudo" name="PseudoPot" source="ion0" wavefunction="psi0" format="xml">
        <pseudo elementType="O" href="O.ccECP.xml" nrule="4"/>
        <pseudo elementType="V" href="V.ccECP.xml" nrule="6"/>
     </pairpot>
  </hamiltonian>
