.. _spin-orbit:

Spin-Orbit Calculations in QMC
==============================

Introduction
------------

In order to introduce relativistic effects in real materials, 
in principle the full Dirac equation must be solved where the resulting wave function is a four-component spinor. 
For the valence electrons that participate in chemistry, 
the single particle spinors can be well approximated by two-component spinors as two of the components are negligible. Note that this is not true for the deeper core electrons, where all four components contribute. 
In light of this fact, relativistic pseudopotentials have been developed to remove the core electrons while providing an effective potential for the valence electrons :cite:`Dolg2012`.
This allows relativistic effects to be studied in QMC methods using two-component spinor wave functions.

In QMCPACK, spin-orbit interactions have been implemented following the methodology described in :cite:`Melton2016-1`
and :cite:`Melton2016-2`.
We briefly describe some of the details below.

Single-Particle Spinors
-----------------------
The single particle spinors used in QMCPACK take the form 

.. math::
  :label: seqn1

    \phi(\mathbf{r},s) = \phi^\uparrow(\mathbf{r}) \chi^\uparrow(s) + \phi^{\downarrow}(\mathbf{r})\chi^\downarrow(s) \\
                       =  \phi^\uparrow(\mathbf{r}) e^{i s} + \phi^{\downarrow}(\mathbf{r}) e^{-i s}

where :math:`s` is the spin variable and using the complex spin representation.
In order to carry out spin-orbit calculations in solids, the single-particle spinors
can be obtained using Quantum ESPRESSO. After carrying out the spin-orbit calculation in QE
(with flags ``noncolin`` = .true., ``lspinorb`` = .true., and a relativistic ``.UPF`` pseudopotential), 
the spinors can be obtained by using the converter *convertpw4qmc*:

::

    convertpw4qmc data-file-schema.xml

where the ``data-file-schema.xml`` file is output from your QE calculation. 
This will produce an ``eshdf.h5`` file which contains the up and down components of the spinors per k-point.

Trial Wavefunction
------------------
Using the generated single particle spinors, we build the many-body wavefunction in a similar fashion to the normal non-relativistic calculations, namely

.. math::
  :label: seqn2

  \Psi_T(\mathbf{R},\mathbf{S}) = e^J \sum\limits_\alpha c_\alpha \det_\alpha\left[ \phi_i(\mathbf{r}_j, s_j) \right]\:,

where we now utilize determinants of spinors, as opposed to the usual product of up and down determinants. An example xml input block for the trial wave function is show below:

.. code-block::
  :caption: wavefunction specification for a single determinant trial wave function
  :name: slisting1

  <?xml version="1.0"?>
  <qmcsystem>
    <wavefunction name="psi0" target="e">
      <sposet_builder name="spo_builder" type="bspline" href="eshdf.h5" tilematrix="100010001" twistnum="0" source="ion0" size="10">
        <sposet type="bspline" name="myspo" size="10">
          <occupation mode="ground"/>
        </sposet>
      </sposet_builder>
      <determinantset>
        <slaterdeterminant>
          <determinant id="det" group="u" sposet="myspo" size="10"/>
        </slaterdeterminant>
      </determinantset>
      <jastrow type="One-Body" name="J1" function="bspline" source="ion0" print="yes">
          <correlation elementType="O" size="8" cusp="0.0">
             <coefficients id="eO" type="Array">                  
             </coefficients>
          </correlation>
       </jastrow> 
      <jastrow type="Two-Body" name="J2" function="bspline" print="yes">
        <correlation speciesA="u" speciesB="u" size="8">
          <coefficients id="uu" type="Array">                  
          </coefficients>
        </correlation> 
      </jastrow> 
    </wavefunction>
  </qmcsystem>
  
We note that we only specify an "up" determinant, since we no longer 
need a product of up and down determinants.
In the Jastrow specification, we only need to provide 
the jastrow terms for the same spin as there is no longer a
distinction between the up and down spins. 

We also make a small modification in the particleset specification:

.. code-block::
  :caption: specification for the electron particle when performing spin-orbit calculations
  :name: slisting2

  <particleset name="e" random="yes" randomsrc="ion0" spinor="yes">
     <group name="u" size="10" mass="1.0">
        <parameter name="charge"              >    -1                    </parameter>
        <parameter name="mass"                >    1.0                   </parameter>
     </group>
  </particleset>

Note that we only provide a single electron group to represent all electrons  in the system, as opposed to the usual separation of up and down electrons. 
The additional keyword ``spinor=yes`` is the *only* required keyword for spinors.
This will be used internally to determine which movers to use in QMC drivers (e.g. VMCUpdatePbyP vs SOVMCUpdatePbyP) and which SPOSets to use in the trial wave function (spinors vs. normal orbitals)

*note*: In the current implementation, spinor wavefunctions are only 
supported at the single determinant level. Multideterminant spinor wave functions will be supported in a future release. 


QMC Methods
-----------
In this formalism, the spin degree of freedom becomes
a continuous variable similar to the spatial degrees of freedom.
In order to sample the spins, we introduce a *spin kinetic energy* operator 

.. math:: 
  :label: seqn3

  T_s = \sum_{i=1}^{N_e} -\frac{1}{2\mu_s} \left[ \frac{\partial^2}{\partial s_i^2} +  1\right]\:, 

where :math:`\mu_s` is a spin mass. This operator vanishes when acting on
an arbitrary spinor or anti-symmetric product of spinors due to the offset.
This avoids any unphysical contribution to the local energy. However, this does contribute to the Green's function in DMC, 

.. math::
  :label: seqn4

  G(\mathbf{R}' \mathbf{S}' \leftarrow \mathbf{R}\mathbf{S}; \tau, \mu_s) \propto G(\mathbf{R}'\leftarrow\mathbf{R}; \tau) \exp\left[ -\frac{\mu_s}{2\tau}\left| \mathbf{S}' - \mathbf{S} - \frac{\tau}{\mu_s}\mathbf{v}_{\mathbf{S}}(\mathbf{S})\right|^2\right] \:,

where :math:`G(\mathbf{R}'\leftarrow\mathbf{R}; \tau)` is the usual 
Green's function for the spatial evolution and the *spin kinetic energy* 
operator introduces a Green's function for the spin variables. 
Note that this includes a contribution from the *spin drift* :math:`\mathbf{v}_{\mathbf{S}}(\mathbf{S}) =  \nabla_{\mathbf{S}} \ln \Psi_T(\mathbf{S})`.

In both the VMC and DMC methods, there are no required changes to a typical input

.. code-block::

  <qmc method="vmc/dmc">
    <parameter name="steps"    >    50   </parameter>
    <parameter name="blocks"   >    50   </parameter>
    <parameter name="walkers"  >    10   </parameter>
    <parameter name="timestep" >  0.01   </parameter>
  </qmc>

Whether or not spin moves are used is determined internally by the ``spinor=yes`` flag in particleset.

By default, the spin mass :math:`\mu_s` (which controls the rate of spin sampling relative to the spatial sampling) is set to 1.0. 
This can be changed by adding an additional parameter to the QMC input

.. code-block::

 <parameter name="spinMass" > 0.25 </parameter>

A larger/smaller spin mass corresponds to slower/faster spin sampling relative to the spatial coordinates.

Spin-Orbit Effective Core Potentials
------------------------------------

The spin-orbit contribution to the Hamiltonian can be introduced 
through the use of Effective Core Potentials (ECPs). 
As described in :cite:`Melton2016-2`, the relativistic (semilocal) ECPs take the general form

.. math::
  :label: seqn5
  
  W^{\rm RECP} = W_{LJ}(r) + \sum_{\ell j m_j} W_{\ell j}(r) | \ell j m_j \rangle \langle \ell j m_j |

where the projectors :math:`|\ell j m_j\rangle` are the so-called spin spherical harmonics. 
An equivalent formulation is to decouple the fully relativistic effective core potential (RECP) into *averaged relativistic* (ARECP)  and *spin-orbit* (SORECP) contributions:

.. math::
  :label: seqn6

  W^{\rm RECP} =  W^{\rm ARECP} + W^{\rm SOECP} \\
  W^{\rm ARECP} =  W^{\rm ARECP}_L(r) + \sum_{\ell m_\ell} W_\ell^{ARECP}(r) | \ell m_\ell \rangle \langle \ell m_\ell| \\
  W^{\rm SORECP} = \sum_\ell \frac{2}{2\ell + 1} \Delta W^{\rm SORECP}_\ell(r) \sum\limits_{m_\ell,m_\ell'} |\ell m_\ell \rangle \langle \ell m_\ell | \vec{\ell} \cdot \vec{s} | \ell m_\ell' \rangle \langle \ell m_\ell'|

Note that the :math:`W^{\rm ARECP}` takes exactly the same form as 
the semilocal pseudopotentials used in standard QMC calculations. 
In the pseudopotential ``.xml`` file format, the :math:`W^{\rm ARECP}_\ell(r)` terms are tabulated as usual.
If spin-orbit terms are included in the ``.xml`` file, the file must tabulate the entire radial spin-orbit prefactor :math:`\frac{2}{2\ell + 1}\Delta  W^{\rm SORECP}_\ell(r)`.
We note the following relations between the two representations of the relativistic potentials

.. math::
  :label: seqn7

  W^{\rm ARECP}_\ell(r) = \frac{\ell+1}{2\ell+1} W^{\rm RECP}_{\ell,j=\ell+1/2}(r) + \frac{\ell}{2\ell+1} W^{\rm RECP}_{\ell,j=\ell-1/2}(r) \\
  \Delta W^{\rm SORECP}_\ell(r) = W^{\rm RECP}_{\ell,j=\ell+1/2}(r) - W^{\rm RECP}_{\ell,j=\ell-1/2}(r)

The structure of the spin-orbit ``.xml`` is 

.. code-block::

  <?xml version="1.0" encoding="UTF-8"?>
  <pseudo>
    <header ... relativistic="yes" ... />
    <grid ... />
    <semilocal units="hartree" format="r*V" npots-down="4" npots-up="0" l-local="3" npots="2">
      <vps l="s" .../>
      <vps l="p" .../>
      <vps l="d" .../>
      <vps l="f" .../>
      <vps_so l="p" .../>
      <vps_so l="d" .../>
    </semilocal>
  </pseudo>

This is included in the Hamiltonian in the same way as the usual pseudopotentials. 
If the ``<vps_so>`` elements are found, the spin-orbit contributions will be present in the calculation. 
By default, the spin-orbit terms *will be* included in the local energy.
In order to accumulate the spin-orbit energy, but exclude it from the local energy (and therefore will not be propogated into the walker weights in DMC for example),
the ``physicalSO`` flag should be set to no in the Hamiltonian input.
A typical application will include the SOC terms in the local energy, and an example input block is given as

.. code-block::
  
  <hamiltonian name="h0" type="generic" target="e">
    <pairpot name="ElecElec" type="coulomb" source="e" target="e" physical="true"/>
    <pairpot name="IonIon" type="coulomb" source=ion0" target="ion0" physical="true"/>
    <pairpot name="PseudoPot" type="pseudo" source="i" wavefunction="psi0" format="xml" algorithm="non-batched">
      <pseudo elementType="Pb" href="Pb.xml"/>
    </pairpot>
  </hamiltonian>

The contribution from the spin-orbit will be printed to the ``.stat.h5`` and ``.scalar.dat`` files for post-processing.
An example output is shown below

::

  LocalEnergy           =           -3.4419 +/-           0.0014
  Variance              =            0.1132 +/-           0.0013
  Kinetic               =            1.1252 +/-           0.0027
  LocalPotential        =           -4.5671 +/-           0.0028
  ElecElec              =            1.6881 +/-           0.0025
  LocalECP              =           -6.5021 +/-           0.0062
  NonLocalECP           =            0.3286 +/-           0.0025
  LocalEnergy_sq        =           11.9601 +/-           0.0086
  SOECP                 =          -0.08163 +/-           0.0003

The ``NonLocalECP`` represents the :math:`W^{\rm ARECP}`, ``SOECP`` represents the :math:`W^{\rm SORECP}`, and the sum is the full :math:`W^{\rm RECP}` contribution.

Note that for now, the default "batched" non-local pseudopotential evaluation is not compatible with dynamical spin QMC calculations.  Therefore, the specification of algorithm="non-batched" in all pseudopotential blocks is required.  


.. bibliography:: /bibs/spin-orbit.bib
