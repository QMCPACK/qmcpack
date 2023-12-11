.. _theory:

QMC Practice in a Nutshell
==========================

The aim of this section is to provide a very brief overview of the
essential concepts undergirding Quantum Monte Carlo calculations
of electronic structure with a particular focus on the key
approximations and quantities to converge to achieve high
accuracy.  The discussion here is not intended to be comprehensive.
For deeper perspectives on QMC, please see the review articles
listed in :ref:`learn-qmc`.

VMC and DMC in the abstract
---------------------------

Ground state QMC methods, such as Variational (VMC) and Diffusion (DMC) Monte
Carlo, attempt to obtain the ground state energy of a many-body quantum system.

.. math::
  :label: eq1

   \begin{aligned}
     E_0 = \frac{\langle{\Psi_0}|\hat{H}|{\Psi_0}\rangle}{\langle{\Psi_0}|{\Psi_0}\rangle}\end{aligned}

The VMC method obtains an upper bound on the ground state energy
(guaranteed by the Variational Principle) by introducing a guess at the
ground state wavefunction, known as the trial wavefunction
:math:`\Psi_T`:

.. math::
  :label: eq2

   \begin{aligned}
     E_{VMC} = \frac{\langle{\Psi_T}|\hat{H}|{\Psi_T}\rangle}{\langle{\Psi_T}|{\Psi_T}\rangle} \ge E_0\end{aligned}

The DMC method improves on this variational bound by projecting out
component eigenstates of the trial wavefunction lying higher in energy
than the ground state. The operator that acts as a projector is the
imaginary time, or thermodynamic, density matrix:

.. math::
  :label: eq3

   \begin{aligned}
     |{\Psi_t}\rangle & = e^{-t\hat{H}}|{\Psi_T} \nonumber \\
                  & = e^{-tE_0}\left( |{\Psi_0}\rangle + \sum_{n>0}e^{-t(E_n-E_0)}|{\Psi_n} \rangle \right) \nonumber \\
                  & \xrightarrow [t\rightarrow\infty]{} e^{-tE_0} |{\Psi_0}\rangle
                  \end{aligned}


The DMC energy approaches the ground state energy from above as the
imaginary time becomes large.

.. math::
  :label: eq4

   \begin{aligned}
     E_{DMC}  & = \lim_{t\rightarrow\infty}\frac{\langle{\Psi_t}|\hat{H}|{\Psi_t}\rangle}{\langle{\Psi_t}|{\Psi_t}\rangle} = E_0\end{aligned}

However from the equations above, one can already anticipate that the
DMC method will struggle in the face of degeneracy or near-degeneracy.

In principle, the DMC method is exact for the ground state, but further
complications arise for systems that are extended, comprised of
fermions, or contain heavy nuclei, pseudized or otherwise.
Approximations arising from the numerical implementation of the method
also require care to keep under control.

From expectation values to random walks
---------------------------------------

Evaluating expectation values of a many-body system involves performing
high dimensional integrals (the dimensionality is at least the
dimensions of the physical space times the number of particles). In VMC,
for example, the expectation value of the total energy is represented
succinctly as:

.. math::
  :label: eq5

   \begin{aligned}
     E_{VMC} &= \int dR |{\Psi_T}|^2 E_L\end{aligned}

where :math:`E_L` is the local energy
:math:`E_L=\Psi_T^{-1}\hat{H}\Psi_T`. The other factor in the integral
:math:`|{\Psi_T}|^2` can clearly be thought of as a probability
distribution and can therefore be sampled by Monte Carlo methods (such
as the Metropolis algorithm) to evaluate the integral exactly.

The sampling procedure takes the form of random walks. A “walker” is
just a set of particle positions, along with a weight, that evolves (or
moves) to new positions according to a set of statistical rules. In VMC
as few walkers are used as possible to reduce the equilibration time
(the number of steps or moves required to lose a memory of the
potentially poor starting guess for particle positions). In DMC, the
walker population is a dynamic feature of the calculation and must be
large enough to avoid introducing bias in expectation values.

The tradeoff of moving to a the sampling procedure for the integration
is that it introduces statistical error into the calculation which
diminishes slowly with the number of samples (it falls off like
:math:`1/(\# of samples)` by the Central Limit Theorem). The good news
for ground state QMC is that this error can be reduced more rapidly
through the discovery of better guesses at the detailed nature of the
many-body wavefunction.

Quality orbitals: planewaves, cutoffs, splines, and meshes
----------------------------------------------------------

Acting on an understanding of perturbation theory, the zeroth order
representation of the wavefunction of an interacting system takes the
form of a Slater determinant of single particle orbitals. In practice,
QMC calculations often obtain a starting guess at these orbitals from
Hartree-Fock or Density Functional Theory calculations (which already
contain non-perturbative contributions from correlation). An important
factor in the generation and use of these orbitals is to ensure that
they are described to high accuracy within the parent theory.

For example, when taking orbitals from a planewave DFT calculation, one
must take care to converge the planewave energy cutoff to a sufficient
level of accuracy (usually far beyond what is required to obtain an
accurate DFT energy). One criterion to use it to converge the kinetic
energy of the Kohn-Sham wavefunction with respect to the planewave
energy cutoff until it is accurate to the energy scale you care about in
your production QMC calcuation. For systems with a small number of
valence electrons, a cutoff of around 200 Ry is often sufficient. To
obtain the kinetic energy from a PWSCF calculation the ``pw2casino.x``
post-processing tool can be used. In Nexus one has the option to compute
the kinetic energy by setting the ``kinetic_E`` flag in the
``standard_qmc`` or ``basic_qmc`` convenience functions.

For efficiency reasons, QMC codes often use a real-space representation
of the wavefunction. It is common to represent the orbitals in terms of
B-splines which have control points, or knots, that fall on a regular
3-D mesh. Analogous to the planewave cutoff, the fineness of the
B-spline mesh controls the quality of the represented orbitals. To
verify that the quality of the orbitals has not been compromised during
the conversion process from planewave to B-spline, one often performs a
VMC calculation with the B-spline Slater determinant wavefunction to
obtain the kinetic energy. This value should agree with the kinetic
energy of the planewave representation within the energy scale of
interest.

In QMCPACK, the B-spline mesh is controlled with the ``meshfactor``
keyword. Larger values correspond to finer meshes. A value of
:math:`1.0` usually gives a similar quality representation as the
original planewave calculation. Control of this parameter is made
available in Nexus through the ``meshfactor`` keyword in the
``standard_qmc`` or ``basic_qmc`` convenience functions.

Quality Jastrows: less variance = more efficient
------------------------------------------------

aking a further cue from perturbation theory, the first order
correction to the Slater determinant wavefunction is the Jastrow
correlation prefactor.

.. math::
  :label: eq6

   \begin{aligned}
     \Psi_T\approx e^{-J}\Psi_{Slater\, Det.}\end{aligned}

In a quantum liquid, an appropriate form for the Jastrow factor is:

.. math::
  :label: eq7

   \begin{aligned}
     J = \sum_{i<j} u_{ij}(|{r_i-r_j}|)\end{aligned}

This form is often used without modification in electronic structure
calculations. Note that the correlation factors :math:`u_{ij}` can be
different for particles of differing species, or, if one of the
particles in the pair is classical (such as a heavy atomic nucleus), the
local electronic environment varies across the system.

The primary role of the Jastrow factor is to increase the efficiency of
the QMC calculation. The variance of the local energy across all samples
of the random walk is directly related to the statistical error of the
final results:

.. math::
  :label: eq8

   \begin{aligned}
     v_{\Psi_T} &= \frac{1}{N_{samples}}\sum_{s\in samples} E_L(s)^2 - \left[\frac{1}{N_{samples}}\sum_{s\in samples} E_L(s)\right]^2 \\
     \sigma_{error} &\approx \sqrt{\frac{v_{\Psi_T}}{N_{samples}}}\end{aligned}

The variance of local energy is usually minimized by performing a
statistical optimization of the Jastrow factor with QMC.

In addition to selecting a good form for the pair correlation functions
:math:`u_{ij}` (which are represented in QMCPACK as 1-D B-spline
functions with a finite cutoff radius), the (iterative) optimization
procedure must be performed with a sufficient number of samples to
converge all the free parameters. Starting with a small number of
samples (:math:`\approx 20,000`) is usually preferable for early
iterations, followed by a larger number for later iterations. This
larger number is something close to
:math:`100,000\times (\#~of~free~parameters)^2`. For B-spline functions,
the number of free parameters is the number of control points, or knots.

The number of samples is controlled with the ``samples`` keyword in
QMCPACK. Control of this parameter is made available in Nexus through
the ``samples`` keyword in the ``linear`` or ``cslinear`` convenience
functions (Which are often used in conjunction with ``standard_qmc`` or
``basic_qmc``). For a B-spline correlation factor, the number of free
parameters/knots is indicated by the ``size`` keyword in either QMCPACK
or Nexus.

Finite size effects: k-points, supercells, and corrections
----------------------------------------------------------

For extended systems, finite size errors are a key consideration. In
addition to the finite size effects that are typically seen in DFT
(k-points related). Correlated, many-body methods such as QMC also must
contend with correlation-related finite size effects. Both types of
finite-size effects are reduced by simply using larger supercells. The
complete elimination of finite size effects using this approach can be
prohibitively costly since the finite size error typically falls off
like :math:`1/\Omega_C`, where :math:`\Omega_C` is the volume of the
supercell. A more sophisticated approach involves a combination of the
supercell size, k-point grid, and additional estimated corrections for
correlation finite size effects.

Although there is no firm rule on the selection of these three elements,
adhering to some general guidelines is usually helpful. For a production
calculation of an extended system, the minimum supercell size is around
50 atoms. The size of the supercell k-point grid can then be determined
by proxy with a DFT calculation (converge the energy down to the scale
of interest). Note that although the cost of a DFT calculation scales
linearly with the number of k-points, the cost of the corresponding QMC
calculation is hardly increased due to the statistical averaging of the
results (the QMC calculation at each separate supercell k-point is
simply performed with fewer samples so that the total number of samples
remains fixed w.r.t. the number of k-points). Finally, corrections for
correlation-related finite size effects are computed during the QMC run
and added to the result by hand in post-processing the data.

In Nexus, the supercell size is controlled through the ``tiling``
parameter in the ``generate_physical_system``, ``generate_structure``,
``Structure``, or ``Crystal`` convenience functions. Supercells can also
be constructed by tiling existing structures through the ``tile`` member
function of ``Structure`` or ``PhysicalSystem`` objects. The k-point
grid is controlled through the ``kgrid`` parameter in the
``generate_physical_system``, ``generate_structure``, ``Structure``, or
``Crystal`` convenience functions. K-point grids can also be added to
existing structures through the ``add_kmesh`` member function of
``Structure`` or ``PhysicalSystem`` objects.

Imaginary time discretization: the DMC timestep
-----------------------------------------------

An analytic form for the imaginary time projection operator is not
known, but real-space approximations to it can be obtained in the small
time limit. With importance sampling included (not covered here), the
short-time projector splits into two parts, known as the drift-diffusion
and branching factors (shown below in atomic units):

.. math::
  :label: eq9

   \begin{aligned}
      \rho(R',R;t)  &= \langle{R'}|{\hat{\Psi_T}e^{-t\hat{H}}\hat{\Psi_T}^{-1}}|{R}\rangle \\
       &= G_d(R',R;t)G_b(R',R,t) +\mathcal{O}(t^2) \\
     G_d(R',R;t) &\equiv \exp{\left(-\tfrac{1}{2t}\left[R'-R-t\nabla_R\log\Psi_T(R)\right]^2\right)} \\
     G_b(R',R;t) &\equiv \exp{\left(\tfrac{1}{2}\left[E_L(R')+E_L(R)\right]\right)}\end{aligned}

The long-time projector is found as the product of many approximate
short-time solutions, which takes the form of a many-body path integral
in real space:

.. math::
  :label: eq10

   \begin{aligned}
     \rho(R_M,R_0; M\tau) = \int dR_1dR_{M-1}\ldots \prod_{m=0}^{M-1}\rho(R_{m+1},R_m;\tau)\end{aligned}

The short-time parameter :math:`\tau` is known as the DMC timestep and
accurate quantities are obtained only in the limit as :math:`\tau`
approaches zero.

Ensuring that the timestep error is sufficiently small usually involves
performing many DMC calculations over a range of timesteps (sometimes on
a smaller supercell than the production calculation). The largest
timestep is chosen that produces a bias smaller than the energy scale of
interest. For very high accuracy, one uses the total energy as a
function of timestep to extrapolate to the zero time limit.

The DMC timestep is made available in Nexus through the ``timestep``
parameter of the ``dmc`` convenience function (which is often used in
conjunction with the ``standard_qmc``, ``basic_qmc``,
``generate_qmcpack``, or ``Qmcpack`` functions).

Population control bias: safety in numbers
------------------------------------------

While the drift-diffusion factor :math:`G_d(R',R;\tau)` can be sampled
exactly using Gaussian distributed random numbers (this generates the
DMC random walk), the branching factor :math:`G_b(R',R;\tau)` is handled
a different way for efficiency. The product of branching factors over an
imaginary time trajectory (random walk) serves as a statistical weight
for each walker. The fluctuations in this weight rapidly become quite
large as the random walk progresses (because it approaches an infinite
product of real numbers). As its name suggests, this weight factor is
used to “branch” walkers every few steps. If the weight is small the
walker is deleted, but if the weight is large the walker is copied many
times (“branched”) with each copy carrying a weight close to unity. This
is more efficient because more walkers are created (and thus more
statistics are gathered) in the high weight regions of phase space that
contribute most to the integral.

The branching process in DMC naturally leads to a fluctuating population
of walkers. The fluctuations in the walker population, if left to its
own dynamics, are unbounded. This means that the walker population can
grow very large, or even become zero. To prevent collapse of the walker
population, population control techniques (not covered here) are added
to the algorithm. The practical upshot of population control is that it
introduces a systematic bias in the DMC results that scales like
:math:`1/(\# of walkers)` (Although note that another route to reduce
the population control bias is to improve the trial wavefunction, since
the fluctuations in the branching weights will become zero for the exact
ground state).

For many production calculations, population control bias is not much of
an issue because the simulations are performed on supercomputers with
thousands of cores per run, and thus tens of thousands of walkers. As a
rule of thumb, the walker population should at least number in the
thousands. One should occasionally explicitly check the magnitude of the
population control bias for the system under study since predictions
have been made that it will eventually diverge exponentially with the
number of particles in the system.

The DMC walker population can be directly controlled in QMCPACK or Nexus
through the ``samples`` (total walker population) or
``samplesperthread`` (walkers per OpenMP thread) keywords in the VMC
block directly proceeding DMC (``vmc`` convenience function in Nexus).
If you opt to use the ``samples`` keyword, check that each thread in the
calculation will have at least a few walkers.

The fixed node/phase approximation: varying the nodes/phase
-----------------------------------------------------------

For every fermionic system, the bosonic ground state lies lower in
energy than the fermionic ground state. This means that projection
methods like DMC will approach the bosonic ground state exponentially
fast in imaginary time if unconstrained (this would show up as an
exponentially diverging statistical error). In order to guarantee that
the projected wavefunction remains in the space of fermionic functions
(and consequently that the projected energy remains an upper bound to
the fermionic ground state energy), the projected wavefunction is
constrained to share the nodes (if it is real-valued) or the phase (if
it is complex-valued) of the trial wavefunction. The fixed node/phase
approximation represents one of the two most important approximations
for electronic structure calculations (the other is the pseudopotential
approximation covered in the next section).

The fixed node/phase error can be reduced, but it cannot be completely
eliminated unless the exact nodes/phase is known. A common approach to
reduce the fixed node/phase error is to perform several DMC calculations
(sometimes on a smaller supercell) with different sets of orbitals
(perhaps generated with different functionals). Another, more expensive
approach, is to include the backflow transformation (this is the second
order correction to the wavefunction; it is not covered in any detail
here) to get a lower bound on how large the fixed node error is in
standard Slater-Jastrow calculations.

To perform a calculation of this type (scanning over orbitals from
different functionals) with Nexus, the DFT functional can be selected
with the ``functional`` keyword in the ``standard_qmc`` or ``basic_qmc``
convenience functions. If you are using pseudopotentials generated for
use in DFT, you should maintain consistency between the functional and
pseudopotential. Even if such consistency is maintained, the impact of
using DFT pseudopotentials (or those made with many other theories) in
QMC can be significant.

Pseudopotentials: theoretical dissonance, the locality approximation, and T-moves
---------------------------------------------------------------------------------

The accurate use of pseudopotentials in electronic structure QMC
calculations remains one of the largest challenges in current practice.
The necessity for pseudopotentials arises from the rapidly increasing
computational cost with increasing nuclear charge (it scales like
:math:`Z^6`, compared with the :math:`N_{electrons}^3` scaling with
:math:`Z` fixed). The challenge in using pseudopotentials in QMC is that
practically no pseudopotentials exist that have been generated
self-consistently with QMC. In other words, QMC is currently reliant on
other theories to provide the pseudopotentials, which can be a critical
source of error.

The current state-of-the-art is not without rigor, however. One source
of Dirac-Fock based pseudopotentials, the Burkatzki-Filippi-Dolg
database (see http://www.burkatzki.com/pseudos/index.2.html), has been
explicitly vetted against quantum chemistry calculations of atoms (a
higher-fidelity proxy for QMC calculations of small systems). It must be
stressed that these pseudopotentials should still be validated for use
in a particular target system. Another collection of Dirac-Fock
pseudopotentials that have been created for use in QMC can be found in
the Trail-Needs database (see
http://www.tcm.phy.cam.ac.uk/~mdt26/casino2_pseudopotentials.html). Many
current calculations also use the OPIUM package (see
http://opium.sourceforge.net/) to generate DFT pseudopotentials and then
port them directly to QMC.

Whatever the source of pseudopotentials (but perhaps especially so for
those derived from DFT), testing and validation remains an important
step preceding production calculations. One option is to perform
parallel pseudopotential and all-electron DMC calculations of atoms with
varying electron count (*i.e.* ionization potential/electron affinity
calculations). As with any electronic structure calculation, it is also
advisable to devise a test in or close to the target host environment.
Validating pseudopotentials remains a difficult task, and while the
suggestions presented here may be of some help, they do not amount to a
panacea for the issue.

Beyond the central approximation of using a pseudopotential at all, two
approximations unique to pseudopotential use in DMC merit discussion.
The direct use of non-local pseudopotentials in DMC leads to a second
sign-problem (akin to the fixed-node issue) in the imaginary time
projector. One solution, devised first, is known as the locality
approximation. In the locality approximation, the non-local
pseudopotential is replaced by a “localized” form:
:math:`V_{NLPP}\rightarrow \Psi_T^{-1}V_{NLPP}\Psi_T`. This
approximation becomes exact as the trial wavefunction approaches the
pseudo ground state, however the Variational Principle of the
pseudo-system is lost (though it should be acknowledged that a
non-variational portion of the energy has been discarded by using
pseudopotentials at all). The Variational Principle for the
pseudo-system can be restored with an advanced sampling technique known
as T-moves (although the first incarnation of the technique reduces to
the locality approximation as the system becomes larger than several
atoms, the second version fixes this oversight).

One can select whether to use the locality approximation or T-moves
(version 1!) in QMCPACK from within Nexus by setting the parameter
``nonlocalmoves`` to True or False in the ``dmc`` convenience function.

Other approximations: what else is missing?
-------------------------------------------

Though a few points could be selected for mention at this point, only one
additional approximation will be highlighted here.  In most modern QMC
calculations of electronic structure, relativistic effects have been neglected
entirely (there have been a few exceptions) or simply assumed to be covered
by the pseudopotential.  Clearly this will become an issue for systems with
large effective core charges.  At present, relativistic corrections are not
available within QMCPACK.
