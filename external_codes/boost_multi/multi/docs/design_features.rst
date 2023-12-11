.. _design-features:

QMCPACK Design and Feature Documentation
========================================

This section contains information on the overall design of QMCPACK.  Also included are detailed explanations/derivations of major features and algorithms present in the code.

QMCPACK design
--------------

TBD.

Feature: Optimized long-range breakup (Ewald)
---------------------------------------------

.. _ Written by Ken Esler as part of the Common codebase used in wfconvert
   Originally titled "Ewald Breakup for Long-Range Potentials in PIMC"
   PIMC-specific portions have been commented out

Consider a group of particles interacting with long-range central
potentials, :math:`v^{\alpha \beta}(|r^{\alpha}_i - r^{\beta}_j|)`,
where the Greek superscripts represent the particle species (e.g.,
:math:`\alpha=\text{electron}`, :math:`\beta=\text{proton}`), and Roman
subscripts refer to particle number within a species. We can then write
the total interaction energy for the system as


.. math::
  :label: eq80

  V = \sum_\alpha \left\{\sum_{i<j} v^{\alpha\alpha}(|\mathbf{r}^\alpha_i - \mathbf{r}^\alpha_j|) +
  \sum_{\beta<\alpha}
  \sum_{i,j} v^{\alpha \beta}(|\mathbf{r}^{\alpha}_i - \mathbf{r}^{\beta}_j|) \right\}

The long-range problem
~~~~~~~~~~~~~~~~~~~~~~

Consider such a system in periodic boundary conditions in a cell defined
by primitive lattice vectors :math:`\mathbf{a}_1`, :math:`\mathbf{a}_2`,
and :math:`\mathbf{a}_3`. Let
:math:`\mathbf{L}\equiv n_1 \mathbf{a}_1 + n_2 \mathbf{a}_2 + n_3\mathbf{a}_3`
be a direct lattice vector. Then the interaction energy per cell for the
periodic system is given by

.. math::
  :label: eq81

  \begin{split}
  V = & \sum_\mathbf{L}\sum_\alpha \left\{
  \overbrace{\sum_{i<j} v^{\alpha\alpha}(|\mathbf{r}^\alpha_i - \mathbf{r}^\alpha_j + \mathbf{L}|)}^{\text{homologous}} +
  \overbrace{\sum_{\beta<\alpha}
  \sum_{i,j} v^{\alpha \beta}(|\mathbf{r}^{\alpha}_i - \mathbf{r}^{\beta}_j+\mathbf{L}|)}^{\text{heterologous}}
  \right\}  \\
  & + \underbrace{\sum_{\mathbf{L}\neq \mathbf{0}} \sum_\alpha N^\alpha v^{\alpha \alpha} (|\mathbf{L}|)}_\text{Madelung}\:.
  \end{split}

where :math:`N^\alpha` is the number particles of species
:math:`\alpha`. If the potentials :math:`v^{\alpha\beta}(r)` are indeed
long-range, the summation over direct lattice vectors will not converge
in this naive form. A solution to the problem was posited by Ewald. We
break the central potentials into two pieces—a short-range and a
long-range part defined by

.. math::
  :label: eq82

  v^{\alpha \beta}(r) = v_s^{\alpha\beta}(r) + v_l^{\alpha \beta}(r)\:.

We will perform the summation over images for the short-range part in
real space, while performing the sum for the long-range part in
reciprocal space. For simplicity, we choose
:math:`v^{\alpha \beta}_s(r)` so that it is identically zero at the
half-the-box length. This eliminates the need to sum over images in real
space.

Reciprocal-space sums
~~~~~~~~~~~~~~~~~~~~~

Heterologous terms
^^^^^^^^^^^^^^^^^^

We begin with :eq:`eq81`, starting with the heterologous terms (i.e., the terms involving particles of different species).  The
short-range terms are trivial, so we neglect them here.

.. math::
  :label: eq83

  \text{heterologous} = \frac{1}{2} \sum_{\alpha \neq \beta} \sum_{i,j} \sum_\mathbf{L}
   v^{\alpha\beta}_l(\mathbf{r}_i^\alpha - \mathbf{r}_j^\beta + \mathbf{L})\:.

We insert the resolution of unity in real space twice:

.. math::
  :label: eq84

  \begin{aligned}
   \text{heterologous} & = & \frac{1}{2}\sum_{\alpha \neq \beta} \int_\text{cell} d\mathbf{r}\, d\mathbf{r}' \, \sum_{i,j}
   \delta(\mathbf{r}_i^\alpha - \mathbf{r}) \delta(\mathbf{r}_j^\beta-\mathbf{r}') \sum_\mathbf{L}
   v^{\alpha\beta}_l(|\mathbf{r}- \mathbf{r}' + \mathbf{L}|)\:, \\
   & = & \frac{1}{2\Omega^2}\sum_{\alpha \neq \beta} \int_\text{cell} d\mathbf{r}\, d\mathbf{r}' \, \sum_{\mathbf{k}, \mathbf{k}', i, j} e^{i\mathbf{k}\cdot(\mathbf{r}_i^\alpha
     - \mathbf{r})} e^{i\mathbf{k}'\cdot(\mathbf{r}_j^\beta - \mathbf{r}')} \sum_\mathbf{L}
   v^{\alpha\beta}_l(|\mathbf{r}- \mathbf{r}' + \mathbf{L}|) \nonumber\:, \\
   & = & \frac{1}{2\Omega^2} \sum_{\alpha \neq \beta} \int_\text{cell} d\mathbf{r}\, d\mathbf{r}'\,
   \sum_{\mathbf{k}, \mathbf{k}', \mathbf{k}'', i, j} e^{i\mathbf{k}\cdot(\mathbf{r}_i^\alpha - \mathbf{r})}
   e^{i\mathbf{k}'\cdot(\mathbf{r}_j^\beta-\mathbf{r}')} e^{i\mathbf{k}''\cdot(\mathbf{r}-\mathbf{r}')}
   v^{\alpha\beta}_{\mathbf{k}''}\nonumber\:.\end{aligned}

Here, the :math:`\mathbf{k}` summations are over reciprocal lattice
vectors given by
:math:`\mathbf{k}= m_1 \mathbf{b}_1 + m_2\mathbf{b}_2 + m_3\mathbf{b}_3`,
where

.. math::
  :label: eq85

  \begin{aligned}
  \mathbf{b}_1 & = & 2\pi \frac{\mathbf{a}_2 \times \mathbf{a}_3}{\mathbf{a}_1 \cdot (\mathbf{a}_2 \times
    \mathbf{a}_3)} \nonumber\:, \\
  \mathbf{b}_2 & = & 2\pi \frac{\mathbf{a}_3 \times \mathbf{a}_1}{\mathbf{a}_1 \cdot (\mathbf{a}_2 \times
    \mathbf{a}_3)}\:, \\
  \mathbf{b}_3 & = & 2\pi \frac{\mathbf{a}_1 \times \mathbf{a}_2}{\mathbf{a}_1 \cdot (\mathbf{a}_2 \times
    \mathbf{a}_3)} \nonumber\:.\end{aligned}

We note that :math:`\mathbf{k}\cdot \mathbf{L}= 2\pi(n_1 m_1 + n_2 m_2 + n_3 m_3)`.

.. math::
  :label: eq86

  \begin{aligned}
  v_{k''}^{\alpha \beta} & = &
  \frac{1}{\Omega} \int_{\text{cell}} d\mathbf{r}'' \sum_\mathbf{L}
  e^{-i\mathbf{k}''\cdot(|\mathbf{r}''+\mathbf{L}|)} v^{\alpha\beta}(|\mathbf{r}''+\mathbf{L}|)\:, \\
  & = & \frac{1}{\Omega} \int_\text{all space} d\tilde{\mathbf{r}} \,
      e^{-i\mathbf{k}'' \cdot \tilde{\mathbf{r}}} v^{\alpha\beta}(\tilde{r})\:, \end{aligned}

where :math:`\Omega` is the volume of the cell. Here we have used the
fact that summing over all cells of the integral over the cell is
equivalent to integrating over all space.

.. math::
  :label: eq87

  \text{hetero} = \frac{1}{2\Omega^2} \sum_{\alpha \neq \beta}
  \int_\text{cell} d\mathbf{r}\, d\mathbf{r}' \, \sum_{\mathbf{k}, \mathbf{k}', \mathbf{k}'', i, j}
  e^{i(\mathbf{k}\cdot \mathbf{r}_i^\alpha + \mathbf{k}' \cdot\mathbf{r}_j^\beta)} e^{i(\mathbf{k}''-\mathbf{k})\cdot \mathbf{r}}
  e^{-i(\mathbf{k}'' + \mathbf{k}')\cdot \mathbf{r}'} v^{\alpha \beta}_{\mathbf{k}''}\:.

We have

.. math::
  :label: eq88

  \frac{1}{\Omega} \int d\mathbf{r}\  e^{i(\mathbf{k}-\mathbf{k}')\cdot \mathbf{r}} =
  \delta_{\mathbf{k},\mathbf{k}'}\:.

Then, performing the integrations we have

.. math::
  :label: eq89

  \begin{aligned}
  \text{hetero} = \frac{1}{2} \sum_{\alpha \neq \beta}
  \sum_{\mathbf{k}, \mathbf{k}', \mathbf{k}'', i, j}
  e^{i(\mathbf{k}\cdot \mathbf{r}_i^\alpha + \mathbf{k}' \cdot\mathbf{r}_j^\beta)} \delta_{\mathbf{k},\mathbf{k}''}
  \delta_{-\mathbf{k}', \mathbf{k}''} v^{\alpha \beta}_{\mathbf{k}''}\:.\end{aligned}

We now separate the summations, yielding

.. math::
  :label: eq90

  \text{hetero} = \frac{1}{2} \sum_{\alpha \neq \beta} \sum_{\mathbf{k}, \mathbf{k}'}
  \underbrace{\left[\sum_i e^{i\mathbf{k}\cdot \mathbf{r}_i^\alpha} \rule{0cm}{0.705cm}
      \right]}_{\rho_\mathbf{k}^\alpha}
  \underbrace{\left[\sum_j e^{i\mathbf{k}' \cdot \mathbf{r}_j^\beta} \right]}_{\rho_{\mathbf{k}'}^\beta}
   \delta_{\mathbf{k},\mathbf{k}''} \delta_{-\mathbf{k}', \mathbf{k}''} v^{\alpha
    \beta}_{\mathbf{k}''}\:.

Summing over :math:`\mathbf{k}` and :math:`\mathbf{k}'`, we have

.. math::
  :label: eq91

  \text{hetero} = \frac{1}{2} \sum_{\alpha \neq \beta} \sum_{\mathbf{k}''}
  \rho_{\mathbf{k}''}^\alpha \, \rho_{-\mathbf{k}''}^\beta v_{k''}^{\alpha \beta}\:.

We can simplify the calculation a bit further by rearranging the
sums over species:

.. math::
  :label: eq92

  \begin{aligned}
  \text{hetero} & = & \frac{1}{2} \sum_{\alpha > \beta} \sum_{\mathbf{k}}
  \left(\rho^\alpha_\mathbf{k}\rho^\beta_{-\mathbf{k}} + \rho^\alpha_{-\mathbf{k}}
  \rho^\beta_\mathbf{k}\right) v_{k}^{\alpha\beta}\:, \\
  & = & \sum_{\alpha > \beta} \sum_\mathbf{k}\mathcal{R}e\left(\rho_\mathbf{k}^\alpha
  \rho_{-\mathbf{k}}^\beta\right)v_k^{\alpha\beta} .\end{aligned}


Homologous terms
^^^^^^^^^^^^^^^^

We now consider the terms involving particles of the same species
interacting with each other.  The algebra is very similar to the
preceding, with the slight difficulty of avoiding the self-interaction term.

.. math::
  :label: eq93

  \begin{aligned}
  \text{homologous} & = & \sum_\alpha \sum_L \sum_{i<j} v_l^{\alpha
    \alpha}(|\mathbf{r}_i^\alpha - \mathbf{r}_j^\alpha + \mathbf{L}|)\:, \\
   & = & \frac{1}{2} \sum_\alpha \sum_L \sum_{i\neq j} v_l^{\alpha
    \alpha}(|\mathbf{r}_i^\alpha - \mathbf{r}_j^\alpha + \mathbf{L}|)\:. \end{aligned}

.. math::
  :label: eq94

  \begin{aligned}
  \text{homologous} & = & \frac{1}{2} \sum_\alpha \sum_L
  \left[
  -N^\alpha v_l^{\alpha \alpha}(|\mathbf{L}|)  + \sum_{i,j} v^{\alpha \alpha}_l(|\mathbf{r}_i^\alpha - \mathbf{r}_j^\alpha + \mathbf{L}|)
    \right]\:, \\
  & = & \frac{1}{2} \sum_\alpha \sum_\mathbf{k}\left(|\rho_k^\alpha|^2 - N
  \right) v_k^{\alpha \alpha}\:.\end{aligned}

Madelung terms
^^^^^^^^^^^^^^

Let us now consider the Madelung term for a single particle of species
:math:`\alpha`. This term corresponds to the interaction of a particle
with all of its periodic images.

.. math::
  :label: eq95

  \begin{aligned}
  v_M^{\alpha} & = & \frac{1}{2} \sum_{\mathbf{L}\neq \mathbf{0}} v^{\alpha
    \alpha}(|\mathbf{L}|)\:, \\
  & = & \frac{1}{2} \left[ -v_l^{\alpha \alpha}(0) + \sum_\mathbf{L}v^{\alpha
    \alpha}(|\mathbf{L}|) \right]\:, \\
  & = & \frac{1}{2} \left[ -v_l^{\alpha \alpha}(0) + \sum_\mathbf{k}v^{\alpha
    \alpha}_\mathbf{k}\right]\:.  \end{aligned}

:math:`\mathbf{k}=\mathbf{0}` terms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Thus far, we have neglected what happens at the special point
:math:`\mathbf{k}=
\mathbf{0}`. For many long-range potentials, such as the Coulomb
potential, :math:`v_k^{\alpha \alpha}` diverges for :math:`k=0`.
However, we recognize that for a charge-neutral system, the divergent
part of the terms cancel each other. If all the potential in the system
were precisely Coulomb, the :math:`\mathbf{k}=\mathbf{0}` terms would
cancel precisely, yielding zero. For systems involving PPs, however, it
may be that the resulting term is finite, but nonzero. Consider the
terms from :math:`\mathbf{k}=\mathbf{0}`:

.. math::
  :label: eq96

  \begin{aligned}
  V_{k=0} & = & \sum_{\alpha>\beta} N^\alpha N^\beta v^{\alpha \beta}_{k=0}
  + \frac{1}{2} \sum_\alpha \left(N^{\alpha}\right)^2 v^{\alpha\alpha}_{k=0}\:, \\
  & = & \frac{1}{2} \sum_{\alpha,\beta} N^\alpha N^\beta v^{\alpha
    \beta}_{k=0}\:.
  \end{aligned}

Next, we must compute :math:`v^{\alpha \beta}_{k=0}`.

.. math::
  :label: eq97

  v^{\alpha \beta}_{k=0} = \frac{4 \pi}{\Omega} \int_0^\infty dr\ r^2
  v_l^{\alpha \beta}(r)\:.

We recognize that this integral will not converge because of the
large-:math:`r` behavior. However, we recognize that when we do the sum
in :eq:`eq96`, the large-:math:`r` parts of the integrals
will cancel precisely. Therefore, we define

.. math::
  :label: eq98

  \tilde{v}^{\alpha \beta}_{k=0} = \frac{4 \pi}{\Omega}
  \int_0^{r_\text{end}} dr\ r^2 v_l^{\alpha \beta}(r)\:,

where :math:`r_{\text{end}}` is some cutoff value after which the
potential tails precisely cancel.

Neutralizing background terms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For systems with a net charge, such as the one-component plasma
(jellium), we add a uniform background charge, which makes the system
neutral.  When we do this, we must add a term that comes from the
interaction of the particle with the neutral background.  It is a
constant term, independent of the particle positions.  In general, we
have a compensating background for each species, which largely cancels
out for neutral systems.

.. math::
  :label: eq99

  V_\text{background} = -\frac{1}{2} \sum_\alpha \left(N^\alpha\right)^2
  v^{\alpha \alpha}_{s\mathbf{0}}
  -\sum_{\alpha > \beta} N_\alpha N_\beta
  v^{\alpha\beta}_{s\mathbf{0}}\:,

where :math:`v^{\alpha \beta}_{s\mathbf{0}}` is given by

.. math::
  :label: eq100

  \begin{aligned}
  v^{\alpha \beta}_{s\mathbf{0}} & = & \frac{1}{\Omega} \int_0^{r_c} d^3 r\
  v^{\alpha \beta}_s(r)\:, \\
  & = & \frac{4 \pi}{\Omega} \int_0^{r_c} r^2 v_s(r) \ dr \nonumber\:.\end{aligned}

Combining terms
~~~~~~~~~~~~~~~

Here, we sum all of the terms we computed in the previous sections:

.. math::
  :label: eq101

  \begin{aligned}
  V & = & \sum_{\alpha > \beta} \left[\sum_{i,j} v_s(|\mathbf{r}_i^\alpha
    -\mathbf{r}_j^\beta|) + \sum_\mathbf{k}\mathcal{R}e\left(\rho_\mathbf{k}^\alpha
    \rho_{-\mathbf{k}}^\beta\right)v^{\alpha\beta}_k  -N^\alpha N^\beta
    v^{\alpha \beta}_{s\mathbf{0}}  \right] \nonumber\:, \\
  & + & \sum_\alpha \left[ N^\alpha v_M^\alpha + \sum_{i>j} v_s(|\mathbf{r}_i^\alpha -
    \mathbf{r}_j^\alpha|) + \frac{1}{2} \sum_\mathbf{k}\left( |\rho_\mathbf{k}^\alpha|^2 -
    N\right) v^{\alpha\alpha}_\mathbf{k}-\frac{1}{2}\left(N_\alpha\right)^2 v_{s\mathbf{0}}^{\alpha\alpha}\right] \nonumber\:, \\
  & = & \sum_{\alpha > \beta} \left[\sum_{i,j} v_s(|\mathbf{r}_i^\alpha
    -\mathbf{r}_j^\beta|) + \sum_\mathbf{k}\mathcal{R}e\left(\rho_\mathbf{k}^\alpha
    \rho_{-\mathbf{k}}^\beta\right) v^{\alpha \beta}_k   -N^\alpha N^\beta
    v^{\alpha \beta}_{s\mathbf{0}}  +\tilde{V}_{k=0} \right]\:, \\
  & + & \sum_\alpha \left[ -\frac{N^\alpha v_l^{\alpha \alpha}(0)}{2}  + \sum_{i>j} v_s(|\mathbf{r}_i^\alpha -
    \mathbf{r}_j^\alpha|) + \frac{1}{2} \sum_\mathbf{k}|\rho_\mathbf{k}^\alpha|^2 v^{\alpha\alpha}_\mathbf{k}- \frac{1}{2}\left(N_\alpha\right)^2
    v_{s\mathbf{0}}^{\alpha\alpha} +\tilde{V}_{k=0}\right]  \nonumber\:.\end{aligned}

Computing the reciprocal potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now we return to :eq:`eq86`. Without loss of generality, we
define for convenience :math:`\mathbf{k}= k\hat{\mathbf{z}}`.

.. math::
  :label: eq102

  v^{\alpha \beta}_k = \frac{2\pi}{\Omega} \int_0^\infty dr \int_{-1}^1
    d\cos(\theta) \ r^2 e^{-i k r \cos(\theta)} v_l^{\alpha \beta}(r)\:.

We do the angular integral first.  By inversion symmetry, the
imaginary part of the integral vanishes, yielding

.. math::
  :label: eq103

  v^{\alpha \beta}_k = \frac{4\pi}{\Omega k}\int _0^\infty dr\ r \sin(kr)
  v^{\alpha \beta}_l(r)\:.

The Coulomb potential
~~~~~~~~~~~~~~~~~~~~~

For the case of the Coulomb potential, the preceding integral is not
formally convergent if we do the integral naively. We may remedy the
situation by including a convergence factor, :math:`e^{-k_0 r}`. For a
potential of the form :math:`v^{\text{coul}}(r) = q_1 q_2/r`, this
yields

.. math::
  :label: eq104

  \begin{aligned}
  v^{\text{screened coul}}_k & = & \frac{4\pi q_1 q_2}{\Omega k} \int_0^\infty dr\ \sin(kr)
  e^{-k_0r}\:, \\
  & = & \frac{4\pi q_1 q_2}{\Omega (k^2 + k_0^2)}\:.\end{aligned}

Allowing the convergence factor to tend to zero, we have

.. math::
  :label: eq105

  v_k^\text{coul} = \frac{4 \pi q_1 q_2}{\Omega k^2}\:.

For more generalized potentials with a Coulomb tail, we cannot
evaluate :eq:`eq103` numerically but must handle the coulomb part
analytically.  In this case, we have

.. math::
  :label: eq106

  v_k^{\alpha \beta} = \frac{4\pi}{\Omega}
  \left\{ \frac{q_1 q_2}{k^2} + \int_0^\infty dr \ r \sin(kr) \left[ v_l^{\alpha \beta}(r) -
    \frac{q_1 q_2}{r} \right] \right\}\:.

Efficient calculation methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fast computation of :math:`\rho_\mathbf{k}`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We wish to quickly calculate the quantity

.. math::
  :label: eq107

  \rho_\mathbf{k}^\alpha \equiv \sum_i e^{i\mathbf{k}\cdot r_i^\alpha}\:.

First, we write

.. math::
  :label: eq108

  \begin{aligned}
  \mathbf{k}& = & m_1 \mathbf{b}_1 + m_2 \mathbf{b}_2 + m_3 \mathbf{b}_3\:, \\
  \mathbf{k}\cdot \mathbf{r}_i^\alpha & = &  m_1 \mathbf{b}_1 \cdot \mathbf{r}_i^\alpha +
  m_2 \mathbf{b}_2 \cdot \mathbf{r}_i^\alpha + m_3 \mathbf{b}_3 \cdot \mathbf{r}_i^\alpha\:, \\
  e^{i\mathbf{k}\cdot r_i^\alpha} & = &
  {\underbrace{\left[e^{i \mathbf{b}_1 \cdot\mathbf{r}_i^\alpha}\right]}_{C^{i\alpha}_1}}^{m_1}
  {\underbrace{\left[e^{i \mathbf{b}_2 \cdot\mathbf{r}_i^\alpha}\right]}_{C^{i\alpha}_2}}^{m_2}
  {\underbrace{\left[e^{i \mathbf{b}_3 \cdot\mathbf{r}_i^\alpha}\right]}_{C^{i\alpha}_3}}^{m_3}\:.\end{aligned}

Now, we note that

.. math::
  :label: eq109

  ^{m_1} = C^{i\alpha}_1 [C^{i\alpha}]^{(m_1-1)}\:.

This allows us to recursively build up an array of the
:math:`C^{i\alpha}`\ s and then compute :math:`\rho_\mathbf{k}` for all
:math:`\mathbf{k}`-vectors by looping over all k-vectors, requiring only
two complex multiplies per particle per :math:`\mathbf{k}`.

.. centered:: Algorithm to quickly calculate :math:`\rho_\mathbf{k}^\alpha`.

|   Create list of :math:`\mathbf{k}`-vectors and corresponding :math:`(m_1, m_2, m_3)` indices.
|   **for all** :math:`\alpha \in` species
|     Zero out :math:`\rho_\mathbf{k}^\alpha`
|     **for all** :math:`i \in` particles **do**
|       **for** :math:`j \in [1\cdots3]` **do**
|         Compute :math:`C^{i \alpha}_j \equiv e^{i \mathbf{b}_j \cdot  \mathbf{r}^{\alpha}_i}`
|         **for** :math:`m \in [-m_{\text{max}}\dots m_{\text{max}}]` **do**
|           Compute :math:`[C^{i \alpha}_j]^m` and store in array
|         **end for**
|       **end for**
|       **for all** :math:`(m_1, m_2, m_3) \in` index list **do**
|         Compute :math:`e^{i \mathbf{k}\cdot r^\alpha_i} = [C^{i\alpha}_1]^{m_1} [C^{i\alpha}_2]^{m_2}[C^{i\alpha}_3]^{m_3}` from array
|       **end for**
|     **end for**
|   **end for**

Gaussian charge screening breakup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This original approach to the short- and long-range breakup adds an
opposite screening charge of Gaussian shape around each point charge.
It then removes the charge in the long-range part of the potential.
In this potential,

.. math::
  :label: eq110

  v_{\text{long}}(r) = \frac{q_1 q_2}{r} \text{erf}(\alpha r)\:,

where :math:`\alpha` is an adjustable parameter used to control how
short ranged the potential should be. If the box size is :math:`L`, a
typical value for :math:`\alpha` might be :math:`7/(Lq_1 q_2)`. We
should note that this form for the long-range potential should also work
for any general potential with a Coulomb tail (e.g., pseudo-Hamiltonian
potentials. For this form of the long-range potential, we have in
:math:`k`-space

.. math::
  :label: eq111

  v_k = \frac{4\pi q_1 q_2 \exp\left[\frac{-k^2}{4\alpha^2}\right]}{\Omega k^2}\:.

Optimized breakup method
~~~~~~~~~~~~~~~~~~~~~~~~

In this section, we undertake the task of choosing a
long-range/short-range partitioning of the potential, which is optimal
in that it minimizes the error for given real and :math:`k`-space
cutoffs :math:`r_c` and :math:`k_c`. Here, we slightly modify the method
introduced by Natoli and Ceperley :cite:`Natoli1995`. We
choose :math:`r_c = \frac{1}{2}\min\{L_i\}` so that we require the nearest image in
real-space summation. :math:`k_c` is then chosen to satisfy our accuracy
requirements.

Here we modify our notation slightly to accommodate details not
previously required. We restrict our discussion to the interaction of
two particle species (which may be the same), and drop our species
indices. Thus, we are looking for short- and long-range potentials
defined by

.. math::
  :label: eq112

  v(r) = v^s(r) + v^\ell(r)\:.

Define :math:`v^s_k` and :math:`v^\ell_k` to be the respective Fourier
transforms of the previous equation. The goal is to choose
:math:`v_s(r)` such that its value and first two derivatives vanish at
:math:`r_c`, while making :math:`v^\ell(r)` as smooth as possible so
that :math:`k`-space components, :math:`v^\ell_k`, are very small for
:math:`k>k_c`. Here, we describe how to do this in an optimal way.

Define the periodic potential, :math:`V_p`, as

.. math::
  :label: eq113

  V_p(\mathbf{r}) = \sum_l v(|\mathbf{r}+ \mathbf{l}|),

where :math:`\mathbf{r}` is the displacement between the two particles
and :math:`\mathbf{l}` is a lattice vector. Let us then define our
approximation to this potential, :math:`V_a`, as

.. math::
  :label: eq114

  V_a(\mathbf{r}) = v^s(r) + \sum_{|\mathbf{k}| < k_c} v^\ell_k e^{i \mathbf{k} \cdot \mathbf{r}}\:.

Now, we seek to minimize the RMS error over the cell,

.. math::
  :label: eq115

  \chi^2 = \frac{1}{\Omega}\int_\Omega d^3 \mathbf{r} \
  \left| V_p(\mathbf{r}) - V_a(\mathbf{r})\right|^2\:.

We may write

.. math::
  :label: eq116

  V_p(\mathbf{r}) = \sum_{\mathbf{k}} v_k e^{i \mathbf{k}\cdot \mathbf{r}}\:,

where

.. math::
  :label: eq117

  v_k = \frac{1}{\Omega} \int d^3\mathbf{r}\ e^{-i\mathbf{k}\cdot\mathbf{r}}v(r)\:.

We now need a basis in which to represent the broken-up potential. We
may choose to represent either :math:`v^s(r)` or :math:`v^\ell(r)` in a
real-space basis. Natoli and Ceperley chose the former in their paper.
We choose the latter for a number of reasons. First, singular potentials
are difficult to represent in a linear basis unless the singularity is
explicitly included. This requires a separate basis for each type of
singularity. The short-range potential may have an arbitrary number of
features for :math:`r<r_c` and still be a valid potential. By
construction, however, we desire that :math:`v^\ell(r)` be smooth in
real-space so that its Fourier transform falls off quickly with
increasing :math:`k`. We therefore expect that, in general,
:math:`v^\ell(r)` should be well represented by fewer basis functions
than :math:`v^s(r)`. Therefore, we define

.. math::
  :label: eq118

  v^\ell(r) \equiv
  \begin{cases}
   \sum_{n=0}^{J-1} t_n h_n(r) & \text{for } r \le r_c \\
   v(r) & \text{for } r > r_c.
  \end{cases}\:,

where the :math:`h_n(r)` are a set of :math:`J` basis functions. We
require that the two cases agree on the value and first two derivatives
at :math:`r_c`. We may then define

.. math::
  :label: eq119

  c_{nk} \equiv \frac{1}{\Omega} \int_0^{r_c} d^3 \mathbf{r}\ e^{-i\mathbf{k}\cdot\mathbf{r}} h_n(r)\:.

Similarly, we define

.. math::
  :label: eq120

  x_k \equiv -\frac{1}{\Omega} \int_{r_c}^\infty d^3\mathbf{r}\ e^{-i\mathbf{k}\cdot\mathbf{r}} v(r)\:.

Therefore,

.. math::
  :label: eq121

  v^\ell_k = -x_k + \sum_{n=0}^{J-1} t_n c_{nk}\:.

Because :math:`v^s(r)` goes identically to zero at the box edge, inside
the cell we may write

.. math::
  :label: eq122

  v^s(\mathbf{r}) = \sum_\mathbf{k}v^s_k e^{i\mathbf{k}\cdot \mathbf{r}}\:.

We then write

.. math::
  :label: eq123

  \chi^2 = \frac{1}{\Omega} \int_\Omega d^3 \mathbf{r}\
  \left| \sum_\mathbf{k}e^{i\mathbf{k}\cdot \mathbf{r}} \left(v_k - v^s_k \right)
  -\sum_{|\mathbf{k}| \le k_c} v^\ell_k \right|^2\:.

We see that if we define

.. math::
  :label: eq124

  v^s(r) \equiv v(r) - v^\ell(r)\:.

Then

.. math::
  :label: eq125

  v^\ell_k + v^s_k = v_k\:,

which then cancels out all terms for :math:`|\mathbf{k}| < k_c`. Then we
have

.. math::
  :label: eq126

  \begin{aligned}
  \chi^2 & = & \frac{1}{\Omega} \int_\Omega d^3 \mathbf{r}\
  \left|\sum_{|\mathbf{k}|>k_c} e^{i\mathbf{k}\cdot\mathbf{r}}
  \left(v_k -v^s_k \right)\right|^2\:, \\
  & = & \frac{1}{\Omega} \int_\Omega d^3 \mathbf{r}\
  \left|\sum_{|\mathbf{k}|>k_c} e^{i\mathbf{k}\cdot\mathbf{r}} v^\ell_k \right|^2\:, \\
  & = &
  \frac{1}{\Omega} \int_\Omega d^3 \mathbf{r}
  \left|\sum_{|\mathbf{k}|>k_c} e^{i\mathbf{k}\cdot\mathbf{r}}\left( -x_k + \sum_{n=0}^{J-1} t_n
  c_{nk}\right) \right|^2\:.\end{aligned}

We expand the summation,

.. math::
  :label: eq127

  \chi^2 = \frac{1}{\Omega} \int_\Omega d^3 \mathbf{r}\negthickspace\negthickspace\negthickspace
  \sum_{\{|\mathbf{k}|,|\mathbf{k}'|\}>k_c} \negthickspace\negthickspace\negthickspace\negthickspace\negthickspace
   e^{i(\mathbf{k}-\mathbf{k}')\cdot \mathbf{r}}
  \left(x_k -\sum_{n=0}^{J-1} t_n c_{nk} \right)
  \left(x_k -\sum_{m=0}^{J-1} t_{m} c_{mk'} \right)\:.

We take the derivative w.r.t. :math:`t_{m}`:

.. math::
  :label: eq128

  \frac{\partial (\chi^2)}{\partial t_{m}} =
  \frac{2}{\Omega}\int_\Omega d^3 \mathbf{r}\negthickspace\negthickspace\negthickspace
  \sum_{\{|\mathbf{k}|,|\mathbf{k}'|\}>k_c} \negthickspace\negthickspace\negthickspace\negthickspace\negthickspace
   e^{i(\mathbf{k}-\mathbf{k}')\cdot \mathbf{r}}
  \left(x_k -\sum_{n=0}^{J-1} t_n c_{nk} \right) c_{mk'}\:.

We integrate w.r.t. :math:`\mathbf{r}`, yielding a Kronecker
:math:`\delta`.

.. math::
  :label: eq129

  \frac{\partial (\chi^2)}{\partial t_{m}} =
  2 \negthickspace\negthickspace\negthickspace\negthickspace\negthickspace\negthickspace\negthickspace
  \sum_{\ \ \ \ \{|\mathbf{k}|,|\mathbf{k}'|\}>k_c} \negthickspace\negthickspace\negthickspace\negthickspace\negthickspace\negthickspace\negthickspace\delta_{\mathbf{k}, \mathbf{k}'}
  \left(x_k -\sum_{n=0}^{J-1} t_n c_{nk} \right) c_{mk'}\:.

Summing over :math:`\mathbf{k}'` and equating the derivative to zero, we
find the minimum of our error function is given by

.. math::
  :label: eq130

  \sum_{n=0}^{J-1} \sum_{|\mathbf{k}|>k_c} c_{mk}c_{nk} t_n =
  \sum_{|\mathbf{k}|>k_c} x_k c_{mk}\:,

which is equivalent in form to Equation 19 in
:cite:`Natoli1995`, where we have :math:`x_k` instead of
:math:`V_k`. Thus, we see that we can optimize the short- or long-range
potential simply by choosing to use :math:`V_k` or :math:`x_k` in the
preceding equation. We now define

.. math::
  :label: eq131

  \begin{aligned}
  A_{mn} & \equiv & \sum_{|\mathbf{k}|>k_c} c_{mk} c_{nk}\:, \\
  b_{m} & \equiv & \sum_{|\mathbf{k}|>k_c} x_k c_{mk}\:.\end{aligned}

Thus, it becomes clear that our minimization equations can be cast in
the canonical linear form

.. math::
  :label: eq132

  \mathbf{A}\mathbf{t} = \mathbf{b}\:.

Solution by SVD
^^^^^^^^^^^^^^^

In practice, we note that the matrix :math:`\mathbf{A}` frequently
becomes singular in practice. For this reason, we use the singular value
decomposition to solve for :math:`t_n`. This factorization decomposes
:math:`A` as

.. math::
  :label: eq133

  \mathbf{A}= \mathbf{U}\mathbf{S}\mathbf{V}^T\:,

where :math:`\mathbf{U}^T\mathbf{U}= \mathbf{V}^T\mathbf{V}= 1` and
:math:`\mathbf{S}` is diagonal. In this form, we have

.. math::
  :label: eq134

  \mathbf{t} = \sum_{i=0}^{J-1} \left( \frac{\mathbf{U}_{(i)} \cdot
    \mathbf{b}}{\mathbf{S}_{ii}} \right) \mathbf{V}_{(i)}\:,

where the parenthesized subscripts refer to columns. The advantage of
this form is that if :math:`\mathbf{S}_{ii}` is zero or very near zero,
the contribution of the :math:`i^{\text{th}}` of :math:`\mathbf{V}` may
be neglected since it represents a numerical instability and has little
physical meaning. It represents the fact that the system cannot
distinguish between two linear combinations of the basis functions.
Using the SVD in this manner is guaranteed to be stable. This
decomposition is available in LAPACK in the DGESVD subroutine.

.. _constraints:

Constraining Values
^^^^^^^^^^^^^^^^^^^

Often, we wish to constrain the value of :math:`t_n` to have a fixed
value to enforce a boundary condition, for example. To do this, we
define

.. math::
  :label: eq135

  \mathbf{b}' \equiv \mathbf{b}- t_n \mathbf{A}_{(n)}\:.

We then define :math:`\mathbf{A}^*` as :math:`\mathbf{A}` with the
:math:`n^{\text{th}}` row and column removed and :math:`\mathbf{b}^*` as
:math:`\mathbf{b}'` with the :math:`n^{\text{th}}` element removed. Then
we solve the reduced equation
:math:`\mathbf{A}^* \mathbf{t}^* = \mathbf{b}^*` and finally insert
:math:`t_n` back into the appropriate place in :math:`\mathbf{t}^*` to
recover the complete, constrained vector :math:`\mathbf{t}`. This may be
trivially generalized to an arbitrary number of constraints.

The LPQHI basis
^^^^^^^^^^^^^^^

The preceding discussion is general and independent of the basis used to
represent :math:`v^\ell(r)`. In this section, we introduce a convenient
basis of localized interpolant functions, similar to those used for
splines, which have a number of properties that are convenient for our
purposes.

First, we divide the region from 0 to :math:`r_c` into :math:`M-1`
subregions, bounded above and below by points we term *knots*, defined
by :math:`r_j
\equiv j\Delta`, where :math:`\Delta \equiv r_c/(M-1)`. We then define
compact basis elements, :math:`h_{j\alpha}`, which span the region
:math:`[r_{j-1},r_{j+1}]`, except for :math:`j=0` and :math:`j=M`. For
:math:`j=0`, only the region :math:`[r_0,r_1]`, while for :math:`j=M`,
only :math:`[r_{M-1}, r_M]`. Thus, the index :math:`j` identifies the
knot the element is centered on, while :math:`\alpha` is an integer from
0 to 2 indicating one of three function shapes. The dual index can be
mapped to the preceding single index by the relation
:math:`n = 3j + \alpha`. The basis functions are then defined as

.. math::
  :label: eq136

  h_{j\alpha}(r) =
  \begin{cases}
  \ \ \ \, \Delta^\alpha \, \, \sum_{n=0}^5 S_{\alpha n}
  \left( \frac{r-r_j}{\Delta}\right)^n,    & r_j < r \le r_{j+1} \\
  (-\Delta)^\alpha \sum_{n=0}^5 S_{\alpha n}
  \left( \frac{r_j-r}{\Delta}\right)^n,    & r_{j-1} < r \le r_j \\
  \quad\quad\quad\quad\quad 0, & \text{otherwise}\:,
  \end{cases}

where the matrix :math:`S_{\alpha n}` is given by

.. math::
  :label: eq137

  S =
  \left[\begin{matrix}
  1 & 0 & 0 & -10 & 15 & -6 \\
  0 & 1 & 0 & -6  &  8 & -3 \\
  0 & 0 & \frac{1}{2} & -\frac{3}{2} & \frac{3}{2} & -\frac{1}{2}
  \end{matrix}\right]\:.

:numref:`fig26` shows plots of these function shapes.

.. _fig26:
.. figure:: /figs/LPQHI.png
  :width: 500
  :align: center

  Basis functions :math:`h_{j0}`, :math:`h_{j1}`, and :math:`h_{j2}` are
  shown. We note that at the left and right extremes, the values and first
  two derivatives of the functions are zero; while at the center,
  :math:`h_{j0}` has a value of 1, :math:`h_{j1}` has a first derivative
  of 1, and :math:`h_{j2}` has a second derivative of 1.

The basis functions have the property that at the left and right
extremes (i.e., :math:`r_{j-1}` and :math:`r_{j+1}`) their values and
first two derivatives are zero. At the center, :math:`r_j`, we have the
properties

.. math::
  :label: eq138

  \begin{aligned}
  h_{j0}(r_j)=1, & h'_{j0}(r_j)=0, & h''_{j0}(r_j)= 0\:, \\
  h_{j1}(r_j)=0, & h'_{j1}(r_j)=1, & h''_{j1}(r_j)= 0\:, \\
  h_{j2}(r_j)=0, & h'_{j2}(r_j)=0, & h''_{j2}(r_j)= 1\:. \end{aligned}

These properties allow the control of the value and first two derivatives
of the represented function at any knot value simply by setting the
coefficients of the basis functions centered around that knot.  Used
in combination with the method described in
:ref:`constraints`, boundary conditions can easily be
enforced.  In our case, we wish require that

.. math::
  :label: eq139

  h_{M0} = v(r_c), \ \ h_{M1} = v'(r_c), \ \ \text{and} \ \  h_{M2} = v''(r_c)\:.

This ensures that :math:`v^s` and its first two derivatives vanish at
:math:`r_c`.

Fourier coefficients
^^^^^^^^^^^^^^^^^^^^

We wish now to calculate the Fourier transforms of the basis
functions, defined as

.. math::
  :label: eq140

  c_{j\alpha k} \equiv \frac{1}{\Omega} \int_0^{r_c} d^3 \mathbf{r}
  e^{-i \mathbf{k}\cdot \mathbf{r}} h_{j\alpha}(r)\:.

We may then write,

.. math::
  :label: eq141

  c_{j\alpha k} =
  \begin{cases}
  \Delta^\alpha \sum_{n=0}^5 S_{\alpha n} D^+_{0 k n}, & j = 0 \\
  \Delta^\alpha \sum_{n=0}^5 S_{\alpha n} (-1)^{\alpha+n} D^-_{M k n}, &
  j = M \\
  \Delta^\alpha \sum_{n=0}^5 S_{\alpha n}
  \left[ D^+_{j k n} + (-1)^{\alpha+n}D^-_{j k n} \right] & \text{otherwise}\:,
  \end{cases}

where

.. math::
  :label: eq142

  D^{\pm}_{jkn} \equiv \frac{1}{\Omega} \int_{r_j}^{r_{j\pm1}} d^3\!\mathbf{r}\
  e^{-i\mathbf{k}\cdot \mathbf{r}} \left( \frac{r-r_j}{\Delta}\right)^n\:.

We then further make the definition that

.. math::
  :label: eq143

  D^{\pm}_{jkn} = \pm \frac{4\pi}{k \Omega}
  \left[ \Delta \text{Im}\left(E^{\pm}_{jk(n+1)}\right) +
  r_j \text{Im}\left(E^{\pm}_{jkn}\right)\right]\:.

It can then be shown that

.. math::
  :label: eq144

  E^{\pm}_{jkn} =
  \begin{cases}
  -\frac{i}{k} e^{ikr_j} \left( e^{\pm i k \Delta} - 1 \right) &
  \text{if } n=0, \\
  -\frac{i}{k}
  \left[ \left(\pm1\right)^n e^{i k (r_j \pm \Delta)} - \frac{n}{\Delta}
  E^\pm_{jk(n-1)}  \right] & \text{otherwise}\:.
  \end{cases}

Note that these equations correct typographical errors present in :cite:`Natoli1995`.

Enumerating :math:`k`-points
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We note that the summations over :math:`k`, which are ubiquitous in this
paper, require enumeration of the :math:`k`-vectors. In particular, we
should sum over all :math:`|\mathbf{k}| > k_c`. In practice, we must
limit our summation to some finite cutoff value
:math:`k_c < |\mathbf{k}| < k_{\text{max}}`, where
:math:`k_{\text{max}}` should be on the order of :math:`3,000/L`, where
:math:`L` is the minimum box dimension. Enumerating these vectors in a
naive fashion even for this finite cutoff would prove quite prohibitive,
as it would require :math:`\sim10^9` vectors.

Our first optimization comes in realizing that all quantities in this
calculation require only :math:`|\mathbf{k}|` and not :math:`\mathbf{k}`
itself. Thus, we may take advantage of the great degeneracy of
:math:`|\mathbf{k}|`. We create a list of :math:`(k,N)` pairs, where
:math:`N` is the number of vectors with magnitude :math:`k`. We make
nested loops over :math:`n_1`, :math:`n_2`, and :math:`n_3`, yielding
:math:`\mathbf{k}= n_1 \mathbf{b}_1 + n_2 \mathbf{b}_2 + n_3
\mathbf{b}_3`. If :math:`|\mathbf{k}|` is in the required range, we
check to see whether there is already an entry with that magnitude on
our list and increment the corresponding :math:`N` if there is, or
create a new entry if not. Doing so typically saves a factor of
:math:`\sim200` in storage and computation.

This reduction is not sufficient for large :math:`k_max` since it
requires that we still look over :math:`10^9` entries. To further reduce
costs, we may pick an intermediate cutoff, :math:`k_{\text{cont}}`,
above which we will approximate the degeneracy assuming a continuum of
:math:`k`-points. We stop our exact enumeration at
:math:`k_{\text{cont}}` and then add :math:`\sim1,000` points,
:math:`k_i`, uniformly spaced between :math:`k_{\text{cont}}` and
:math:`k_{\text{max}}`. We then approximate the degeneracy by

.. math::
  :label: eq145

  N_i = \frac{4 \pi}{3} \frac{\left( k_b^3 -k_a^3\right)}{(2\pi)^3/\Omega}\:,

where :math:`k_b = (k_i + k_{i+1})/2` and :math:`k_a = (k_i + k_{i-1})`.
In doing so, we typically reduce our total number of k-points to sum
more than :math:`\sim2,500` from the :math:`10^9` we had to start.

Calculating :math:`x_k`'s
^^^^^^^^^^^^^^^^^^^^^^^^^

The Coulomb potential
.....................

For :math:`v(r) = \frac{1}{r}`, :math:`x_k` is given by

.. math::
  :label: eq146

  x_k^{\text{coulomb}} = -\frac{4 \pi}{\Omega k^2} \cos(k r_c)\:.

The :math:`1/r^2` potential
...........................

For :math:`v(r) = \frac{1}{r^2}`, :math:`x_k` is given by

.. math::
  :label: eq147

  x_k^{1/r^2} = \frac{4 \pi}{\omega k}
  \left[ \text{Si}(k r_c) -\frac{\pi}{2}\right],

where the *sin integral*, :math:`\text{Si}(z)`, is given by

.. math::
  :label: eq148

  \text{Si}(z) \equiv \int_0^z \frac{\sin \ t}{t} dt\:.

The :math:`1/r^3` potential
...........................

For :math:`v(r) = \frac{1}{r^3}`, :math:`x_k` is given by

.. math::
  :label: eq149

  x_k^{1/r^2} = \frac{4 \pi}{\omega k}
  \left[ \text{Si}(k r_c) -\frac{\pi}{2}\right],

where the *cosine integral*, :math:`\text{Ci}(z)`, is given by

.. math::
  :label: eq150

  \text{Ci}(z) \equiv -\int_z^\infty \frac{\cos t}{t} dt\:.

The :math:`1/r^4` potential
...........................

For :math:`v(r) = \frac{1}{r^4}`, :math:`x_k` is given by

.. math::
  :label: eq151

  x_k^{1/r^4} = -\frac{4 \pi}{\Omega k}
  \left\{
  \frac{k \cos(k r_c)}{2 r_c} + \frac{\sin(k r_c)}{2r_c^2} + \frac{k^2}{2} \left[ \text{Si}(k r_c) - \frac{\pi}{2}\right]\right\}\:.

Feature: Optimized long-range breakup (Ewald) 2
-----------------------------------------------

Given a lattice of vectors :math:`\mathbf{L}`, its associated reciprocal
lattice of vectors :math:`\mathbf{k}` and a function
:math:`\psi(\mathbf{r})` periodic on the lattice we define its Fourier
transform :math:`\widetilde{\psi}(\mathbf{k})` as

.. math::
  :label: eq152

  \widetilde{\psi}(\mathbf{k})=\frac{1}{\Omega}\int_\Omega d\mathbf{r}\psi(\mathbf{r}) e^{-i\mathbf{k}\mathbf{r}}\:,

where we indicated both the cell domain and the cell volume by
:math:`\Omega`. :math:`\psi(\mathbf{r})` can then be expressed as

.. math::
  :label: eq153

  \psi(\mathbf{r})=\sum_{\mathbf{k}} \widetilde{\psi}(\mathbf{k})e^{i\mathbf{k}\mathbf{r}}\:.

The potential generated by charges sitting on the lattice positions at a
particular point :math:`\mathbf{r}` inside the cell is given by

.. math::
  :label: eq154

  V(\mathbf{r})=\sum_{\mathbf{L}}v(|\mathbf{r}+\mathbf{L}|)\:,

and its Fourier transform can be explicitly written as a function of
:math:`V` or :math:`v`

.. math::
  :label: eq155

  \widetilde{V}(\mathbf{k})=\frac{1}{\Omega}\int_\Omega d\mathbf{r}V(\mathbf{r}) e^{-i\mathbf{k}\mathbf{r}}=
  \frac{1}{\Omega}\int_{\mathbb{R}^3} d\mathbf{r}v(\mathbf{r}) e^{-i\mathbf{k}\mathbf{r}}\:,

where :math:`\mathbb{R}^3` denotes the whole 3D space. We now want to
find the best (“best” to be defined later) approximate potential of the
form

.. math::
  :label: eq156

  V_a(\mathbf{r})=\sum_{k\le k_c} \widetilde{Y}(k) e^{i\mathbf{k}\mathbf{r}} + W(r)\:,

where :math:`W(r)` has been chosen to go smoothly to :math:`0` when
:math:`r=r_c`, being :math:`r_c` lower or equal to the Wigner-Seitz
radius of the cell. Note also the cutoff :math:`k_c` on the momentum
summation.

The best form of :math:`\widetilde{Y}(k)` and :math:`W(r)` is given by
minimizing

.. math::
  :label: eq157

  \chi^2=\frac{1}{\Omega}\int d\mathbf{r}\left(V(\mathbf{r})-W(\mathbf{r})-
    \sum_{k\le k_c}\widetilde{Y}(k)e^{i\mathbf{k}\mathbf{r}}\right)^2
    \:,

or the reciprocal space equivalent

.. math::
  :label: eq158

  \chi^2=\sum_{k\le k_c}(\widetilde{V}(k)-\widetilde{W}(k)-\widetilde{Y}(k))^2+\sum_{k>k_c}(\widetilde{V}(k)-\widetilde{W}(k))^2
    \:.

:eq:`eq158` follows from :eq:`eq157` and the unitarity
(norm conservation) of the Fourier transform.

This last condition is minimized by

.. math::
  :label: eq159

  \widetilde{Y}(k)=\widetilde{V}(k)-\widetilde{W}(k)\qquad \min_{\widetilde{W}(k)}\sum_{k>k_c}(\widetilde{V}(k)-\widetilde{W}(k))^2.

We now use a set of basis function :math:`c_i(r)` vanishing smoothly at
:math:`r_c` to expand :math:`W(r)`; that is,

.. math::
  :label: eq160

  W(r)=\sum_i t_i c_i(r)\qquad\text{or}\qquad \widetilde{W}(k)=\sum_i t_i \widetilde{c}_i(k)\:.

Inserting the reciprocal space expansion of :math:`\widetilde{W}` in the
second condition of :eq:`eq159` and minimizing with respect to
:math:`t_i` leads immediately to the linear system
:math:`\mathbf{A}\mathbf{t}=\mathbf{b}` where

.. math::
  :label: eq161

   \begin{aligned}
   A_{ij}=\sum_{k>k_c}\widetilde{c}_i(k)\widetilde{c}_j(k)\qquad b_j=\sum_{k>k_c} V(k) \widetilde{c}_j(k)
   \:.\end{aligned}

Basis functions
~~~~~~~~~~~~~~~

The basis functions are splines. We define a uniform grid with
:math:`N_{\text{knot}}` uniformly spaced knots at position
:math:`r_i=i\frac{r_c}{N_{\text{knot}}}`, where
:math:`i\in[0,N_{\text{knot}}-1]`. On each knot we center :math:`m+1`
piecewise polynomials :math:`c_{i\alpha}(r)` with
:math:`\alpha\in[0,m]`, defined as

.. math::
  :label: eq162

   \begin{aligned}
   c_{i\alpha}(r)=\begin{cases}
   \Delta^\alpha \sum_{n=0}^\mathcal{N} S_{\alpha n}(\frac{r-r_i}{\Delta})^n & r_i<r\le r_{i+1} \\
   \Delta^{-\alpha} \sum_{n=0}^\mathcal{N} S_{\alpha n}(\frac{r_i-r}{\Delta})^n & r_{i-1}<r\le r_i \\
   0 & |r-r_i| > \Delta
   \end{cases}
   \:.\end{aligned}

These functions and their derivatives are, by construction, continuous
and odd (even) (with respect to :math:`r-r_i\rightarrow r_i-r`) when
:math:`\alpha` is odd (even). We further ask them to satisfy

.. math::
  :label: eq163

   \begin{aligned}
   \left.\frac{d^\beta}{dr^\beta} c_{i\alpha}(r)\right|_{r=r_i}=
   \delta_{\alpha\beta} \quad \beta\in[0,m]\:,\\
   \left.\frac{d^{\beta}}{dr^{\beta}} c_{i\alpha}(r)\right|_{r=r_{i+1}}=0\quad \beta\in[0,m]
   \:.\end{aligned}

(The parity of the functions guarantees that the second constraint is
satisfied at :math:`r_{i-1}` as well). These constraints have a simple
interpretation: the basis functions and their first :math:`m`
derivatives are :math:`0` on the boundary of the subinterval where they
are defined; the only function to have a nonzero :math:`\beta`-th
derivative in :math:`r_i` is :math:`c_{i\beta}`. These :math:`2(m+1)`
constraints therefore impose :math:`\mathcal{N}=2m+1`. Inserting the
definitions of :eq:`eq162` in the constraints of
:eq:`eq163` leads to the set of :math:`2(m+1)` linear equation
that fixes the value of :math:`S_{\alpha n}`:

.. math::
  :label: eq164

   \begin{aligned}
   \Delta^{\alpha-\beta} S_{\alpha\beta} \beta!=\delta_{\alpha\beta}
   \\
   \Delta^{\alpha-\beta}\sum_{n=\beta}^{2m+1} S_{\alpha n} \frac{n!}{(n-\beta)!}=0\:.\end{aligned}

We can further simplify inserting the first of these equations into the
second and write the linear system as

.. math::
  :label: eq165

   \sum_{n=m+1}^{2m+1} S_{\alpha n} \frac{n!}{(n-\beta)!}=\begin{cases}
   -\frac{1}{(\alpha-\beta)!}& \alpha\ge \beta \\
   0 & \alpha < \beta
   \end{cases}
   \:.

Fourier components of the basis functions in 3D
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`k\ne 0`, non-Coulomb case
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We now need to evaluate the Fourier transform
:math:`\widetilde{c}_{i\alpha}(k)`. Let us start by writing the
definition

.. math::
  :label: eq166

  \widetilde{c}_{i\alpha}(k)=\frac{1}{\omega}\int_\Omega d\mathbf{r}e^{-i\mathbf{k}\mathbf{r}} c_{i\alpha}(r)\:.

Because :math:`c_{i\alpha}` is different from zero only inside the
spherical crown defined by :math:`r_{i-1}<r<r_i`, we can conveniently
compute the integral in spherical coordinates as

.. math::
  :label: eq167

   \begin{aligned}
   \widetilde{c}_{i\alpha}(k)=\Delta^\alpha\sum_{n=0}^\mathcal{N} S_{\alpha n} \left[
   D_{in}^+(k) +w_{\text{knot}}(-1)^{\alpha+n}D_{in}^-(k)\right]\:,
   \end{aligned}

where we used the definition :math:`w_{\text{knot}}=1-\delta_{i0}` and

.. math::
  :label: eq168

   D_{in}^\pm(k)=\pm\frac{4\pi}{k\Omega}\text{Im}\left[\int_{r_i}^{r_i\pm\Delta}
   dr\left(\frac{r-r_i}{\Delta}\right)^n r e^{ikr}\right]\:,

obtained by integrating the angular part of the Fourier transform. Using
the identity

.. math::
  :label: eq169

  \left(\frac{r-r_i}{\Delta}\right)^n r=\Delta\left(\frac{r-r_i}{\Delta}\right)^{n+1}+\left(\frac{r-r_i}{\Delta}\right)^n r_i

and the definition

.. math::
  :label: eq170

   E_{in}^\pm(k)=\int_{r_i}^{r_i\pm\Delta}
   dr\left(\frac{r-r_i}{\Delta}\right)^n e^{ikr}\:,

we rewrite Equation :eq:`eq168` as

.. math::

   \begin{aligned}
   D_{in}^\pm(k)=\pm\frac{4\pi}{k\Omega}\text{Im}\left[\Delta E_{i(n+1)}^\pm(k)+
   r_i E_{in}^\pm(k)\right]\:.
   \end{aligned}

Finally, using integration by part, we can define :math:`E^\pm_{in}`
recursively as

.. math::
  :label: eq171

   \begin{aligned}
   E^\pm_{in}(k)=\frac{1}{ik}\left[(\pm)^ne^{ik(r_i\pm\Delta)}-\frac{n}{\Delta}
   E^\pm_{i(n-1)}(k)\right]\:.
   \end{aligned}

Starting from the :math:`n=0` term,

.. math::
  :label: eq172

   \begin{aligned}
   E^\pm_{i0}(k)=\frac{1}{ik}e^{ikr_i}\left(e^{\pm ik\Delta}-1\right)\:.
   \end{aligned}

:math:`k\ne 0`, Coulomb case
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To efficiently treat the Coulomb divergence at the origin, it is
convenient to use a basis set :math:`c_{i\alpha}^{\text{coul}}` of the
form

.. math::
  :label: eq173

  c_{i\alpha}^{\text{coul}}=\frac{c_{i\alpha}}{r}\:.

An equation identical to :eq:`eq168` holds but with the modified
definition

.. math::
  :label: eq174

   D_{in}^\pm(k)=\pm\frac{4\pi}{k\Omega}\text{Im}\left[\int_{r_i}^{r_i\pm\Delta}
   dr\left(\frac{r-r_i}{\Delta}\right)^n e^{ikr}\right]\:,

which can be simply expressed using :math:`E^\pm_{in}(k)` as

.. math::
  :label: eq175

   \begin{aligned}
   D_{in}^\pm(k)=\pm\frac{4\pi}{k\Omega}\text{Im}\left[E_{in}^\pm(k)\right]\:.
   \end{aligned}

:math:`k=0` Coulomb and non-Coulomb case
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The definitions of :math:`D_{in}(k)` given so far are clearly
incompatible with the choice :math:`k=0` (they involve division by
:math:`k`). For the non-Coulomb case, the starting definition is

.. math::
  :label: eq176

   D^\pm_{in}(0)=\pm\frac{4\pi}{\Omega}\int_{r_i}^{r_i\pm\Delta}r^2
   \left(\frac{r-r_i}{\Delta}\right)^ndr\:.

Using the definition :math:`I_n^\pm=(\pm)^{n+1}\Delta/(n+1)`, we can
express this as

.. math::
  :label: eq177

   \begin{aligned}
   D^\pm_{in}(0)=\pm\frac{4\pi}{\Omega}\left[\Delta^2 I_{n+2}^\pm
   +2r_i\Delta I_{n+1}^\pm+2r_i^2I_n^\pm\right]\:.
   \end{aligned}

For the Coulomb case, we get

.. math::
  :label: eq178

   \begin{aligned}
   D^\pm_{in}(0)=\pm\frac{4\pi}{\Omega}\left(
   \Delta I^\pm_{n+1} + r_i I^\pm_n\right)\:.
   \end{aligned}

Fourier components of the basis functions in 2D
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:eq:`eq167` still holds provided we define

.. math::
  :label: eq179

   D^\pm_{in}(k)=\pm\frac{2\pi}{\Omega \Delta^n} \sum_{j=0}^n \binom{n}{j}
   (-r_i)^{n-j}\int_{r_i}^{r_i\pm \Delta}\negthickspace \negthickspace
   \negthickspace \negthickspace \negthickspace \negthickspace \negthickspace
   dr r^{j+1-C} J_0(kr)\:,


where :math:`C=1(=0)` for the Coulomb(non-Coulomb) case.
:eq:`eq179` is obtained using the integral definition of the
zero order Bessel function of the first kind:

.. math::
  :label: eq180

  J_0(z)=\frac{1}{\pi}\int_0^\pi e^{iz\cos\theta}d\theta\:,

and the binomial expansion for :math:`(r-r_i)^n`. The integrals can be
computed recursively using the following identities:

.. math::
  :label: eq181

   \begin{aligned}
   &\int dz J_0(z)=\frac{z}{2}\left[\pi J_1(z)H_0(z)+J_0(z)(2-\pi H_1(z))\right]
   \:,\\
   &\int dz z J_0(z)= z J_1(z)
   \:,\\
   &\int dz z^n J_0(z)= z^nJ_1(z)+(n-1)x^{n-1}J_0(z)
   -(n-1)^2\int dz z^{n-2} J_0(z)\:.
   \end{aligned}

The bottom equation of :eq:`eq181` is obtained using the second equation in the same set, integration by part, and
the identity :math:`\int J_1(z) dz =-J_0(z)``. In the top equation, :math:`H_0` and :math:`H_1` are Struve functions.

Construction of the matrix elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the previous equations, we can construct the matrix elements in
:eq:`eq161` and proceed solving for :math:`t_i`. It is
sometimes desirable to put some constraints on the value of :math:`t_i`.
For example, when the Coulomb potential is concerned, we might want to
set :math:`t_{0}=1`. If the first :math:`g` variable is constrained by
:math:`t_{m}=\gamma_m` with :math:`m=[1,g]`, we can simply redefine
:eq:`eq161` as

.. math::
  :label: 182

   \begin{split}
   A_{ij}=&\sum_{k>k_c} \widetilde{c}_i(k)\widetilde{c}_j(k)  \quad i,j\notin[1,g]\:, \\
   b_j=&\sum_{k>k_c} \left(\widetilde{V}(k)-\sum_{m=1}^g \gamma_m \widetilde{c}_m(k)\right)\widetilde{c}_j(k)\quad j\notin[1,g]\:.
   \end{split}

Feature: Cubic spline interpolation
-----------------------------------

.. _ Written by Kenneth P .Esler Jr.
 Originally titled ``Cubic Spline Interpolation in 1, 2 and 3 Dimensions''

We present the basic equations and algorithms necessary to
construct and evaluate cubic interpolating splines in one, two, and
three dimensions.  Equations are provided for both natural and
periodic boundary conditions.

One dimension
~~~~~~~~~~~~~

Let us consider the problem in which we have a function :math:`y(x)`
specified at a discrete set of points :math:`x_i`, such that
:math:`y(x_i) = y_i`. We wish to construct a piecewise cubic polynomial
interpolating function, :math:`f(x)`, which satisfies the following
conditions:

-  :math:`f(x_i) = y_i`.

-  :math:`f'(x_i^-) = f'(x_i^+)`.

-  :math:`f''(x_i^-) = f''(x_i+)`.

Hermite interpolants
^^^^^^^^^^^^^^^^^^^^

In our piecewise representation, we wish to store only the values
:math:`y_i` and first derivatives, :math:`y'_i`, of our function at each
point :math:`x_i`, which we call *knots*. Given this data, we wish to
construct the piecewise cubic function to use between :math:`x_i` and
:math:`x_{i+1}`, which satisfies the preceding conditions. In
particular, we wish to find the unique cubic polynomial, :math:`P(x)`,
satisfying

.. math::
  :label: eq183

   \begin{aligned}
   P(x_i)      & = & y_i      \:, \\
   P(x_{i+1})  & = & y_{i+1}  \:, \\
   P'(x_i)     & = & y'_i     \:, \\
   P'(x_{i+1}) & = & y'_{i+1} \:.\end{aligned}

.. math::
  :label: eq184

   \begin{aligned}
   h_i & \equiv & x_{i+1} - {x_i}\:, \\
   t & \equiv & \frac{x-x_i}{h_i}\:.\end{aligned}

We then define the basis functions,

.. math::
  :label: eq185

   \begin{aligned}
   p_1(t) & = & (1+2t)(t-1)^2  \:, \\
   q_1(t) & = & t (t-1)^2\:,       \\
   p_2(t) & = & t^2(3-2t)\:,       \\
   q_2(t) & = & t^2(t-1)\:.       \end{aligned}

On the interval, :math:`(x_i, x_{i+1}]`, we define the interpolating
function

.. math::
  :label: eq186

  P(x) = y_i p_1(t) + y_{i+1}p_2(t) + h\left[y'_i q_1(t) + y'_{i+1} q_2(t)\right]\:.

It can be easily verified that :math:`P(x)` satisfies conditions of
equations 1 through 3 of :eq:`eq183`. It is now left to determine the
proper values for the :math:`y'_i\,`\ s such that the continuity
conditions given previously are satisfied.

By construction, the value of the function and derivative will match at
the knots; that is,

.. math:: P(x_i^-) = P(x_i^+), \ \ \ \ P'(x_i^-) = P'(x_i^+)\:.

Then we must now enforce only the second derivative continuity:

.. math::
  :label: eq187

   \begin{aligned}
   P''(x_i^-) & = & P''(x_i^+)\:,  \\
   \frac{1}{h_{i-1}^2}\left[\rule{0pt}{0.3cm}6 y_{i-1} -6 y_i + h_{i-1}\left(2 y'_{i-1} +4 y'_i\right) \right]& = &
   \frac{1}{h_i^2}\left[\rule{0pt}{0.3cm}-6 y_i + 6 y_{i+1} +h_i\left( -4 y'_i -2 y'_{i+1} \right)\right] \nonumber\:. \end{aligned}

Let us define

.. math::
  :label: eq188

   \begin{aligned}
   \lambda_i & \equiv & \frac{h_i}{2(h_i+h_{i-1})}\:,  \\
   \mu_i & \equiv & \frac{h_{i-1}}{2(h_i+h_{i-1})}  = \frac{1}{2} - \lambda_i\:. \end{aligned}

Then we may rearrange

.. math::
  :label: eq189

   \lambda_i y'_{i-1} + y'_i + \mu_i y'_{i+1} = \underbrace{3 \left[\lambda_i \frac{y_i - y_{i-1}}{h_{i-1}} + \mu_i \frac{y_{i+1}
       - y_i}{h_i} \right] }_{d_i}\:.

This equation holds for all :math:`0<i<(N-1)`, so we have a tridiagonal
set of equations. The equations for :math:`i=0` and :math:`i=N-1` depend
on the boundary conditions we are using.

Periodic boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For periodic boundary conditions, we have

.. math::
  :label: eq190

   \begin{matrix}
   y'_0           & +  & \mu_0 y'_1     &   &                   &            & \dots                   & +  \lambda_0 y'_{N-1} & = & d_0\:,  \\
   \lambda_1 y'_0 & +  & y'_1           & + &  \mu_1 y'_2       &            & \dots                   &                       & = & d_1\:,  \\
                  &    & \lambda_2 y'_1 & + &  y'_2           + & \mu_2 y'_3 & \dots                   &                       & = & d_2\:,  \\
                  &    &                &   &  \vdots           &            &                         &                       &   &     \\
   \mu_{N-1} y'_0 &    &                &   &                   &            & +\lambda_{N-1} y'_{N-1} & +  y'_{N-2}           & = & d_3\:.
   \end{matrix}

Or, in matrix form, we have

.. math::
  :label: eq191

     \begin{pmatrix}
     1         & \mu_0     &    0   &   0           & \dots         &      0        & \lambda_0 \\
     \lambda_1 &  1        & \mu_1  &   0           & \dots         &      0        &     0     \\
     0         & \lambda_2 &   1    & \mu_2         & \dots         &      0        &     0     \\
     \vdots    & \vdots    & \vdots & \vdots        & \ddots        &   \vdots      &  \vdots   \\
     0         &   0       &   0    & \lambda_{N-3} &      1        & \mu_{N-3}     &    0      \\
     0         &   0       &   0    &   0           & \lambda_{N-2} &      1        & \mu_{N-2} \\
     \mu_{N-1} &   0       &   0    &   0           &   0           & \lambda_{N-1} &  1
     \end{pmatrix}
     \begin{pmatrix} y'_0 \\ y'_1 \\ y'_2 \\ \vdots \\ y'_{N-3} \\ y'_{N-2} \\ y'_{N-1} \end{pmatrix} =
     \begin{pmatrix} d_0  \\  d_1 \\  d_2 \\ \vdots \\  d_{N-3} \\  d_{N-2} \\  d_{N-1} \end{pmatrix} .

The system is tridiagonal except for the two elements in the upper right
and lower left corners. These terms complicate the solution a bit,
although it can still be done in :math:`\mathcal{O}(N)` time. We first
proceed down the rows, eliminating the the first non-zero term in each
row by subtracting the appropriate multiple of the previous row. At the
same time, we eliminate the first element in the last row, shifting the
position of the first non-zero element to the right with each iteration.
When we get to the final row, we will have the value for
:math:`y'_{N-1}`. We can then proceed back upward, backsubstituting
values from the rows below to calculate all the derivatives.

Complete boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If we specify the first derivatives of our function at the end points,
we have what is known as *complete* boundary conditions.  The
equations in that case are trivial to solve:

.. math::
  :label: eq192

   \begin{pmatrix}
   1         &  0        &    0   &   0           & \dots         &      0        &     0     \\
   \lambda_1 &  1        & \mu_1  &   0           & \dots         &      0        &     0     \\
   0         & \lambda_2 &   1    & \mu_2         & \dots         &      0        &     0     \\
   \vdots    & \vdots    & \vdots & \vdots        & \ddots        &   \vdots      &  \vdots   \\
   0         &   0       &   0    & \lambda_{N-3} &      1        & \mu_{N-3}     &    0      \\
   0         &   0       &   0    &   0           & \lambda_{N-2} &      1        & \mu_{N-2} \\
   0         &   0       &   0    &   0           &   0           &      0        &  1
   \end{pmatrix}
   \begin{pmatrix} y'_0 \\ y'_1 \\ y'_2 \\ \vdots \\ y'_{N-3} \\ y'_{N-2} \\ y'_{N-1} \end{pmatrix} =
   \begin{pmatrix} d_0  \\  d_1 \\  d_2 \\ \vdots \\  d_{N-3} \\  d_{N-2} \\  d_{N-1} \end{pmatrix} .

This system is completely tridiagonal, and we may solve trivially by
performing row eliminations downward, then proceeding upward as before.

Natural boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

If we do not have information about the derivatives at the boundary
conditions, we may construct a *natural spline*, which assumes the
second derivatives are zero at the end points of our spline. In this
case our system of equations is the following:

.. math::
  :label: eq193

   \begin{pmatrix}
   1         & \frac{1}{2} &    0   &   0           & \dots         &      0        &     0     \\
   \lambda_1 &  1          & \mu_1  &   0           & \dots         &      0        &     0     \\
   0         & \lambda_2   &   1    & \mu_2         & \dots         &      0        &     0     \\
   \vdots    & \vdots      & \vdots & \vdots        & \ddots        &   \vdots      &  \vdots   \\
   0         &   0         &   0    & \lambda_{N-3} &      1        & \mu_{N-3}     &    0      \\
   0         &   0         &   0    &   0           & \lambda_{N-2} &      1        & \mu_{N-2} \\
   0         &   0         &   0    &   0           &   0           &  \frac{1}{2}  &  1
   \end{pmatrix}
   \begin{pmatrix} y'_0 \\ y'_1 \\ y'_2 \\ \vdots \\ y'_{N-3} \\ y'_{N-2} \\ y'_{N-1} \end{pmatrix} =
   \begin{pmatrix} d_0  \\  d_1 \\  d_2 \\ \vdots \\  d_{N-3} \\  d_{N-2} \\  d_{N-1} \end{pmatrix} ,

with

.. math::
  :label: eq194

  d_0 = \frac{3}{2} \frac{y_1-y_1}{h_0}\:,  \ \ \ \ \ d_{N-1} = \frac{3}{2} \frac{y_{N-1}-y_{N-2}}{h_{N-1}}\:.

Bicubic splines
~~~~~~~~~~~~~~~

It is possible to extend the cubic spline interpolation method to
functions of two variables, that is, :math:`F(x,y)`. In this case, we
have a rectangular mesh of points given by
:math:`F_{ij} \equiv F(x_i,y_j)`. In the case of 1D splines, we needed
to store the value of the first derivative of the function at each
point, in addition to the value. In the case of *bicubic splines*, we
need to store four quantities for each mesh point:

.. math::
  :label: eq195

   \begin{aligned}
   F_{ij}    & \equiv & F(x_i, y_i)\:,             \\
   F^x_{ij}  & \equiv & \partial_x F(x_i, y_i)\:,  \\
   F^y_{ij}  & \equiv & \partial_y F(x_i, y_i)\:,  \\
   F^{xy}    & \equiv & \partial_x \partial_y F(x_i, y_i)\:. \end{aligned}

Consider the point :math:`(x,y)` at which we wish to interpolate
:math:`F`. We locate the rectangle that contains this point, such that
:math:`x_i <= x <
x_{i+1}` and :math:`y_i <= x < y_{i+1}`. Let

.. math::
  :label: eq196

   \begin{aligned}
   h & \equiv & x_{i+1}-x_i\:,  \\
   l & \equiv & y_{i+1}-y_i\:,  \\
   u & \equiv & \frac{x-x_i}{h}\:,  \\
   v & \equiv & \frac{y-y_i}{l}\:. \end{aligned}

Then, we calculate the interpolated value as

.. math::
  :label: eq197

   F(x,y) =
   \begin{pmatrix}
   p_1(u) \\ p_2(u) \\ h q_1(u) \\ h q_2(u)
   \end{pmatrix}^T
   \begin{pmatrix}({*{20}{c}})
   F_{i,j}     & F_{i+1,j}     & F^y_{i,j}      & F^y_{i,j+1}     \\
   F_{i+1,j}   & F_{i+1,j+1}   & F^y_{i+1,j}    & F^y_{i+1,j+1}   \\
   F^x_{i,j}   & F^x_{i,j+1}   & F^{xy}_{i,j}   & F^{xy}_{i,j+1}  \\
   F^x_{i+1,j} & F^x_{i+1,j+1} & F^{xy}_{i+1,j} & F^{xy}_{i+1,j+1}
   \end{pmatrix}
   \begin{pmatrix}
   p_1(v)\\ p_2(v)\\ k q_1(v) \\ k q_2(v)
   \end{pmatrix}\:.

Construction bicubic splines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We now address the issue of how to compute the derivatives that are
needed for the interpolation. The algorithm is quite simple. For every
:math:`x_i`, we perform the tridiagonal solution as we did in the 1D
splines to compute :math:`F^y_{ij}`. Similarly, we perform a tridiagonal
solve for every value of :math:`F^x_{ij}`. Finally, to compute the
cross-derivative we may *either* to the tridiagonal solve in the
:math:`y` direction of :math:`F^x_{ij}`, *or* solve in the :math:`x`
direction for :math:`F^y_{ij}` to obtain the cross-derivatives
:math:`F^{xy}_{ij}`. Hence, only minor modifications to the :math:`1D`
interpolations are necessary.

Tricubic splines
~~~~~~~~~~~~~~~~

Bicubic interpolation required two 4-component vectors and a
:math:`4 \times 4` matrix. By extension, tricubic interpolation requires
three 4-component vectors and a :math:`4 \times 4 \times 4` tensor. We
summarize the forms of these vectors in the following:

.. math::
  :label: eq198

   \begin{aligned}
   h & \equiv & x_{i+1}-x_i\:, \\
   l & \equiv & y_{i+1}-y_i\:, \\
   m & \equiv & z_{i+1}-z_i\:, \\
   u & \equiv & \frac{x-x_i}{h}\:, \\
   v & \equiv & \frac{y-y_i}{l}\:, \\
   w & \equiv & \frac{z-z_i}{m}\:.\end{aligned}

.. math::
  :label: eq199

   \begin{aligned}
   \vec{a} & = &
   \begin{pmatrix}
   p_1(u) & p_2(u) & h q_1(u) & h q_2(u)
   \end{pmatrix}^T\:, \\
   \vec{b} & = &
   \begin{pmatrix}
   p_1(v) & p_2(v) & k q_1(v) & k q_2(v)
   \end{pmatrix}^T\:, \\
   \vec{c} & = &
   \begin{pmatrix}
   p_1(w) & p_2(w) & l q_1(w) & l q_2(w)
   \end{pmatrix}^T\:. \end{aligned}

.. math::
  :label: eq200

   \begin{pmatrix}
   A_{000} = F_{i,j,k}     & A_{001}=F_{i,j,k+1}     & A_{002}=F^z_{i,j,k}      & A_{003}=F^z_{i,j,k+1}      \\
   A_{010} = F_{i,j+1,k}   & A_{011}=F_{i,j+1,k+1}   & A_{012}=F^z_{i,j+1,k}    & A_{013}=F^z_{i,j+1,k+1}    \\
   A_{020} = F^y_{i,j,k}   & A_{021}=F^y_{i,j,k+1}   & A_{022}=F^{yz}_{i,j,k}   & A_{023}=F^{yz}_{i,j,k+1}   \\
   A_{030} = F^y_{i,j+1,k} & A_{031}=F^y_{i,j+1,k+1} & A_{032}=F^{yz}_{i,j+1,k} & A_{033}=F^{yz}_{i,j+1,k+1} \\
                           &                         &                          &                            \\
   A_{100} = F_{i+1,j,k}     & A_{101}=F_{i+1,j,k+1}     & A_{102}=F^z_{i+1,j,k}      & A_{103}=F^z_{i+1,j,k+1}      \\
   A_{110} = F_{i+1,j+1,k}   & A_{111}=F_{i+1,j+1,k+1}   & A_{112}=F^z_{i+1,j+1,k}    & A_{113}=F^z_{i+1,j+1,k+1}    \\
   A_{120} = F^y_{i+1,j,k}   & A_{121}=F^y_{i+1,j,k+1}   & A_{122}=F^{yz}_{i+1,j,k}   & A_{123}=F^{yz}_{i+1,j,k+1}   \\
   A_{130} = F^y_{i+1,j+1,k} & A_{131}=F^y_{i+1,j+1,k+1} & A_{132}=F^{yz}_{i+1,j+1,k} & A_{133}=F^{yz}_{i+1,j+1,k+1} \\
                           &                         &                          &                            \\
   A_{200} = F^x_{i,j,k}      & A_{201}=F^x_{i,j,k+1}      & A_{202}=F^{xz}_{i,j,k}      & A_{203}=F^{xz}_{i,j,k+1}    \\
   A_{210} = F^x_{i,j+1,k}    & A_{211}=F^x_{i,j+1,k+1}    & A_{212}=F^{xz}_{i,j+1,k}    & A_{213}=F^{xz}_{i,j+1,k+1}  \\
   A_{220} = F^{xy}_{i,j,k}   & A_{221}=F^{xy}_{i,j,k+1}   & A_{222}=F^{xyz}_{i,j,k}     & A_{223}=F^{xyz}_{i,j,k+1}   \\
   A_{230} = F^{xy}_{i,j+1,k} & A_{231}=F^{xy}_{i,j+1,k+1} & A_{232}=F^{xyz}_{i,j+1,k}   & A_{233}=F^{xyz}_{i,j+1,k+1} \\
                           &                         &                          &                                      \\
   A_{300} = F^x_{i+1,j,k}      & A_{301}=F^x_{i+1,j,k+1}      & A_{302}=F^{xz}_{i+1,j,k}    & A_{303}=F^{xz}_{i+1,j,k+1}   \\
   A_{310} = F^x_{i+1,j+1,k}    & A_{311}=F^x_{i+1,j+1,k+1}    & A_{312}=F^{xz}_{i+1,j+1,k}  & A_{313}=F^{xz}_{i+1,j+1,k+1} \\
   A_{320} = F^{xy}_{i+1,j,k}   & A_{321}=F^{xy}_{i+1,j,k+1}   & A_{322}=F^{xyz}_{i+1,j,k}   & A_{323}=F^{xyz}_{i+1,j,k+1}  \\
   A_{330} = F^{xy}_{i+1,j+1,k} & A_{331}=F^{xy}_{i+1,j+1,k+1} & A_{332}=F^{xyz}_{i+1,j+1,k} & A_{333}=F^{xyz}_{i+1,j+1,k+1}
   \end{pmatrix}\:.

Now, we can write

.. math::
  :label: eq201

  F(x,y,z) = \sum_{i=0}^3 a_i \sum_{j=0}^3 b_j \sum_{k=0}^3 c_k \ A_{i,j,k}\:.

The appropriate derivatives of :math:`F` may be computed by a
generalization of the previous method used for bicubic splines.

Feature: B-spline orbital tiling (band unfolding)
-------------------------------------------------

In continuum QMC simulations, it is necessary to evaluate the electronic
orbitals of a system at real-space positions hundreds of millions of
times. It has been found that if these orbitals are represented in a
localized, B-spline basis, each evaluation takes a small, constant time
that is independent of system size.

Unfortunately, the memory required for storing the B-spline grows with
the second power of the system size. If we are studying perfect
crystals, however, this can be reduced to linear scaling if we *tile*
the primitive cell. In this approach, a supercell is constructed by
tiling the primitive cell :math:`N_1 \times N_2 \times N_3` in the three
lattice directions. The orbitals are then represented in real space only
in the primitive cell and an :math:`N_1 \times N_2 \times N_3` k-point
mesh. To evaluate an orbital at any point in the supercell, it is only
necessary to wrap that point back into the primitive cell, evaluate the
spline, and then multiply the phase factor,
:math:`e^{-i\mathbf{k}\cdot\mathbf{r}}`.

Here, we show that this approach can be generalized to a tiling
constructed with a :math:`3\times 3` nonsingular matrix of integers, of
which the preceding approach is a special case. This generalization
brings with it a number of advantages. The primary reason for performing
supercell calculations in QMC is to reduce finite-size errors. These
errors result from three sources: (1) the quantization of the crystal
momentum, (2) the unphysical periodicity of the exchange-correlation
(XC) hole of the electron, and (3) the kinetic-energy contribution from
the periodicity of the long-range Jastrow correlation functions. The
first source of error can be largely eliminated by twist averaging. If
the simulation cell is large enough that XC hole does not “leak” out of
the simulation cell, the second source can be eliminated either through
use of the MPC interaction or the *a postiori* correction of Chiesa et
al.

The satisfaction of the leakage requirement is controlled by whether the
minimum distance, :math:`L_{\text{min}}`, from one supercell image to
the next is greater than the width of the XC hole. Therefore, given a
choice, it is best to use a cell that is as nearly cubic as possible
since this choice maximizes :math:`L_{\text{min}}` for a given number of
atoms. Most often, however, the primitive cell is not cubic. In these
cases, if we wish to choose the optimal supercell to reduce finite-size
effects, we cannot use the simple primitive tiling scheme. In the
generalized scheme we present, it is possible to choose far better
supercells (from the standpoint of finite-size errors), while retaining
the storage efficiency of the original tiling scheme.

The mathematics
~~~~~~~~~~~~~~~

Consider the set of primitive lattice vectors,
:math:`\{\mathbf{a}^{\text{p}}_1, \mathbf{a}^{\text{p}}_2,
\mathbf{a}^{\text{p}}_3\}`. We may write these vectors in a matrix,
:math:`\mathbf{L}_p`, whose rows are the primitive lattice vectors.
Consider a nonsingular matrix of integers, :math:`\mathbf{S}`. A
corresponding set of supercell lattice vectors,
:math:`\{\mathbf{a}^{\text{s}}_1, \mathbf{a}^{\text{s}}_2, \mathbf{a}^{\text{s}}_3\}`,
can be constructed by the matrix product

.. math::
  :label: eq202

  \mathbf{a}^{\text{s}}_i = S_{ij} \mathbf{a}^{\text{p}}_j\:.

If the primitive cell contains :math:`N_p` atoms, the supercell will
then contain :math:`N_s = |\det(\mathbf{S})| N_p` atoms.

Example: FeO
~~~~~~~~~~~~

As an example, consider the primitive cell for antiferromagnetic FeO
(wustite) in the rocksalt structure.  The primitive vectors, given in
units of the lattice constant, are given by

.. math::
  :label: eq203

   \begin{aligned}
   \mathbf{a}^{\text{p}}_1 & = & \frac{1}{2}\hat{\mathbf{x}}+ \frac{1}{2}\hat{\mathbf{y}}+      \ \   \hat{\mathbf{z}}\:, \\
   \mathbf{a}^{\text{p}}_2 & = & \frac{1}{2}\hat{\mathbf{x}}+      \ \   \hat{\mathbf{y}}+ \frac{1}{2}\hat{\mathbf{z}}\:, \\
   \mathbf{a}^{\text{p}}_3 & = &   \ \      \hat{\mathbf{x}}+ \frac{1}{2}\hat{\mathbf{y}}+ \frac{1}{2}\hat{\mathbf{z}}\:. \end{aligned}

This primitive cell contains two iron atoms and two oxygen atoms. It is
a very elongated cell with acute angles and, thus, has a short minimum
distance between adjacent images.

The smallest cubic cell consistent with the AFM ordering can be
constructed with the matrix

.. math::
  :label: eq204

   \mathbf{S}= \left[\begin{array}{rrr}
     -1 & -1 &  3 \\
     -1 &  3 & -1 \\
      3 & -1 & -1
     \end{array}\right]\:.

This cell has :math:`2|\det(\mathbf{S})| = 32` iron atoms and 32 oxygen
atoms. In this example, we may perform the simulation in the 32-iron
supercell, while storing the orbitals only in the 2-iron primitive cell,
for a savings of a factor of 16.

The k-point mesh
^^^^^^^^^^^^^^^^

To be able to use the generalized tiling scheme, we need to have the
appropriate number of bands to occupy in the supercell. This may be
achieved by appropriately choosing the k-point mesh. In this section, we
explain how these points are chosen.

For simplicity, let us assume that the supercell calculation will be
performed at the :math:`\Gamma`-point. We can easily lift this
restriction later. The fact that supercell calculation is performed at
:math:`\Gamma` implies that the k-points used in the primitive-cell
calculation must be :math:`\mathbf{G}`-vectors of the superlattice. This
still leaves us with an infinite set of vectors. We may reduce this set
to a finite number by considering that the orbitals must form a linearly
independent set. Orbitals with k-vectors :math:`\mathbf{k}^p_1` and
:math:`\mathbf{k}^p_2` will differ by at most a constant factor of
:math:`\mathbf{k}^p_1 - \mathbf{k}^p_2 = \mathbf{G}^p`, where
:math:`\mathbf{G}^p` is a reciprocal lattice vector of the primitive
cell.

Combining these two considerations gives us a prescription for
generating our k-point mesh. The mesh may be taken to be the set of
k-point which are G-vectors of the superlattice, reside within the first
Brillouin zone (FBZ) of the primitive lattice, whose members do not
differ a G-vector of the primitive lattice. Upon constructing such a
set, we find that the number of included k-points is equal to
:math:`|\det(\mathbf{S})|`, precisely the number we need. This can by
considering the fact that the supercell has a volume
:math:`|\det(\mathbf{S})|` times that of the primitive cell. This
implies that the volume of the supercell’s FBZ is
:math:`|\det(\mathbf{S})|^{-1}` times that of the primitive cell. Hence,
:math:`|\det(\mathbf{S})|` G-vectors of the supercell will fit in the
FBZ of the primitive cell. Removing duplicate k-vectors, which differ
from another by a reciprocal lattice vector, avoids double-counting
vectors that lie on zone faces.

Formulae
^^^^^^^^

Let :math:`\mathbf{A}` be the matrix whose rows are the direct lattice
vectors, :math:`\{\mathbf{a}_i\}`. Then, let the matrix
:math:`\mathbf{B}` be defined as :math:`2\pi(\mathbf{A}^{-1})^\dagger`.
Its rows are the primitive reciprocal lattice vectors. Let
:math:`\mathbf{A}_p` and :math:`\mathbf{A}_s` represent the primitive
and superlattice matrices, respectively, and similarly for their
reciprocals. Then we have

.. math::
  :label: eq205

   \begin{aligned}
   \mathbf{A}_s & = & \mathbf{S}\mathbf{A}_p\:, \\
   \mathbf{B}_s & = & 2\pi\left[(\mathbf{S}\mathbf{A}_p)^{-1}\right]^\dagger\:, \\
           & = & 2\pi\left[\mathbf{A}_p^{-1} \mathbf{S}^{-1}\right]^\dagger\:, \\
           & = & 2\pi(\mathbf{S}^{-1})^\dagger (\mathbf{A}_p^{-1})^\dagger\:, \\
           & = & (\mathbf{S}^{-1})^\dagger \mathbf{B}_p\:.\end{aligned}

Consider a k-vector, :math:`\mathbf{k}`. It may alternatively be written
in basis of reciprocal lattice vectors as :math:`\mathbf{t}`.

.. math::
  :label: eq206

   \begin{aligned}
   \mathbf{k}& = & (\mathbf{t}^\dagger \mathbf{B})^\dagger\:, \\
       & = & \mathbf{B}^\dagger \mathbf{t}\:,           \\
   \mathbf{t}& = & (\mathbf{B}^\dagger)^{-1} \mathbf{k}\:,    \\
       & = & (\mathbf{B}^{-1})^\dagger \mathbf{k}\:,    \\
       & = & \frac{\mathbf{A}\mathbf{k}}{2\pi}\:.\end{aligned}

We may then express a twist vector of the primitive lattice,
:math:`\mathbf{t}_p`, in terms of the superlattice.

.. math::
  :label: eq207

   \begin{aligned}
   \mathbf{t}_s & = & \frac{\mathbf{A}_s \mathbf{k}}{2\pi}\:,                           \\
         & = & \frac{\mathbf{A}_s \mathbf{B}_p^\dagger \mathbf{t}_p}{2\pi}\:,         \\
         & = & \frac{\mathbf{S}\mathbf{A}_p \mathbf{B}_p^\dagger \mathbf{t}_p}{2\pi}\:,   \\
         & = & \frac{2\pi \mathbf{S}\mathbf{A}_p \mathbf{A}_p^{-1} \mathbf{t}_p}{2\pi}\:, \\
         & = & \mathbf{S}\mathbf{t}_p\:.\end{aligned}

This gives the simple result that twist vectors transform in precisely
the same way as direct lattice vectors.

Feature: Hybrid orbital representation
--------------------------------------

.. _ Written by Kenneth P. Esler, Jr.
  Document originally included in QMCPACK at src/QMCWaveFunctions/AtomicOrbital.tex
  Originally titled "Hybrid orbital representation"

.. math::
  :label: eq208

   \phi(\mathbf{r}) = \sum_{\ell=0}^{\ell_\text{max}} \sum_{m=-\ell}^\ell Y_\ell^m (\hat{\Omega})
   u_{\ell m}(r)\:,

where :math:`u_{lm}(r)` are complex radial functions represented in some
radial basis (e.g., splines).

Real spherical harmonics
~~~~~~~~~~~~~~~~~~~~~~~~

If :math:`\phi(\mathbf{r})` can be written as purely real, we can change
the representation so that

.. math::
  :label: eq209

   \phi(\mathbf{r}) = \sum_{l=0}^{l_\text{max}} \sum_{m=-\ell}^\ell Y_{\ell m}(\hat{\Omega})
   \bar{u}_{lm}(r)\:,

where :math:`\bar{Y}_\ell^m` are the *real* spherical harmonics defined
by

.. math::
  :label: eq210

   Y_{\ell m} = \begin{cases}
   Y_\ell^0 & \mbox{if } m=0\\
   {1\over 2}\left(Y_\ell^m+(-1)^m \, Y_\ell^{-m}\right) \ = \rm Re\left[Y_\ell^m\right]
   %\sqrt{2} N_{(\ell,m)} P_\ell^m(\cos \theta) \cos m\varphi
   & \mbox{if } m>0 \\
   {1\over i 2}\left(Y_\ell^{-m}-(-1)^{m}\, Y_\ell^{m}\right) = \rm Im\left[Y_\ell^{-m}\right]
   %\sqrt{2} N_{(\ell,m)} P_\ell^{-m}(\cos \theta) \sin m\varphi
   &\mbox{if } m<0\:.
   \end{cases}


We need then to relate :math:`\bar{u}_{\ell m}` to :math:`u_{\ell m}`.
We wish to express

.. math::
  :label: eq211

   \rm Re\left[\phi(\mathbf{r})\right] = \sum_{\ell=0}^{\ell_\text{max}} \sum_{m=-\ell}^\ell
   \rm Re\left[Y_\ell^m (\hat{\Omega}) u_{\ell m}(r)\right]

in terms of :math:`\bar{u}_{\ell m}(r)` and :math:`Y_{\ell m}`.

.. math::
  :label: eq212

   \begin{aligned}
   \rm Re\left[Y_\ell^m u_{\ell m}\right] & = & \rm Re\left[Y_\ell^m\right]
   \rm Re\left[u_{\ell m}\right] - \rm Im\left[Y_\ell^m\right] \rm Im\left[u_{\ell m}\right]\:.\end{aligned}

For :math:`m>0`,

.. math::
  :label: eq213

  \rm Re\left[Y_\ell^m\right] = Y_{\ell m} \qquad \text{and} \qquad \rm Im\left[Y_\ell^m\right] = Y_{\ell\,-m}\:.

For :math:`m<0`,

.. math::
  :label: eq214

  \rm Re\left[Y_\ell^m\right] = (-1)^m Y_{\ell\, -m} \qquad \text and \qquad \rm Im\left[Y_\ell^m\right] = -(-1)^m Y_{\ell m}\:.

Then for :math:`m > 0`,

.. math::
  :label: eq215

   \begin{aligned}
   \bar{u}_{\ell m} & = & \rm Re\left[u_{\ell m}\right] + (-1)^m \rm Re\left[u_{\ell\,-m}\right]\:, \\
   \bar{u}_{\ell\, -m} & = & -\rm Im\left[u_{\ell m}\right] + (-1)^m \rm Im\left[u_{\ell\,-m}\right]\:.\end{aligned}


Projecting to atomic orbitals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Inside a muffin tin, orbitals are represented as products of spherical
harmonics and 1D radial functions, primarily represented by splines. For
a muffin tin centered at :math:`\mathbf{I}`,

.. math::
  :label: eq216

   \phi_n(\mathbf{r}) = \sum_{\ell,m} Y_\ell^m(\hat{\mathbf{r}-\mathbf{I}})
   u_{lm}\left(\left|\mathbf{r}- \mathbf{I}\right|\right) \:.

Let use consider the case that our original representation for
:math:`\phi(\mathbf{r})` is of the form

.. math::
  :label: eq217

  \phi_{n,\mathbf{k}}(\mathbf{r}) = \sum_\mathbf{G}c_{\mathbf{G}+\mathbf{k}}^n e^{i(\mathbf{G}+ \mathbf{k})\cdot \mathbf{r}}\:.

Recall that

.. math::
  :label: eq218

   e^{i\mathbf{k}\cdot\mathbf{r}} = 4\pi \sum_{\ell,m} i^\ell j_\ell(|\mathbf{r}||\mathbf{k}|)
   Y_\ell^m(\hat{\mathbf{k}}) \left[Y_\ell^m(\hat{\mathbf{r}})\right]^*\:.

Conjugating,

.. math::
  :label: eq219

   e^{-i\mathbf{k}\cdot\mathbf{r}} = 4\pi\sum_{\ell,m} (-i)^\ell j_\ell(|\mathbf{r}||\mathbf{k}|)
   \left[Y_\ell^m(\hat{\mathbf{k}})\right]^* Y_\ell^m(\hat{\mathbf{r}})\:.

Setting :math:`\mathbf{k}\rightarrow -k`,

.. math::
  :label: eq220

   e^{i\mathbf{k}\cdot\mathbf{r}} = 4\pi\sum_{\ell,m} i^\ell j_\ell(|\mathbf{r}||\mathbf{k}|)
   \left[Y_\ell^m(\hat{\mathbf{k}})\right]^* Y_\ell^m(\hat{\mathbf{r}})\:.

Then,

.. math::
  :label: eq221

   e^{i\mathbf{k}\cdot(\mathbf{r}-\mathbf{I})} = 4\pi\sum_{\ell,m} i^\ell j_\ell(|\mathbf{r}-\mathbf{I}||\mathbf{k}|)
   \left[Y_\ell^m(\hat{\mathbf{k}})\right]^* Y_\ell^m(\hat{\mathbf{r}-\mathbf{I}})\:.

.. math::
  :label: eq222

   e^{i\mathbf{k}\cdot\mathbf{r}} = 4\pi e^{i\mathbf{k}\cdot\mathbf{I}} \-\sum_{\ell,m} i^\ell j_\ell(|\mathbf{r}-\mathbf{I}||\mathbf{k}|)
   \left[Y_\ell^m(\hat{\mathbf{k}})\right]^* Y_\ell^m(\hat{\mathbf{r}-\mathbf{I}})\:.

Then

.. math::
  :label: eq223

   \phi_{n,\mathbf{k}}(\mathbf{r}) =  \sum_\mathbf{G}4\pi c_{\mathbf{G}+\mathbf{k}}^n
   e^{i(\mathbf{G}+\mathbf{k})\cdot\mathbf{I}} \sum_{\ell,m}
     i^\ell j_\ell(|\mathbf{G}+\mathbf{k}||\mathbf{r}-\mathbf{I}|)
     \left[Y_\ell^m(\hat{\mathbf{G}+\mathbf{k}})\right]^*
   Y_\ell^m(\hat{\mathbf{r}- \mathbf{I}})\:.

Comparing with :eq:`eq216`,

.. math::
  :label: eq224

   u_{\ell m}^n(r) = 4\pi i^\ell \sum_G c_{\mathbf{G}+\mathbf{k}}^n e^{i(\mathbf{G}+\mathbf{k})\cdot\mathbf{I}}  j_\ell\left(|\mathbf{G}+ \mathbf{k}|r|\right)
   \left[Y_\ell^m(\hat{\mathbf{G}+ \mathbf{k}})\right]^*\:.

If we had adopted the opposite sign convention for Fourier transforms
(as is unfortunately the case in wfconvert), we would have

.. math::
  :label: eq225

   u_{\ell m}^n(r) = 4\pi (-i)^\ell \sum_G c_{\mathbf{G}+\mathbf{k}}^n e^{-i(\mathbf{G}+\mathbf{k})\cdot\mathbf{I}}  j_\ell\left(|\mathbf{G}+ \mathbf{k}|r|\right)
   \left[Y_\ell^m(\hat{\mathbf{G}+ \mathbf{k}})\right]^*\:.


Feature: Electron-electron-ion Jastrow factor
---------------------------------------------

.. _ Written by Kenneth P. Esler, Jr.
  Document originally included in QMCPACK at src/QMCWaveFunctions/Jastrow/eeI_Jastrow.tex
  Originally titled ``Electron-electron-ion Jastrow factor''

The general form of the 3-body Jastrow we describe here depends on the
three interparticle distances, :math:`(r_{ij}, r_{iI}, r_{jI})`.

.. math::
  :label: eq226

   J_3 = \sum_{I\in\text{ions}} \sum_{i,j \in\text{elecs};i\neq j} U(r_{ij}, r_{iI},
   r_{jI})\:.

Note that we constrain the form of :math:`U` such that
:math:`U(r_{ij}, r_{iI},r_{jI}) = U(r_{ij}, r_{jI},r_{iI})` to preserve
the particle symmetry of the wavefunction. We then compute the gradient
as

.. math::
  :label: eq227

   \nabla_i J_3 =  \sum_{I\in\text{ions}} \sum_{j \neq i}
   \left[\frac{\partial U(r_{ij}, r_{iI},r_{jI})}{\partial r_{ij}}
     \frac{\mathbf{r}_i - \mathbf{r}_j}{|\mathbf{r}_i - \mathbf{r}_j|}
   + \frac{\partial U(r_{ij}, r_{iI},r_{jI})}{\partial r_{iI}}
     \frac{\mathbf{r}_i - \mathbf{I}}{|\mathbf{r}_i - \mathbf{I}|}  \right]\:.

To compute the Laplacian, we take

.. math::
  :label: eq228

   \begin{aligned}
   \nabla_i^2 J_3 & = & \nabla_i \cdot \left(\nabla_i J_3\right)\:, \\
   & = & \sum_{I\in\text{ions}} \sum_{j\neq i } \left[
   \frac{\partial^2 U}{\partial r_{ij}^2} + \frac{2}{r_{ij}} \frac{\partial
     U}{\partial r_{ij}} + 2 \frac{\partial^2 U}{\partial r_{ij}\partial
     r_{iI}}\frac{\mathbf{r}_{ij}\cdot\mathbf{r}_{iI}}{r_{ij}r_{iI}} +\frac{\partial^2 U}{\partial
     r_{iI}^2}
   + \frac{2}{r_{iI}}\frac{\partial U}{\partial r_{iI}} \nonumber
   \right]\:.\end{aligned}

We now wish to compute the gradient of these terms w.r.t. the ion
position, :math:`I`.

.. math::
  :label: eq229

   \nabla_I J_3 = -\sum_{j\neq i} \left[ \frac{\partial U(r_{ij}, r_{iI},r_{jI})}{\partial r_{iI}}
     \frac{\mathbf{r}_i - \mathbf{I}}{|\mathbf{r}_i - \mathbf{I}|}
   +\frac{\partial U(r_{ij}, r_{iI},r_{jI})}{\partial r_{jI}}
     \frac{\mathbf{r}_j - \mathbf{I}}{|\mathbf{r}_j - \mathbf{I}|} \right]\:.

For the gradient w.r.t. :math:`i` of the gradient w.r.t. :math:`I`, the
result is a tensor:

.. math::
  :label: eq230

   \begin{aligned}
   \nabla_I \nabla_i J_3 & = & \nabla_I \sum_{j \neq i}
   \left[\frac{\partial U(r_{ij}, r_{iI},r_{jI})}{\partial r_{ij}}
     \frac{\mathbf{r}_i - \mathbf{r}_j}{|\mathbf{r}_i - \mathbf{r}_j|}
   + \frac{\partial U(r_{ij}, r_{iI},r_{jI})}{\partial r_{iI}}
     \frac{\mathbf{r}_i - \mathbf{I}}{|\mathbf{r}_i - \mathbf{I}|}  \right]\:, \\\nonumber \\\nonumber
   & = & -\sum_{j\neq i} \left[
   \frac{\partial^2 U}{\partial r_{ij}r_{iI}} \hat{\mathbf{r}}_{ij} \otimes
   \hat{\mathbf{r}}_{iI} + \left(\frac{\partial^2 U}{\partial r_{iI}^2} -
   \frac{1}{r_{iI}} \frac{\partial U}{\partial r_{iI}}\right)
   \hat{\mathbf{r}}_{iI} \otimes \hat{\mathbf{r}}_{iI} \right. + \\\nonumber
   & & \left. \qquad \ \ \  \frac{\partial^U}{\partial r_{ij}r_{jI}} \hat{\mathbf{r}}_{ij} \otimes \hat{\mathbf{r}}_{jI} + \frac{\partial^2 U}{\partial r_{iI}\partial r_{jI}}
   \hat{\mathbf{r}}_{iI}\otimes \hat{\mathbf{r}}_{jI}  +
   \frac{1}{r_{iI}} \frac{\partial U}{\partial r_{iI}} \overleftrightarrow{\mathbf{1}}\right]\:.\end{aligned}

.. math::
  :label: eq231

   \begin{aligned}
   \nabla_I \nabla_i J_3 & = & \nabla_I \sum_{j \neq i}
   \left[\frac{\partial U(r_{ij}, r_{iI},r_{jI})}{\partial r_{ij}}
     \frac{\mathbf{r}_i - \mathbf{r}_j}{|\mathbf{r}_i - \mathbf{r}_j|}
   + \frac{\partial U(r_{ij}, r_{iI},r_{jI})}{\partial r_{iI}}
     \frac{\mathbf{r}_i - \mathbf{I}}{|\mathbf{r}_i - \mathbf{I}|}  \right]\:, \\\nonumber
   & = & \sum_{j\neq i} \left[ -\frac{\partial^2 U}{\partial r_{ij}\partial r_{iI}} \hat{\mathbf{r}}_{ij} \otimes \hat{\mathbf{r}}_{iI} +
   \left(-\frac{\partial^2 U}{\partial r_{iI}^2}  + \frac{1}{r_{iI}}\frac{\partial U}{\partial r_{iI}} \right)
   \hat{\mathbf{r}}_{iI} \otimes \hat{\mathbf{r}}_{iI} - \frac{1}{r_{iI}}\frac{\partial U}{\partial r_{iI}} \overleftrightarrow{\mathbf{1}}
   \right]\:.\end{aligned}

For the Laplacian,

.. math::
  :label: eq232

   \begin{aligned}
   \nabla_I \nabla_i^2 J_3 & = & \nabla_I\left[\nabla_i \cdot \left(\nabla_i J_3\right)\right]\:, \\
   & = & \nabla_I \sum_{j\neq i } \left[
   \frac{\partial^2 U}{\partial r_{ij}^2} + \frac{2}{r_{ij}} \frac{\partial
     U}{\partial r_{ij}} + 2 \frac{\partial^2 U}{\partial r_{ij}\partial
     r_{iI}}\frac{\mathbf{r}_{ij}\cdot\mathbf{r}_{iI}}{r_{ij}r_{iI}} +\frac{\partial^2 U}{\partial
     r_{iI}^2}
   + \frac{2}{r_{iI}}\frac{\partial U}{\partial r_{iI}} \nonumber
   \right]\:, \\
   & = & \sum_{j\neq i }
   \left[ \frac{\partial^3 U}{\partial r_{iI} \partial^2 r_{ij}} +
   \frac{2}{r_{ij}} \frac{\partial^2 U}{\partial r_{iI} \partial r_{ij}}
   + 2\left(\frac{\partial^3 U}{\partial r_{ij}\partial^2 r_{iI}} -\frac{1}{r_{iI}} \frac{\partial^2 U}{\partial r_{ij}\partial r_{iI}}\right)\frac{\mathbf{r}_{ij}\cdot\mathbf{r}_{iI}}{r_{ij}r_{iI}} + \frac{\partial^3 U}{\partial^3 r_{iI}} - \frac{2}{r_{iI}^2} \frac{\partial U}{ \partial r_{iI}} + \frac{2}{r_{iI}} \frac{\partial^2 U}{\partial^2 r_{iI}}
   \right] \frac{\mathbf{I} - \mathbf{r}_i}{|\mathbf{I} - \mathbf{r}_i|} + \nonumber \\\nonumber
    & & \sum_{j\neq i } \left[ \frac{\partial^3U}{\partial r_{ij}^2 \partial r_{jI}} + \frac{2}{r_{ij}}\frac{\partial^2 U}{\partial r_{jI}\partial r_{ij}}
   + 2\frac{\partial^3 U}{\partial r_{ij}\partial r_{iI}\partial r_{jI}}\frac{\mathbf{r}_{ij}\cdot\mathbf{r}_{iI}}{r_{ij}r_{iI}}
   +\frac{\partial^3 U}{\partial r_{iI}^2 \partial r_{jI}} + \frac{2}{r_{iI}}\frac{\partial^2 U}{\partial r_{iI}\partial r_{jI}} \right]
   \frac{\mathbf{I} - \mathbf{r}_j}{|\mathbf{r}_j - \mathbf{I}|} + \\\nonumber
   & & \sum_{j\neq i } \left[ -\frac{2}{r_{iI}}\frac{\partial^2 U}{\partial r_{ij}\partial r_{iI}}\right] \frac{\mathbf{r}_{ij}}{r_{ij}}\:.\end{aligned}

.. _feature-kspace-jastrow:

Feature: Reciprocal-space Jastrow factors
-----------------------------------------

.. _ Written by Kenneth P. Esler, Jr.
  Document originally included in QMCPACK at src/QMCWaveFunctions/Jastrow/kSpaceJastrowNotes.tex
  Originally titled "Notes on Reciprocal-Space Jastrow Factors"

Two-body Jastrow
~~~~~~~~~~~~~~~~

.. math::
  :label: eq233

  J_2 = \sum_{\mathbf{G}\neq \mathbf{0}}\sum_{i\neq j} a_\mathbf{G}e^{i\mathbf{G}\cdot(\mathbf{r}_i-\mathbf{r}_j)}\:.

This may be rewritten as

.. math::
  :label: eq234

   \begin{aligned}
   J_2 & = & \sum_{\mathbf{G}\neq \mathbf{0}}\sum_{i\neq j} a_\mathbf{G}e^{i\mathbf{G}\cdot\mathbf{r}_i}e^{-i\mathbf{G}\cdot\mathbf{r}_j}\:, \\
   & = & \sum_{\mathbf{G}\neq \mathbf{0}} a_\mathbf{G}\left\{
   \underbrace{\left[\sum_i e^{i\mathbf{G}\cdot\mathbf{r}_i} \right]}_{\rho_\mathbf{G}}
   \underbrace{\left[\sum_j e^{-i\mathbf{G}\cdot\mathbf{r}_j} \right]}_{\rho_{-\mathbf{G}}}  -1 \right\}\:.\end{aligned}

The :math:`-1` is just a constant term and may be subsumed into the
:math:`a_\mathbf{G}` coefficient by a simple redefinition. This leaves a
simple, but general, form:

.. math::
  :label: eq235

  J_2 = \sum_{\mathbf{G}\neq\mathbf{0}} a_\mathbf{G}\rho_\mathbf{G}\rho_{-\mathbf{G}}\:.

We may now further constrain this on physical grounds. First, we
recognize that :math:`J_2` should be real. Since
:math:`\rho_{-\mathbf{G}} =
\rho_\mathbf{G}^*`, it follows that
:math:`\rho_{\mathbf{G}}\rho_{-\mathbf{G}} = |\rho_\mathbf{G}|^2` is
real, so that :math:`a_\mathbf{G}` must be real. Furthermore, we group
the :math:`\mathbf{G}`\ ’s into :math:`(+\mathbf{G}, -\mathbf{G})` pairs
and sum over only the positive vectors to save time.

One-body Jastrow
~~~~~~~~~~~~~~~~

The 1-body Jastrow has a similar form but depends on the displacement
from the electrons to the ions in the system.

.. math::
  :label: eq236

   J_1 = \sum_{\mathbf{G}\neq\mathbf{0}} \sum_{\alpha}
   \sum_{i\in\mathbf{I}^\alpha}\sum_{j\in\text{elec.}} b^{\alpha}_\mathbf{G}
     e^{i\mathbf{G}\cdot(\mathbf{I}^{\alpha}_i - \mathbf{r}_j)}\:,

where :math:`\alpha` denotes the different ionic species. We may rewrite
this in terms of :math:`\rho^{\alpha}_\mathbf{G}`:

.. math::
  :label: eq237

   J_1 = \sum_{\mathbf{G}\neq\mathbf{0}} \left[\sum_\alpha b^\alpha_\mathbf{G}
     \rho_\mathbf{G}^\alpha\right] \rho_{-\mathbf{G}}\:,

where

.. math::
  :label: eq238

  \rho^\alpha_\mathbf{G}= \sum_{i\in\mathbf{I}^\alpha} e^{i\mathbf{G}\cdot\mathbf{I}^\alpha_i}\:.

We note that in the preceding equation, for a single configuration of
the ions, the sum in brackets can be rewritten as a single constant.
This implies that the per-species 1-body coefficients,
:math:`b^\alpha_\mathbf{G}`, are underdetermined for single
configuration of the ions. In general, if we have :math:`N` species, we
need :math:`N` linearly independent ion configurations to uniquely
determine :math:`b^{\alpha}_\mathbf{G}`. For this reason, we will drop
the :math:`\alpha` superscript of :math:`b_\mathbf{G}` for now.

If we do desire to find a reciprocal space 1-body Jastrow that is
transferable to systems with different ion positions and :math:`N` ionic
species, we must perform compute :math:`b_\mathbf{G}` for :math:`N`
different ion configurations. We may then construct :math:`N` equations
at each value of :math:`\mathbf{G}` to solve for the :math:`N` unknown
values, :math:`b^\alpha_\mathbf{G}`.

In the 2-body case, :math:`a_\mathbf{G}` was constrained to be real by
the fact that :math:`\rho_\mathbf{G}\rho_{-\mathbf{G}}` was real.
However, in the 1-body case, there is no such guarantee about
:math:`\rho^\alpha_\mathbf{G}\rho_\mathbf{G}`. Therefore, in general,
:math:`b_\mathbf{G}` may be complex.

Symmetry considerations
~~~~~~~~~~~~~~~~~~~~~~~

For a crystal, many of the :math:`\mathbf{G}`-vectors will be equivalent
by symmetry. It is useful then to divide the :math:`\mathbf{G}`-vectors
into symmetry-related groups and then require that they share a common
coefficient. Two vectors, :math:`\mathbf{G}` and :math:`\mathbf{G}'`,
may be considered to be symmetry related if, for all :math:`\alpha` and
:math:`\beta`

.. math::
  :label: eq239

  \rho^\alpha_\mathbf{G}\rho^\beta_{-\mathbf{G}} = \rho^\alpha_{\mathbf{G}'} \rho^\beta_{-\mathbf{G}'}\:.

For the 1-body term, we may also omit from our list of
:math:`\mathbf{G}`-vectors those for which all species structure factors
are zero. This is equivalent to saying that if we are tiling a primitive
cell we should include only the :math:`\mathbf{G}`-vectors of the
primitive cell and not the supercell. Note that this is not the case for
the 2-body term since the XC hole should not have the periodicity of the
primitive cell.

Gradients and Laplacians
~~~~~~~~~~~~~~~~~~~~~~~~

.. math::
  :label: eq240

   \begin{aligned}
   \nabla_{\mathbf{r}_i} J_2 & = & \sum_{\mathbf{G}\neq 0} a_\mathbf{G}\left[\left(\nabla_{\mathbf{r}_i}\rho_\mathbf{G}\right) \rho_{-\mathbf{G}} + \text{c.c.}\right]\:, \\
   & = & \sum_{\mathbf{G}\neq \mathbf{0}} 2\mathbf{G}a_\mathbf{G}\mathbf{Re}\left(i e^{i\mathbf{G}\cdot\mathbf{r}_i} \rho_{-\mathbf{G}} \right)\:, \\
   & = & \sum_{\mathbf{G}\neq \mathbf{0}} -2\mathbf{G}a_\mathbf{G}\mathbf{Im}\left(e^{i\mathbf{G}\cdot\mathbf{r}_i} \rho_{-\mathbf{G}} \right)\:.\end{aligned}

The Laplacian is then given by

.. math::
  :label: eq241

   \begin{aligned}
     \nabla^2 J_2 & = & \sum_{\mathbf{G}\neq\mathbf{0}} a_\mathbf{G}\left[\left(\nabla^2 \rho_\mathbf{G}\right) \rho_{-\mathbf{G}} + \text{c.c.}
     + 2\left(\nabla \rho_\mathbf{G})\cdot(\nabla \rho_{-\mathbf{G}}\right)\right]\:, \\
   & = & \sum_{\mathbf{G}\neq\mathbf{0}} a_\mathbf{G}\left[ -2G^2\mathbf{Re}(e^{i\mathbf{G}\cdot\mathbf{r}_i}\rho_{-\mathbf{G}}) +
       2\left(i\mathbf{G}e^{i\mathbf{G}\cdot\mathbf{r}_i}\right) \cdot \left(-i\mathbf{G}e^{-i\mathbf{G}\cdot\mathbf{r}_i}\right)
   \right]\:, \\
   & = & 2 \sum_{\mathbf{G}\neq\mathbf{0}} G^2 a_\mathbf{G}\left[-\mathbf{Re}\left(e^{i\mathbf{G}\cdot\mathbf{r}_i}\rho_{-\mathbf{G}}\right) + 1\right]\:. \end{aligned}

.. bibliography:: /bibs/design_features.bib
