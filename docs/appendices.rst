.. _appendices:

Appendices
==========

.. _appendix-a:

Appendix A: Derivation of twist averaging efficiency
----------------------------------------------------

In this appendix we derive the relative statistical efficiency of twist
averaging with an irreducible (weighted) set of k-points versus using
uniform weights over an unreduced set of k-points (e.g., a full
Monkhorst-Pack mesh).

Consider the weighted average of a set of statistical variables
:math:`\{x_m\}` with weights :math:`\{w_m\}`:

.. math::
  :label: eq266

   \begin{aligned}
     x_{TA} = \frac{\sum_mw_mx_m}{\sum_mw_m}\:.\end{aligned}

If produced by a finite QMC run at a set of twist angles/k-points
:math:`\{k_m\}`, each variable mean :math:`\langle{x_m}\rangle` has a statistical
error bar :math:`\sigma_m`, and we can also obtain the statistical error
bar of the mean of the twist-averaged quantity :math:`\langle{x_{TA}\rangle}`:

.. math::
  :label: eq267

   \begin{aligned}
     \sigma_{TA} = \frac{\left(\sum_mw_m^2\sigma_m^2\right)^{1/2}}{\sum_mw_m}\:.\end{aligned}

The error bar of each individual twist :math:`\sigma_m` is related to
the autocorrelation time :math:`\kappa_m`, intrinsic variance
:math:`v_m`, and the number of postequilibration MC steps
:math:`N_{step}` in the following way:

.. math::
  :label: eq268

   \begin{aligned}
     \sigma_m^2=\frac{\kappa_mv_m}{N_{step}}\:.\end{aligned}

In the setting of twist averaging, the autocorrelation time and variance
for different twist angles are often very similar across twists, and we
have

.. math::
  :label: eq269

   \begin{aligned}
     \sigma_m^2=\sigma^2=\frac{\kappa v}{N_{step}}\:.\end{aligned}

If we define the total weight as :math:`W`, that is,
:math:`W\equiv\sum_{m=1}^Mw_m`, for the weighted case with :math:`M`
irreducible twists, the error bar is

.. math::
  :label: eq270

   \begin{aligned}
     \sigma_{TA}^{weighted}=\frac{\left(\sum_{m=1}^Mw_m^2\right)^{1/2}}{W}\sigma\:.\end{aligned}

For uniform weighting with :math:`w_m=1`, the number of twists is
:math:`W` and we have

.. math::
  :label: eq271

   \begin{aligned}
     \sigma_{TA}^{uniform}=\frac{1}{\sqrt{W}}\sigma\:.\end{aligned}

We are interested in comparing the efficiency of choosing weights
uniformly or based on the irreducible multiplicity of each twist angle
for a given target error bar :math:`\sigma_{target}`. The number of MC
steps required to reach this target for uniform weighting is

.. math::
  :label: eq272

   \begin{aligned}
     N_{step}^{uniform} = \frac{1}{W}\frac{\kappa v}{\sigma_{target}^2}\:,\end{aligned}

while for nonuniform weighting we have

.. math::
  :label: eq273

   \begin{aligned}
     N_{step}^{weighted} &= \frac{\sum_{m=1}^Mw_m^2}{W^2}\frac{\kappa v}{\sigma_{target}^2} \nonumber\:,\\
                     &=\frac{\sum_{m=1}^Mw_m^2}{W}N_{step}^{uniform}\:.\end{aligned}

The MC efficiency is defined as

.. math::
  :label: eq274

   \begin{aligned}
     \xi = \frac{1}{\sigma^2t}\:,\end{aligned}

where :math:`\sigma` is the error bar and :math:`t` is the total CPU
time required for the MC run.

The main advantage made possible by irreducible twist weighting is to
reduce the equilibration time overhead by having fewer twists and,
hence, fewer MC runs to equilibrate. In the context of twist averaging,
the total CPU time for a run can be considered to be

.. math::
  :label: eq275

   \begin{aligned}
     t=N_{twist}(N_{eq}+N_{step})t_{step}\:,\end{aligned}

where :math:`N_{twist}` is the number of twists, :math:`N_{eq}` is the
number of MC steps required to reach equilibrium, :math:`N_{step}` is
the number of MC steps included in the statistical averaging as before,
and :math:`t_{step}` is the wall clock time required to complete a
single MC step. For uniform weighting :math:`N_{twist}=W`; while for
irreducible weighting :math:`N_{twist}=M`.

We can now calculate the relative efficiency (:math:`\eta`) of
irreducible vs. uniform twist weighting with the aim of obtaining a
target error bar :math:`\sigma_{target}`:

.. math::
  :label: eq276

   \begin{aligned}
     \eta &= \frac{\xi_{TA}^{weighted}}{\xi_{TA}^{uniform}} \nonumber\:, \\
          &= \frac{\sigma_{target}^2t_{TA}^{uniform}}{\sigma_{target}^2t_{TA}^{weighted}} \nonumber\:, \\
          &= \frac{W(N_{eq}+N_{step}^{uniform})}{M(N_{eq}+N_{step}^{weighted})} \nonumber\:, \\
          &= \frac{W(N_{eq}+N_{step}^{uniform})}{M(N_{eq}+\frac{\sum_{m=1}^Mw_m^2}{W}N_{step}^{uniform})} \nonumber\:, \\
          &= \frac{W}{M}\frac{1+f}{1+\frac{\sum_{m=1}^Mw_m^2}{W}f}\:.\end{aligned}

In this last expression, :math:`f` is the ratio of the number of usable
MC steps to the number that must be discarded during equilibration
(:math:`f=N_{step}^{uniform}/N_{eq}`); and as before,
:math:`W=\sum_mw_m`, which is the number of twist angles in the uniform
weighting case. It is important to recall that
:math:`N_{step}^{uniform}` in :math:`f` is defined relative to uniform
weighting and is the number of MC steps required to reach a target
accuracy in the case of uniform twist weights.

The formula for :math:`\eta` in the preceding can be easily changed with
the help of :eq:`eq273` to reflect the number of MC
steps obtained in an irreducibly weighted run instead. A good exercise
is to consider runs that have already completed with either uniform or
irreducible weighting and calculate the expected efficiency change had
the opposite type of weighting been used.

The break even point :math:`(\eta=1)` can be found at a usable step
fraction of

.. math::
  :label: eq277

   \begin{aligned}
     f=\frac{W-M}{M\frac{\sum_{m=1}^Mw_m^2}{W}-W}\:.\end{aligned}

The relative efficiency :math:`(\eta)` is useful to consider in view of
certain scenarios. An important case is where the number of required
sampling steps is no larger than the number of equilibration steps
(i.e., :math:`f\approx 1`). For a very simple case with eight uniform
twists with irreducible multiplicities of :math:`w_m\in\{1,3,3,1\}`
(:math:`W=8`, :math:`M=4`), the relative efficiency of irreducible vs.
uniform weighting is
:math:`\eta=\frac{8}{4}\frac{2}{1+20/8}\approx 1.14`. In this case,
irreducible weighting is about :math:`14`\ % more efficient than uniform
weighting.

Another interesting case is one in which the number of sampling steps
you can reach with uniform twists before wall clock time runs out is
small relative to the number of equilibration steps
(:math:`f\rightarrow 0`). In this limit, :math:`\eta\approx W/M`. For
our eight-uniform-twist example, this would result in a relative
efficiency of :math:`\eta=8/4=2`, making irreducible weighting twice as
efficient.

A final case of interest is one in which the equilibration time is short
relative to the available sampling time :math:`(f\rightarrow\infty)`,
giving :math:`\eta\approx W^2/(M\sum_{m=1}^Mw_m^2)`. Again, for our
simple example we find :math:`\eta=8^2/(4\times 20)\approx 0.8`, with
uniform weighting being :math:`25`\ % more efficient than irreducible
weighting. For this example, the crossover point for irreducible
weighting being more efficient than uniform weighting is :math:`f<2`,
that is, when the available sampling period is less than twice the
length of the equilibration period. The expected efficiency ratio and
crossover point should be checked for the particular case under
consideration to inform the choice between twist averaging methods.
