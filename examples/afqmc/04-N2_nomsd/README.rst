Example 4: NOMSD Trial
----------------------

In this example we will show how to format trial different wavefunctions in such a way
that qmcpack can read them.

Rather than use the `pyscf_to_afqmc.py`, script we will break up the process to allow for
more flexibility and show what is going on under the hood.

The qmcpack input can be generated with the scf.py script. See the comments in scf.py for
a breakdown of the steps involved.

Currently QMCPACK can deal with trial wavefunctions in two forms: Non-orthogonal multi
slater determinant trial wavefunctions (NOMSD) and particle-hole style trial wavefunctions
(PHMSD). The NOMSD trial wavefunctions are the most general form and expect Slater
determinants in the form of M X N matrices of molecular orbital coefficients, where N is
the number of electrons and M is the number of orbitals, along with a list of ci
coefficients.  Importantly the Slater determinants must be non-orthogonal.
