Nexus QE Example 2: Self-consistent DFT+U+V calculations on diamond primitive cell (DFT only)
=============================================================

In this example,we show how to perform self-consistent DFT+U+V calculations 
using the updated Hubbard card input in QE>7.1. 

The new Hubbard input card in the QE input file of this example is the following:
HUBBARD atomic
V C-2p C-2p 1 2 1e-08

Detailed information about the syntax here is provided atomic
https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm1745

Basically, this input initializes +V values near zero while the integers provide the pair
atomic index for the V parameter for the zeroth step of the self-consistent U calculations.
The "atomic" manifold is also used to project the hubbard operators. Ideally, the projection 
manifold need to be explored as it can lead to significantly different magnetic moments. 

In the nexus script this is simply defined as:
hubbard      = {'V' : {('C-2p', 'C-2p'): 1e-8}}

Here, nexus considers all the atomic pair combinations automatically. However, if the user
wants to define the pairs specifically, then in nexus it can be set up as:
hubbard      = {'V' : {('C-2p', 'C-2p'): [{'indices':(1,1), 'value':1e-8},
                                          {'indices':(1,2), 'value':1e-8}
                                          ]}}

which will produce the Hubbard card in the QE input file as:
HUBBARD atomic
V C-2p C-2p 1 1 1e-08
V C-2p C-2p 1 2 1e-08

The nexus input described here only works between the interactions of standard orbitals. 
Interactions involving background orbitals are not fully supported. For more info:
https://www.quantum-espresso.org/Doc/user_guide_PDF/Hubbard_input.pdf


After the initial SCF calculation, hp.x runs the linear response calculations. 
All the parameters are supported in https://www.quantum-espresso.org/Doc/INPUT_HP.html. 

The nexus input of hp is the following:

hp = generate_hp(
    nq1          = 2,
    nq2          = 2,
    nq3          = 2,
    lmin         = 0,
    job          = job(cores=16,app='/Users/ksu/Software/qe-7.1/build_mpi/bin/hp.x'),
    path         = 'diamond/scf_step_{}'.format(step),
    dependencies = (sims[-1], 'other')
)

Here, the important part is the "dependencies = (sims[-1], 'other')" line which ensures that the last sim
prior to the hp.x calculation (which is always the SCF calculation here) is used as a dependency. 

As a final warning, the new Hubbard input is very specific about the angular momentum channels defined in the
pseudopotential. For the C-diamond example to work, I have modified the C.BFD.upf file converting the "0s"->"2S" and 
"0p"->"2P" in the PP_HEADER and PP_PSWFC sections. You may need to update your pps by hand to make it work. 