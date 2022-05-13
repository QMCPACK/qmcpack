#!/usr/bin/env python3

# This file provided just for reference, not executed in the test

from ase.build import bulk
from ase.units import Ry
from gpaw import GPAW,PW,FermiDirac,setup_paths

silicon = bulk('Si', 'diamond', a=5.459)

econv = 1e-9*2*Ry
ecut  = 30*Ry
silicon.calc = GPAW(xc='PBE',
        setups      = {'default':'sg15'}, # Based on standard library of soft SG15 ecps
        mode        = PW(ecut,force_complex_dtype=True), # give all modes explicitly
        txt         = 'silicon.txt',
        kpts        = (2,2,2),
        occupations = FermiDirac(0.1),
        convergence = {'energy':econv},
        symmetry    = 'off'  
        )
silicon.get_potential_energy()
silicon.calc.write('restart.gpw',mode='all') # important to request mode='all'
