#! /usr/bin/env python3

import numpy
from functools import reduce
from pyscf.pbc import gto, scf 
from pyscf.pbc import tools as pbctools

alat0 = 3.6

cell = gto.Cell()
cell.a = (numpy.ones((3,3))-numpy.eye(3))*alat0/2.0
cell.atom = (('C',0,0,0),('C',numpy.array([0.25,0.25,0.25])*alat0))
cell.basis = 'gth-dzvp'
cell.pseudo = 'gth-pade'
cell.gs = [10]*3  
cell.verbose = 4
cell.build()

mf = scf.RHF(cell)
mf.chkfile = 'scf.gamma.dump'
ehf = mf.kernel()

from pyscf import tools 

c = mf.mo_coeff

h1e = reduce(numpy.dot, (c.T, mf.get_hcore(), c))
eri = mf.with_df.ao2mo(c,compact=True)

madelung = pbctools.pbc.madelung(cell, numpy.zeros(3))
e0 = cell.energy_nuc() + madelung*cell.nelectron * -.5

tools.fcidump.from_integrals('fcidump.gamma.dat', h1e, eri, c.shape[1],
                             cell.nelectron, nuc = e0, ms=0, tol=1e-10)
