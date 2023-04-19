#!/usr/bin/env python

'''
Gamma point Hartree-Fock/DFT

The 2-electron integrals are computed using Poisson solver with FFT by default.
In most scenario, it should be used with pseudo potential.
'''

# Note import path which is different to molecule code
from pyscf.pbc import gto, scf, dft
from pyscf import gto as Mgto
import numpy as np

l_to_char = 'spdfghi'

def bas_str(l):
    return f'''
    O {l.upper()}
     10 0.2
     3 0.8
     1 0.4
     '''
def bas_str_n(n):
    return bas_str(l_to_char[n])

def mkcell(l):
    cell = gto.Cell()
    alat = 6.0
    clat = 20.0
    cell.a = np.array([
        [alat, 0, 0],
        [-0.5*alat, 0.5*alat*np.sqrt(3.0), 0],
        [0, 0, clat]])
    cell.atom = [
            ['He', [0.59*alat, 0.1, 0.5*clat]],
            ['He', [0.35*alat, 0.25*np.sqrt(3.0)*alat, 0.5*clat]]]

    cell.basis = {'He': Mgto.parse(''.join(bas_str(i) for i in l))}
    cell.verbose = 0
    cell.cart = False
    cell.build()
    return cell

def runscf(l):
    cell = mkcell(l)
    mf = scf.RHF(cell).density_fit(auxbasis='weigend')
    mf.conv_tol = 1e-12
    ehf = mf.kernel()
    #import sys
    #sys.path.append(f'{QMCROOT}/src/QMCTools')
    from PyscfToQmcpack import savetoqmcpack
    savetoqmcpack(cell,mf,title=f'he_{l}_sph')
    #from PyscfSphToCart import qmch5
    #qh5 = qmch5(f'he_{l}_sph.h5')
    #qh5.generate_h5(newpath = f'he_{l}_cart.h5')
    return

# add higher l when supported in QMCPACK
#for lbas in 'pdfghi':
for lbas in 'pd':
    # run scf with 'sp', 'sd', 'sf', ... basis
    # need 2+ shells at a time to verify correct relative normalization between shells of different l
    runscf('s'+lbas)
