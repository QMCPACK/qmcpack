#!/usr/bin/env python3
######################################################################################
## This file is distributed under the University of Illinois/NCSA Open Source License.
## See LICENSE file in top directory for details.
##
## Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
##
## File developed by: Kevin Gasperich, kgasperich@anl.gov, Argonne National Laboratory
##
## File created by:   Kevin Gasperich, kgasperich@anl.gov, Argonne National Laboratory
#######################################################################################

import h5py
import numpy as np
import itertools
import shutil
import functools


@functools.lru_cache()
def doublefactorial(n):
    if n>=2:
        return n * doublefactorial(n-2)
    else:
        return 1

@functools.lru_cache()
def cartnorm(nxyz):
    return np.prod([doublefactorial(2*i-1) for i in nxyz])

def cartnorm_str(s):
    return cartnorm(tuple(s.count(i) for i in 'xyz'))

@functools.lru_cache()
def gms_cart_order(l):
    def sort1(s):
        # return single xyz string in gamess order
        return ''.join(sorted((i*s.count(i) for i in 'xyz'), key=lambda x:(-len(x),x)))

    # unsorted, unordered xyz strings
    sxyz0 = map(''.join, itertools.combinations_with_replacement('xyz', r=l))

    # unsorted, ordered xyz strings
    sxyz1 = map(sort1, sxyz0)

    # sorted, ordered xyz strings
    return sorted(sxyz1, key=lambda s: (sorted((-s.count(i) for i in 'xyz')), s))

def pyscf_cart_order(l):
    return sorted(map(''.join, itertools.combinations_with_replacement('xyz', r=l)))

def pyscf_to_gms_order(l):
    xg = gms_cart_order(l)
    xp = pyscf_cart_order(l)
    m = np.zeros((len(xp),len(xg)))
    #p2g = [(xp.index(''.join(sorted(xi))),i) for i,xi in enumerate(xg)]
    for i,xi in enumerate(xg):
        m[i,xp.index(''.join(sorted(xi)))] = 1
    return m.T

def debug_pyscf_to_gms_order(l):
    xg = gms_cart_order(l)
    xp = pyscf_cart_order(l)
    print(xg)
    print(xp)
    m = np.zeros((len(xp),len(xg)))
    #p2g = [(xp.index(''.join(sorted(xi))),i) for i,xi in enumerate(xg)]
    for i,xi in enumerate(xg):
        m[i,xp.index(''.join(sorted(xi)))] = 1
    print(m)
    return m

# this data was generated using pyscf.gto.cart2sph
# for each ((i,j), v) in one of these lists, the corresponding matrix element of the cart2sph array is
#     sign(v) * sqrt(|v|*(2*l+1)/(4*pi*2^(l)))
#
# equivalent to:
#     import pyscf
#     import numpy as np
#     
#     def s2c_scaled_nonzero(l,tol=1E-12):
#         coefs = pyscf.gto.cart2sph(l)
#         scaled_coefs = coefs**2 *((4*np.pi)*2**l / (2*l+1)) * np.sign(coefs)
#         return [((i,j),cij.round(8)) 
#                  for i,row in enumerate(scaled_coefs.T) 
#                    for j,cij in enumerate(row) 
#                      if abs(cij)>tol]
#     
#     s2cdict = {l:s2c_scaled_nonzero(l) for l in range(2,7)}

s2cdict = {
        2: [((0, 1), 12.0), ((1, 4), 12.0), ((2, 0), -1.0), ((2, 3), -1.0), ((2, 5), 4.0), ((3, 2), 12.0), ((4, 0), 3.0), ((4, 3), -3.0)],
        3: [((0, 1), 45.0), ((0, 6), -5.0), ((1, 4), 120.0), ((2, 1), -3.0), ((2, 6), -3.0), ((2, 8), 48.0), ((3, 2), -18.0), ((3, 7), -18.0), ((3, 9), 8.0), ((4, 0), -3.0), ((4, 3), -3.0), ((4, 5), 48.0), ((5, 2), 30.0), ((5, 7), -30.0), ((6, 0), 5.0), ((6, 3), -45.0)],
        4: [((0, 1), 140.0), ((0, 6), -140.0), ((1, 4), 630.0), ((1, 11), -70.0), ((2, 1), -20.0), ((2, 6), -20.0), ((2, 8), 720.0), ((3, 4), -90.0), ((3, 11), -90.0), ((3, 13), 160.0), ((4, 0), 2.25), ((4, 3), 9.0), ((4, 5), -144.0), ((4, 10), 2.25), ((4, 12), -144.0), ((4, 14), 16.0), ((5, 2), -90.0), ((5, 7), -90.0), ((5, 9), 160.0), ((6, 0), -5.0), ((6, 5), 180.0), ((6, 10), 5.0), ((6, 12), -180.0), ((7, 2), 70.0), ((7, 7), -630.0), ((8, 0), 8.75), ((8, 3), -315.0), ((8, 10), 8.75)],
        5: [((0, 1), 393.75), ((0, 6), -1575.0), ((0, 15), 15.75), ((1, 4), 2520.0), ((1, 11), -2520.0), ((2, 1), -78.75), ((2, 6), -35.0), ((2, 8), 5040.0), ((2, 15), 8.75), ((2, 17), -560.0), ((3, 4), -840.0), ((3, 11), -840.0), ((3, 13), 3360.0), ((4, 1), 7.5), ((4, 6), 30.0), ((4, 8), -1080.0), ((4, 15), 7.5), ((4, 17), -1080.0), ((4, 19), 480.0), ((5, 2), 112.5), ((5, 7), 450.0), ((5, 9), -800.0), ((5, 16), 112.5), ((5, 18), -800.0), ((5, 20), 32.0), ((6, 0), 7.5), ((6, 3), 30.0), ((6, 5), -1080.0), ((6, 10), 7.5), ((6, 12), -1080.0), ((6, 14), 480.0), ((7, 2), -210.0), ((7, 9), 840.0), ((7, 16), 210.0), ((7, 18), -840.0), ((8, 0), -8.75), ((8, 3), 35.0), ((8, 5), 560.0), ((8, 10), 78.75), ((8, 12), -5040.0), ((9, 2), 157.5), ((9, 7), -5670.0), ((9, 16), 157.5), ((10, 0), 15.75), ((10, 3), -1575.0), ((10, 10), 393.75)],
        6: [((0, 1), 1039.5), ((0, 6), -11550.0), ((0, 15), 1039.5), ((1, 4), 8662.5), ((1, 11), -34650.0), ((1, 22), 346.5), ((2, 1), -252.0), ((2, 8), 25200.0), ((2, 15), 252.0), ((2, 17), -25200.0), ((3, 4), -4252.5), ((3, 11), -1890.0), ((3, 13), 30240.0), ((3, 22), 472.5), ((3, 24), -3360.0), ((4, 1), 52.5), ((4, 6), 210.0), ((4, 8), -13440.0), ((4, 15), 52.5), ((4, 17), -13440.0), ((4, 19), 13440.0), ((5, 4), 525.0), ((5, 11), 2100.0), ((5, 13), -8400.0), ((5, 22), 525.0), ((5, 24), -8400.0), ((5, 26), 1344.0), ((6, 0), -6.25), ((6, 3), -56.25), ((6, 5), 2025.0), ((6, 10), -56.25), ((6, 12), 8100.0), ((6, 14), -3600.0), ((6, 21), -6.25), ((6, 23), 2025.0), ((6, 25), -3600.0), ((6, 27), 64.0), ((7, 2), 525.0), ((7, 7), 2100.0), ((7, 9), -8400.0), ((7, 16), 525.0), ((7, 18), -8400.0), ((7, 20), 1344.0), ((8, 0), 13.125), ((8, 3), 13.125), ((8, 5), -3360.0), ((8, 10), -13.125), ((8, 14), 3360.0), ((8, 21), -13.125), ((8, 23), 3360.0), ((8, 25), -3360.0), ((9, 2), -472.5), ((9, 7), 1890.0), ((9, 9), 3360.0), ((9, 16), 4252.5), ((9, 18), -30240.0), ((10, 0), -15.75), ((10, 3), 393.75), ((10, 5), 1575.0), ((10, 10), 393.75), ((10, 12), -56700.0), ((10, 21), -15.75), ((10, 23), 1575.0), ((11, 2), 346.5), ((11, 7), -34650.0), ((11, 16), 8662.5), ((12, 0), 28.875), ((12, 3), -6496.875), ((12, 10), 6496.875), ((12, 21), -28.875)],
        7: [((0, 1), 2627.625), ((0, 6), -65690.625), ((0, 15), 23648.625), ((0, 28), -53.625), ((1, 4), 27027.0), ((1, 11), -300300.0), ((1, 22), 27027.0), ((2, 1), -721.875), ((2, 6), 721.875), ((2, 8), 103950.0), ((2, 15), 2338.875), ((2, 17), -415800.0), ((2, 28), -28.875), ((2, 30), 4158.0), ((3, 4), -16632.0), ((3, 13), 184800.0), ((3, 22), 16632.0), ((3, 24), -184800.0), ((4, 1), 212.625), ((4, 6), 590.625), ((4, 8), -85050.0), ((4, 15), 23.625), ((4, 17), -37800.0), ((4, 19), 151200.0), ((4, 28), -23.625), ((4, 30), 9450.0), ((4, 32), -16800.0), ((5, 4), 4725.0), ((5, 11), 18900.0), ((5, 13), -134400.0), ((5, 22), 4725.0), ((5, 24), -134400.0), ((5, 26), 48384.0), ((6, 1), -21.875), ((6, 6), -196.875), ((6, 8), 12600.0), ((6, 15), -196.875), ((6, 17), 50400.0), ((6, 19), -50400.0), ((6, 28), -21.875), ((6, 30), 12600.0), ((6, 32), -50400.0), ((6, 34), 3584.0), ((7, 2), -612.5), ((7, 7), -5512.5), ((7, 9), 22050.0), ((7, 16), -5512.5), ((7, 18), 88200.0), ((7, 20), -14112.0), ((7, 29), -612.5), ((7, 31), 22050.0), ((7, 33), -14112.0), ((7, 35), 128.0), ((8, 0), -21.875), ((8, 3), -196.875), ((8, 5), 12600.0), ((8, 10), -196.875), ((8, 12), 50400.0), ((8, 14), -50400.0), ((8, 21), -21.875), ((8, 23), 12600.0), ((8, 25), -50400.0), ((8, 27), 3584.0), ((9, 2), 1181.25), ((9, 7), 1181.25), ((9, 9), -33600.0), ((9, 16), -1181.25), ((9, 20), 12096.0), ((9, 29), -1181.25), ((9, 31), 33600.0), ((9, 33), -12096.0), ((10, 0), 23.625), ((10, 3), -23.625), ((10, 5), -9450.0), ((10, 10), -590.625), ((10, 12), 37800.0), ((10, 14), 16800.0), ((10, 21), -212.625), ((10, 23), 85050.0), ((10, 25), -151200.0), ((11, 2), -1039.5), ((11, 7), 25987.5), ((11, 9), 11550.0), ((11, 16), 25987.5), ((11, 18), -415800.0), ((11, 29), -1039.5), ((11, 31), 11550.0), ((12, 0), -28.875), ((12, 3), 2338.875), ((12, 5), 4158.0), ((12, 10), 721.875), ((12, 12), -415800.0), ((12, 21), -721.875), ((12, 23), 103950.0), ((13, 2), 750.75), ((13, 7), -168918.75), ((13, 16), 168918.75), ((13, 29), -750.75), ((14, 0), 53.625), ((14, 3), -23648.625), ((14, 10), 65690.625), ((14, 21), -2627.625)],
        8: [((0, 1), 6435.0), ((0, 6), -315315.0), ((0, 15), 315315.0), ((0, 28), -6435.0), ((1, 4), 78828.75), ((1, 11), -1970718.75), ((1, 22), 709458.75), ((1, 37), -1608.75), ((2, 1), -1930.5), ((2, 6), 10510.5), ((2, 8), 378378.0), ((2, 15), 10510.5), ((2, 17), -4204200.0), ((2, 28), -1930.5), ((2, 30), 378378.0), ((3, 4), -56306.25), ((3, 11), 56306.25), ((3, 13), 900900.0), ((3, 22), 182432.25), ((3, 24), -3603600.0), ((3, 37), -2252.25), ((3, 39), 36036.0), ((4, 1), 693.0), ((4, 6), 693.0), ((4, 8), -399168.0), ((4, 15), -693.0), ((4, 19), 1108800.0), ((4, 28), -693.0), ((4, 30), 399168.0), ((4, 32), -1108800.0), ((5, 4), 23388.75), ((5, 11), 64968.75), ((5, 13), -1039500.0), ((5, 22), 2598.75), ((5, 24), -462000.0), ((5, 26), 665280.0), ((5, 37), -2598.75), ((5, 39), 115500.0), ((5, 41), -73920.0), ((6, 1), -157.5), ((6, 6), -1417.5), ((6, 8), 141750.0), ((6, 15), -1417.5), ((6, 17), 567000.0), ((6, 19), -1008000.0), ((6, 28), -157.5), ((6, 30), 141750.0), ((6, 32), -1008000.0), ((6, 34), 161280.0), ((7, 4), -2756.25), ((7, 11), -24806.25), ((7, 13), 176400.0), ((7, 22), -24806.25), ((7, 24), 705600.0), ((7, 26), -254016.0), ((7, 37), -2756.25), ((7, 39), 176400.0), ((7, 41), -254016.0), ((7, 43), 9216.0), ((8, 0), 19.140625), ((8, 3), 306.25), ((8, 5), -19600.0), ((8, 10), 689.0625), ((8, 12), -176400.0), ((8, 14), 176400.0), ((8, 21), 306.25), ((8, 23), -176400.0), ((8, 25), 705600.0), ((8, 27), -50176.0), ((8, 36), 19.140625), ((8, 38), -19600.0), ((8, 40), 176400.0), ((8, 42), -50176.0), ((8, 44), 256.0), ((9, 2), -2756.25), ((9, 7), -24806.25), ((9, 9), 176400.0), ((9, 16), -24806.25), ((9, 18), 705600.0), ((9, 20), -254016.0), ((9, 29), -2756.25), ((9, 31), 176400.0), ((9, 33), -254016.0), ((9, 35), 9216.0), ((10, 0), -39.375), ((10, 3), -157.5), ((10, 5), 35437.5), ((10, 12), 35437.5), ((10, 14), -252000.0), ((10, 21), 157.5), ((10, 23), -35437.5), ((10, 27), 40320.0), ((10, 36), 39.375), ((10, 38), -35437.5), ((10, 40), 252000.0), ((10, 42), -40320.0), ((11, 2), 2598.75), ((11, 7), -2598.75), ((11, 9), -115500.0), ((11, 16), -64968.75), ((11, 18), 462000.0), ((11, 20), 73920.0), ((11, 29), -23388.75), ((11, 31), 1039500.0), ((11, 33), -665280.0), ((12, 0), 43.3125), ((12, 3), -693.0), ((12, 5), -24948.0), ((12, 10), -4331.25), ((12, 12), 623700.0), ((12, 14), 69300.0), ((12, 21), -693.0), ((12, 23), 623700.0), ((12, 25), -2494800.0), ((12, 36), 43.3125), ((12, 38), -24948.0), ((12, 40), 69300.0), ((13, 2), -2252.25), ((13, 7), 182432.25), ((13, 9), 36036.0), ((13, 16), 56306.25), ((13, 18), -3603600.0), ((13, 29), -56306.25), ((13, 31), 900900.0), ((14, 0), -53.625), ((14, 3), 10510.5), ((14, 5), 10510.5), ((14, 12), -2364862.5), ((14, 21), -10510.5), ((14, 23), 2364862.5), ((14, 36), 53.625), ((14, 38), -10510.5), ((15, 2), 1608.75), ((15, 7), -709458.75), ((15, 16), 1970718.75), ((15, 29), -78828.75), ((16, 0), 100.546875), ((16, 3), -78828.75), ((16, 10), 492679.6875), ((16, 21), -78828.75), ((16, 36), 100.546875)]}


def lm_cart_rescale(l):
    '''
    rescaling of certain functions within shells
    these arise due to different normalization factors for the different Cartesian functions
    a Cartesian GTO:
       x^(n_x) y^(n_y) z^(n_z) exp(-a r^2)
    is relatively normalized (i.e. relative to other function with the same shell)
    by a factor of:
       product((2*n_k-1)!! for n_k in (n_x, n_y, n_z))
    '''
    if l<=1:
        return np.eye(shsze(l,False))
    else:
        n0 = cartnorm_str('x'*l)
        return np.diag(np.sqrt([cartnorm_str(si)/n0 for si in gms_cart_order(l)])) * (2.0)**(-0.5*l)

def transform_block_s2c(l):
    '''
    h5 mo coefs are stored as [mo,ao(s)]
    we want [mo,ao(c)], so multiply by [ao(s),ao(c)]
    '''
    ncart = shsze(l,False)
    nsph  = shsze(l,True)
    if l<=1:
        return np.eye(nsph,dtype=np.float64)
    else:
        b = np.zeros((nsph,ncart),dtype=np.float64)
        for (i,j),v in s2cdict[l]:
            b[i,j] = v
        b = np.sign(b)*np.sqrt(np.abs(b))
        return (b @ pyscf_to_gms_order(l) @ lm_cart_rescale(l))


def shsze(l,is_sph):
    '''
    number of AOs in shell with angular momentum l
    '''
    return 2*l+1 if is_sph else ((l+1)*(l+2))//2

def atom_shell_table(atom_l_sph):
    '''
    input is list of (l, is_sph) (one for each shell)
    output is list of [l, i0_cart, i1_cart, i0_sph, i1_sph] (one for each shell)
        ao indices (start/end+1) for each shell
    '''
    ntot_sph=0
    ntot_cart=0
    dat = []
    for i,(il,isph) in enumerate(atom_l_sph):
        nao_sph = shsze(il,True)
        nao_cart = shsze(il,False)
        dat.append((il, ntot_cart, ntot_cart+nao_cart, ntot_sph, ntot_sph+nao_sph))
        ntot_cart += nao_cart
        ntot_sph  += nao_sph
    return dat

def get_h5f_dsets(h5f):
    '''
    return list of paths to all dataset objects in an hdf5 file
    '''
    sets = []
    find_dsets(h5f, sets=sets)
    return sets


def find_dsets(h5f, sets=[], pwd=''):
    '''
    h5f: h5py file object
    sets: list for accumulating h5paths to datasets
    pwd: current location in h5 file

    recurse through hdf5 file and generate list of paths to all dataset objects
    '''
    for key in h5f.keys():
        if type(h5f[key]) == h5py._hl.dataset.Dataset:
            sets.append(pwd+'/'+key)
        else:
            find_dsets(h5f[key], sets, pwd+'/'+key)
    return

class qmch5:
    '''
    class for converting qmcpack h5 file from spherical to cartesian AO basis
        modified data includes:
        /basisset/atomicBasisSet[iat]/angular:   "spherical" -> "cartesian"
        /basisset/atomicBasisSet[iat]/expandYlm: "pyscf" -> "Gamess"
        /*/eigenset*: transform d,f,g,etc. sub-blocks
        /parameters/numAO


        relevant data:
        /atoms/number_of_atoms
        /atoms/number_of_species
        /atoms/species_ids
        /basisset/atomicBasisSet[iat]/NbBasisGroups
        /basisset/atomicBasisSet[iat]/basisGroup[igrp]/l
    '''
    def __init__(self,h5path):
        self.h5path = h5path
        with h5py.File(h5path,'r') as f:
            self.dset_paths = get_h5f_dsets(f)
        self.set_params()
        return

    def set_params(self):
        '''
        set data members of object
        n_atoms
        n_species
        species_ids:  list of ids, one for each atom
        n_bas_el:     number of atomic basis sets
        n_bas_groups: number of shells(?) in each atomic basis
        bas_l_sph:    2-level list of 2-tuples over [atoms, [shells]]
                      contains (l, is_spherical) for each shell
        '''
        with h5py.File(self.h5path,'r') as f:
            self.n_atoms = f['atoms/number_of_atoms'][0]
            self.n_species = f['atoms/number_of_species'][0]
            self.species_ids = f['atoms/species_ids'][()]
            assert(len(self.species_ids) == self.n_atoms)
            assert(max(self.species_ids) < self.n_species)
            self.n_bas_el = f['basisset/NbElements'][0]
            self.n_bas_groups = [self.ao_bs(iat,f)['NbBasisGroups'][0] for iat in range(self.n_bas_el)]
            # nested list of 2-tuples (l, is_spherical); dims are atoms, shells
            self.bas_l_sph = [
                    [(self.ao_grp(iat,igrp,f)['l'][0],
                        self.ao_bs(iat,f)['angular'][0].decode()=="spherical") 
                        for igrp in range(ngrp)] 
                    for iat,ngrp in enumerate(self.n_bas_groups)]
        self.set_shell_ao_table()
        self.set_ao_idx_table()
        return

    def set_shell_ao_table(self):
        '''
        3-level list over [atom, [shell, [data]]] where data is: 
        [l, i0c, i1c, i0s, i1s]
        where i0/i1 is index of first/last+1 AO in cart/sph basis
        '''
        self.shell_ao_table = [atom_shell_table(i) for i in self.bas_l_sph]
        return

    def set_ao_idx_table(self):
        '''
        similar to shell_ao_table but flattened/accumulated over all atoms in system
        list of [l, i0c, i1c, i0s, i1s] for each shell
        useful to identify nonzero sub-blocks of sph-to-cart transformation matrix
        '''
        ntot_c = 0
        ntot_s = 0
        self.ao_idx_table = []
        for ispec in self.species_ids:
            for l, _ in self.bas_l_sph[ispec]:
                nao_c = shsze(l,False)
                nao_s = shsze(l,True)
                self.ao_idx_table.append([l,
                    ntot_c, ntot_c + nao_c,
                    ntot_s, ntot_s + nao_s])
                ntot_c += nao_c
                ntot_s += nao_s
        return


    # utility functions to get basisset and group objects
    def ao_bs(self,i,h5f):
        return h5f[f'basisset/atomicBasisSet{i}']
    def ao_grp(self, iat, igrp, h5f):
        return self.ao_bs(iat,h5f)[f'basisGroup{igrp}']


    def s2c_transformation_matrix(self):
        '''
        block-diagonal (rectangular) transformation matrix for sph to cart
        '''
        ntot_c = self.ao_idx_table[-1][2]
        ntot_s = self.ao_idx_table[-1][4]
        A = np.zeros((ntot_s,ntot_c),dtype=np.float64)
        for l,i0c,i1c,i0s,i1s in self.ao_idx_table:
            A[i0s:i1s,i0c:i1c] = transform_block_s2c(l)
        return A

    def get_eigenset_paths(self):
        '''
        location of each object in h5 file with 'eigenset' somewhere in the path
        '''
        return [i for i in self.dset_paths if 'eigenset' in i]

    def generate_h5(self, newpath=None):
        '''
        make copy of h5 file and do s2c conversion
        move from 'name.h5' to 'name.s2c.h5'
        '''
        # copy file
        if not newpath:
            if self.h5path.endswith('.h5'):
                newpath = self.h5path.replace('.h5','.s2c.h5')
            else:
                newpath = self.h5path+'.s2c.h5'
        self.newpath = newpath
        print(f'\tperforming sph-to-cart conversion on {self.h5path}; saving as {self.newpath}')
        shutil.copy(self.h5path,self.newpath)
        with h5py.File(self.newpath,'r+') as h5f:
            # change basis info
            for ibset in range(self.n_bas_el):
                bset = self.ao_bs(ibset,h5f)
                assert(bset['angular'][()] == b'spherical')
                assert(bset['expandYlm'][()] == b'pyscf')
                del bset['angular']
                del bset['expandYlm']
                bset.create_dataset('angular', (1,),'S9',np.bytes_("cartesian"))
                bset.create_dataset('expandYlm', (1,),'S6',np.bytes_("Gamess"))

            # generate transformation matrix (over all AOs)
            A = self.s2c_transformation_matrix()

            # transform anything with "eigenset" in its name
            for eigenset in self.get_eigenset_paths():
                print("\t\ttransforming eigenset:          ", h5f[eigenset])
                coef_sph = h5f[eigenset][()]
                coef_cart = coef_sph @ A
                del h5f[eigenset]
                h5f.create_dataset(eigenset, data=coef_cart)
                print("\t\tfinished transforming eigenset: ", h5f[eigenset])

            # set numAO to correct value
            print("\t\told numAO: ", h5f['parameters/numAO'][()])
            assert(h5f['parameters/numAO'][()] == A.shape[0])
            h5f['parameters/numAO'][()] = A.shape[1]
            print("\t\tnew numAO: ", h5f['parameters/numAO'][()])
        return

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser('convert MO coefs from spherical to cartesian')

    parser.add_argument('input_h5path', type=str, help='h5 file to convert')
    parser.add_argument('-o','--output', type=str, help='output h5 file')
    #parser.add_argument('--debug', action='store_true', help='debug')

    args = parser.parse_args()

    qh5 = qmch5(args.input_h5path)
    qh5.generate_h5(newpath=args.output)
