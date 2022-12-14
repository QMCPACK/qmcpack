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
    return m

s2cdict = {
         2:[((0,1),12),
            ((1,4),12),
            ((2,0),-1),
            ((2,3),-1),
            ((2,5), 4),
            ((3,2),12),
            ((4,0), 3),
            ((4,3),-3)],
         3:[((0, 1), 45),   
            ((0, 6), -5),                
            ((1, 4), 120),  
            ((2, 1), -3), 
            ((2, 6), -3), 
            ((2, 8), 48),   
            ((3, 2), -18),               
            ((3, 7), -18),               
            ((3, 9), 8), 
            ((4, 0), -3),
            ((4, 3), -3),
            ((4, 5), 48),
            ((5, 2), 30),
            ((5, 7), -30),
            ((6, 0), 5),                
            ((6, 3), -45)],
		4:[((0, 1), 140.0),
           ((0, 6), -140.0),
           ((1, 4), 630.0),
           ((1, 11), -70.0),
           ((2, 1), -20.0),
           ((2, 6), -20.0),
           ((2, 8), 720.0),
           ((3, 4), -90.0),
           ((3, 11), -90.0),
           ((3, 13), 160.0),
           ((4, 0), 2.25),
           ((4, 3), 9.0),
           ((4, 5), -144.0),
           ((4, 10), 2.25),
           ((4, 12), -144.0),
           ((4, 14), 16.0),
           ((5, 2), -90.0),
           ((5, 7), -90.0),
           ((5, 9), 160.0),
           ((6, 0), -5.0),
           ((6, 5), 180.0),
           ((6, 10), 5.0),
           ((6, 12), -180.0),
           ((7, 2), 70.0),
           ((7, 7), -630.0),
           ((8, 0), 8.75),
           ((8, 3), -315.0),
           ((8, 10), 8.75)],
		5:[((0, 1), 393.75),
           ((0, 6), -1575.0),
           ((0, 15), 15.75),
           ((1, 4), 2520.0),
           ((1, 11), -2520.0),
           ((2, 1), -78.75),
           ((2, 6), -35.0),
           ((2, 8), 5040.0),
           ((2, 15), 8.75),
           ((2, 17), -560.0),
           ((3, 4), -840.0),
           ((3, 11), -840.0),
           ((3, 13), 3360.0),
           ((4, 1), 7.5),
           ((4, 6), 30.0),
           ((4, 8), -1080.0),
           ((4, 15), 7.5),
           ((4, 17), -1080.0),
           ((4, 19), 480.0),
           ((5, 2), 112.5),
           ((5, 7), 450.0),
           ((5, 9), -800.0),
           ((5, 16), 112.5),
           ((5, 18), -800.0),
           ((5, 20), 32.0),
           ((6, 0), 7.5),
           ((6, 3), 30.0),
           ((6, 5), -1080.0),
           ((6, 10), 7.5),
           ((6, 12), -1080.0),
           ((6, 14), 480.0),
           ((7, 2), -210.0),
           ((7, 9), 840.0),
           ((7, 16), 210.0),
           ((7, 18), -840.0),
           ((8, 0), -8.75),
           ((8, 3), 35.0),
           ((8, 5), 560.0),
           ((8, 10), 78.75),
           ((8, 12), -5040.0),
           ((9, 2), 157.5),
           ((9, 7), -5670.0),
           ((9, 16), 157.5),
           ((10, 0), 15.75),
           ((10, 3), -1575.0),
           ((10, 10), 393.75)],
		6:[((0, 1), 1039.5),
           ((0, 6), -11550.0),
           ((0, 15), 1039.5),
           ((1, 4), 8662.5),
           ((1, 11), -34650.0),
           ((1, 22), 346.5),
           ((2, 1), -252.0),
           ((2, 8), 25200.0),
           ((2, 15), 252.0),
           ((2, 17), -25200.0),
           ((3, 4), -4252.5),
           ((3, 11), -1890.0),
           ((3, 13), 30240.0),
           ((3, 22), 472.5),
           ((3, 24), -3360.0),
           ((4, 1), 52.5),
           ((4, 6), 210.0),
           ((4, 8), -13440.0),
           ((4, 15), 52.5),
           ((4, 17), -13440.0),
           ((4, 19), 13440.0),
           ((5, 4), 525.0),
           ((5, 11), 2100.0),
           ((5, 13), -8400.0),
           ((5, 22), 525.0),
           ((5, 24), -8400.0),
           ((5, 26), 1344.0),
           ((6, 0), -6.25),
           ((6, 3), -56.25),
           ((6, 5), 2025.0),
           ((6, 10), -56.25),
           ((6, 12), 8100.0),
           ((6, 14), -3600.0),
           ((6, 21), -6.25),
           ((6, 23), 2025.0),
           ((6, 25), -3600.0),
           ((6, 27), 64.0),
           ((7, 2), 525.0),
           ((7, 7), 2100.0),
           ((7, 9), -8400.0),
           ((7, 16), 525.0),
           ((7, 18), -8400.0),
           ((7, 20), 1344.0),
           ((8, 0), 13.125),
           ((8, 3), 13.125),
           ((8, 5), -3360.0),
           ((8, 10), -13.125),
           ((8, 14), 3360.0),
           ((8, 21), -13.125),
           ((8, 23), 3360.0),
           ((8, 25), -3360.0),
           ((9, 2), -472.5),
           ((9, 7), 1890.0),
           ((9, 9), 3360.0),
           ((9, 16), 4252.5),
           ((9, 18), -30240.0),
           ((10, 0), -15.75),
           ((10, 3), 393.75),
           ((10, 5), 1575.0),
           ((10, 10), 393.75),
           ((10, 12), -56700.0),
           ((10, 21), -15.75),
           ((10, 23), 1575.0),
           ((11, 2), 346.5),
           ((11, 7), -34650.0),
           ((11, 16), 8662.5),
           ((12, 0), 28.875),
           ((12, 3), -6496.875),
           ((12, 10), 6496.875),
           ((12, 21), -28.875)]}



def transform_block_s2c(l):
    '''
    h5 mo coefs are stored as [mo,ao(s)]
    we want [mo,ao(c)], so 
    '''
    ncart = shsze(l,False)
    nsph  = shsze(l,True)
    lfac = np.sqrt((2*l+1)/(np.pi*2**(l+2)))
    if l<=0:
        return np.eye(nsph,dtype=np.float64)
    else:
        b = np.zeros((nsph,ncart),dtype=np.float64)
        for (i,j),v in s2cdict[l]:
            b[i,j] = v
        b = np.sign(b)*np.sqrt(np.abs(b)*lfac)
        return b @ pyscf_to_gms_order(l)


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
                bset.create_dataset('angular', (1,),'S9',b'cartesian')
                bset.create_dataset('expandYlm', (1,),'S6',b'Gamess')

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

    args = parser.parse_args()

    qh5 = qmch5(args.input_h5path)
    qh5.generate_h5(newpath=args.output)


