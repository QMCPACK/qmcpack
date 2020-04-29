#! /usr/bin/env python3
# Useful routines

import h5py
import numpy
from pyscf import lib
from pyscf.lib.chkfile import load, load_mol
from pyscf.pbc.lib.chkfile import load_cell
from afqmctools.utils.linalg import get_ortho_ao_mol

def load_from_pyscf_chk(chkfile,hcore=None,orthoAO=False):

    cell = load_cell(chkfile)
    assert(cell is not None)

    kpts=None
    singleK = False
    if lib.chkfile.load(chkfile, 'scf/kpt') is not None :
        kpts = numpy.asarray(lib.chkfile.load(chkfile, 'scf/kpt'))
        singleK = True
    else:
        kpts = numpy.asarray(lib.chkfile.load(chkfile, 'scf/kpts'))
        assert(kpts is not None)
    kpts = numpy.reshape(kpts,(-1,3))
    nkpts = len(kpts)
    nao = cell.nao_nr()
    nao_tot = nao*nkpts

    Xocc = lib.chkfile.load(chkfile, 'scf/mo_occ')
    mo_energy = lib.chkfile.load(chkfile, 'scf/mo_energy')
    mo_coeff = lib.chkfile.load(chkfile, 'scf/mo_coeff')
    fock = numpy.asarray(lib.chkfile.load(chkfile, 'scf/fock'))
    assert(fock is not None)
    if isinstance(Xocc,list):
        # 3 choices:
        if isinstance(Xocc[0],list):
            # KUHF
            isUHF = True
            assert(len(Xocc[0])==nkpts)
        elif singleK:
            # UHF
            isUHF = True
            assert(len(Xocc) == 2)
            Xocc = numpy.asarray(Xocc)
        else:
            # KRHF
            isUHF = False
            assert(len(Xocc) == nkpts)
    else:
        assert(singleK)
        if len(Xocc) == 2:
            isUHF = True
        else:
            # single kpoint RHF
            isUHF = False
            Xocc = ([Xocc])
            assert(len(Xocc)==nkpts)

    if hcore is None:
        hcore = numpy.asarray(lib.chkfile.load(chkfile, 'scf/hcore'))
        assert(hcore is not None)
    hcore = numpy.reshape(hcore,(-1,nao,nao))
    assert(hcore.shape[0]==nkpts)

    if(cell.spin!=0 and isUHF==False):
        print(" cell.spin!=0 only allowed with UHF calculation \n")
        comm.abort()

    if(orthoAO==False and isUHF==True):
        print(" orthoAO=True required with UHF calculation \n")
        quit()

    if orthoAO:
        X_ = numpy.asarray(lib.chkfile.load(chkfile, 'scf/orthoAORot')).reshape(nkpts,nao,-1)
        assert(X_ is not None)
        nmo_pk = numpy.asarray(lib.chkfile.load(chkfile, 'scf/nmo_per_kpt'))
        # do this properly!!!
        if len(nmo_pk.shape) == 0:
            nmo_pk = numpy.asarray([nmo_pk])
        X = []
        for k in range(len(nmo_pk)):
            X.append(X_[k][:,0:nmo_pk[k]])
        assert(nmo_pk is not None)
    else:
        # can safely assume isUHF == False
        X = lib.chkfile.load(chkfile, 'scf/mo_coeff')
        if singleK:
            assert(len(X.shape) == 2)
            assert(X.shape[0] == nao)
            X = ([X])
        assert(len(X) == nkpts)
        nmo_pk = numpy.zeros(nkpts,dtype=numpy.int32)
        for ki in range(nkpts):
            nmo_pk[ki]=X[ki].shape[1]
            assert(nmo_pk[ki] == Xocc[ki].shape[0])

    if singleK:
        assert(nkpts==1)
        if isUHF:
            assert len(fock.shape) == 3
            assert fock.shape[0] == 2
            assert fock.shape[1] == nao
            fock = fock.reshape((2,1,fock.shape[1],fock.shape[2]))
            assert len(Xocc.shape) == 2
            Xocc = Xocc.reshape((2,1,Xocc.shape[1]))
            assert len(mo_energy.shape) == 2
            mo_energy = mo_energy.reshape((2,1,mo_energy.shape[1]))
        else:
            assert len(fock.shape) == 2
            assert fock.shape[0] == nao
            fock = fock.reshape((1,1)+fock.shape)
            mo_energy = mo_energy.reshape((1,-1))
    if len(fock.shape) == 3:
        fock = fock.reshape((1,)+fock.shape)
    scf_data = {'cell': cell, 'kpts': kpts,
                'Xocc': Xocc, 'isUHF': isUHF,
                'hcore': hcore, 'X': X, 'nmo_pk': nmo_pk,
                'mo_coeff': mo_coeff,
                'nao': nao, 'fock': fock,
                'mo_energy': mo_energy}
    return scf_data

def load_from_pyscf_chk_mol(chkfile, base='scf'):
    mol = load_mol(chkfile)
    with h5py.File(chkfile, 'r') as fh5:
        try:
            hcore = fh5['/scf/hcore'][:]
        except KeyError:
            hcore = mol.intor_symmetric('int1e_nuc')
            hcore += mol.intor_symmetric('int1e_kin')
            if len(mol._ecpbas) > 0:
                hcore += mol.intor_symmetric('ECPScalar')
        try:
            X = fh5['/scf/orthoAORot'][:]
        except KeyError:
            s1e = mol.intor('int1e_ovlp_sph')
            X = get_ortho_ao_mol(s1e)
        try:
            df_ints = fh5['j3c'][:]
        except KeyError:
            df_ints = None
    mo_occ = numpy.array(lib.chkfile.load(chkfile, base+'/mo_occ'))
    mo_coeff = numpy.array(lib.chkfile.load(chkfile, base+'/mo_coeff'))
    uhf = len(mo_coeff.shape) == 3
    if mol.nelec[0] != mol.nelec[1] and not uhf:
        rohf = True
    else:
        rohf = False
    scf_data = {'mol': mol, 'mo_occ': mo_occ, 'hcore': hcore,
                'X': X, 'mo_coeff': mo_coeff,
                'isUHF': uhf, 'df_ints': df_ints,
                'rohf': rohf}
    return scf_data
