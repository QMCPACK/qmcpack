from h5py import File
import numpy
from pyscf.pbc import gto, tools
from pyscf.pbc.dft import numint
from pyscf import gto as molgto
import os
import numpy
from mpi4py import MPI
from gto_basis_utils import extend_gto_id
from pyscf_driver import (pyscf_driver_init, pyscf_driver_get_info, pyscf_driver_end, 
                pyscf_driver_mp2,pyscf_driver_hamil,pyscf_driver_mp2no)

def make_cell(latt,sp_label,atid,atpos,basis_='gthdzvp',pseudo_='gthpbe',mesh=None,prec=1e-8):
    assert(len(atid) == len(atpos))
    assert(atpos.ndim == 2)
    assert(atpos.shape[1] == 3)
    cell = gto.Cell()
    cell.a = '''
       {} {} {}
       {} {} {}
       {} {} {}'''.format(latt[0,0],latt[0,1],latt[0,2],
                          latt[1,0],latt[1,1],latt[1,2],
                          latt[2,0],latt[2,1],latt[2,2])
    atom = ''
    for i in range(len(atid)):
        atom += '''{} {} {} {} \n'''.format(sp_label[atid[i]-1],
                                            atpos[i,0],atpos[i,1],atpos[i,2])
    cell.atom = atom
    cell.basis = basis_
    cell.pseudo = pseudo_
    if mesh is not None:
        cell.mesh = mesh
    cell.verbose = 5
    cell.unit = 'B'
    cell.precision = prec
    cell.build()
    return cell

def write_esh5_orbitals(cell, name):
    def to_qmcpack_complex(array):
        shape = array.shape
        return array.view(numpy.float64).reshape(shape+(2,))
    nao = cell.nao_nr()
    kpts = numpy.zeros((1,3),dtype=numpy.float64)

    fh5 = File(name,'w')
    coords = cell.gen_uniform_grids(cell.mesh)

    norbs = numpy.zeros((1,),dtype=int)
    norbs[0] = nao

    norbs = numpy.zeros((1,),dtype=int)
    norbs[0] = nao
    grp = fh5.create_group("OrbsG")
    dset = grp.create_dataset("reciprocal_vectors", data=cell.reciprocal_vectors())
    dset = grp.create_dataset("number_of_kpoints", data=len(kpts))
    dset = grp.create_dataset("kpoints", data=kpts)
    dset = grp.create_dataset("number_of_orbitals", data=norbs)
    dset = grp.create_dataset("fft_grid", data=cell.mesh)
    dset = grp.create_dataset("grid_type", data=int(0))
    nnr = cell.mesh[0]*cell.mesh[1]*cell.mesh[2]
    orb = numpy.zeros((nnr,2),dtype=numpy.float64)
    aoi_G = numpy.empty((1,len(coords)), dtype=numpy.complex128)
    # loop over kpoints later
    for (ik,k) in enumerate(kpts):
        ao = numint.KNumInt().eval_ao(cell, coords, k)[0]
        fac = numpy.exp(-1j * numpy.dot(coords, k))
        for i in range(norbs[ik]):
            aoi = fac * numpy.asarray(ao[:,i].T, order='C')
            aoi_G = tools.fft(aoi, cell.mesh)
            aoi_G = aoi_G.reshape(cell.mesh).transpose(2,1,0).reshape(nnr)
            dset = grp.create_dataset('kp'+str(ik)+'_b'+str(i), data=to_qmcpack_complex(aoi_G))
    fh5.close()

def make_image_comm(nimage, comm=MPI.COMM_WORLD):
    parent_nproc = comm.size
    parent_mype = comm.rank
    assert( parent_nproc%nimage == 0 )
    nproc_image = parent_nproc / nimage
    my_image_id = parent_mype / nproc_image
    me_image    = parent_mype%nproc_image 
    intra_image = comm.Split(my_image_id,comm.rank)
    inter_image = comm.Split(me_image,comm.rank)
    return intra_image, inter_image

# put these in modules later
def qe_driver_init(norb, qe_prefix, qe_outdir, atm_labels,
                   intra_image=MPI.COMM_WORLD, inter_image=None,
                   npools=1, outdir='./qedrv', #remove_dir=True, 
                   set_soft_links=True, verbose=True,add_image_tag=True):
#    if remove_dir:
#        assert(set_soft_links)
    assert(intra_image.size%npools==0)
    intra_rank = intra_image.rank
    intra_size = intra_image.size

    if inter_image is not None:
        inter_rank = inter_image.rank
        inter_size = inter_image.size
    else:
        inter_rank = 0
        inter_size = 1

    fname = outdir
    if add_image_tag: 
        fname += '.'+str(inter_rank)+'/'
    else:
        fname += '/'
    if  intra_rank == 0:
#        if remove_dir:
#            os.system('rm -rf '+fname+'\n')
        if set_soft_links:
            os.system('mkdir '+fname)
            os.system('ln -s ./'+qe_outdir+'/'+qe_prefix+'.xml '+fname+'/'+qe_prefix+'.xml')
            os.system('ln -s ./'+qe_outdir+'/'+qe_prefix+'.save/ '+fname+'/'+qe_prefix+'.save')
    MPI.COMM_WORLD.barrier()

    # initialize driver: 
    nkpts, nat, nsp, npwx, ngm, mesh = pyscf_driver_init(
            inter_size,npools,intra_size/npools,norb,qe_prefix,fname,verbose)
    atms = numpy.array(atm_labels)  # don't know how to return an array of strings
    atom_ids,atom_pos,kpts,latt = pyscf_driver_get_info(nat,nsp,nkpts)#,atms)
    atom_pos=atom_pos.T
    kpts=kpts.T
    latt=latt.T
    qe_info = {'species' : atms,
               'nsp' : nsp,
               'nat' : nat,
               'at_id' : atom_ids, 
               'at_pos' : atom_pos, 
               'nkpts' : nkpts, 
               'kpts' : kpts, 
               'latt' : latt, 
               'npwx' : npwx, 
               'mesh' : mesh, 
               'ngm' : ngm, 
               'outdir' : fname
               } 
    if verbose and (MPI.COMM_WORLD.rank==0):
        print("# species = {}".format(qe_info['nsp']))
        print("# atoms = {}".format(qe_info['nat']))
        print("# kpts = {}".format(qe_info['nkpts']))
        print("FFT mesh = {} {} {}".format(qe_info['mesh'][0],qe_info['mesh'][1],qe_info['mesh'][2]))
        print(" Atom species: ")
        print(qe_info['species'])
        print(" Atom positions: ")
        print(qe_info['at_pos'])
        print(" Lattice: ")
        print(qe_info['latt'])
        print(" K-points: ")
        print(qe_info['kpts'])

    return qe_info 

def qe_driver_end():
    if(MPI.COMM_WORLD.rank==0):
        print(" Closing QE driver.")
    pyscf_driver_end()

def gen_qe_gto(qe_info,gto_files={},x=None,basis_map={},  
               fname='pyscf.orbitals.h5',prec=1e-12):
    basis = {}
    n0=0
    n1=0
    for I, atm in enumerate(qe_info['species']):
        gtof = ''
        if atm in gto_files:
            gtof = gto_files[atm]
        if atm in basis_map:
            n1 += len(basis_map[atm])
            assert(x is not None)
            assert(len(x) >= n1)
            basis.update({atm: molgto.parse(extend_gto_id(gtof,atm,
                        x[n0:n1],basis_map[atm]))})
            n0 += len(basis_map[atm])
        else:
            assert(os.path.exists(gtof))
            basis.update({atm: molgto.parse(open(gtof,"r").read())})

    cell = make_cell(qe_info['latt'],
                 qe_info['species'],
                 qe_info['at_id'],
                 qe_info['at_pos'],
                 basis_=basis,
                 mesh=qe_info['mesh'],prec=prec)
    nao = cell.nao_nr()
    write_esh5_orbitals(cell, qe_info['outdir']+fname)
    return nao

def qe_driver_MP2(qe_info,out_prefix='pyscf_drv',
                        diag_type='keep_occ',
                        nread_from_h5=0,h5_add_orbs='',
                        eigcut=1e-3,nextracut=1e-6,kappa=0.0,regp=0):
    if diag_type=='fullpw':
        emp2=pyscf_driver_mp2(out_prefix,True,diag_type,
                     0,'',0.0,
                     0.0,kappa,regp)
    else:
        emp2=pyscf_driver_mp2(out_prefix,True,diag_type,
                     nread_from_h5,h5_add_orbs,eigcut,
                     nextracut,kappa,regp)
    return emp2

def qe_driver_MP2NO(qe_info,out_prefix='pyscf_drv',
                        appnos=False,
                        diag_type='keep_occ',
                        nread_from_h5=0,h5_add_orbs='',nskip=0,
                        eigcut=1e-3,nextracut=1e-6,mp2noecut=1e-6,kappa=0.0,regp=0):
    if diag_type=='fullpw':
        pyscf_driver_mp2no(out_prefix,True,diag_type,appnos,
                     0,'',nskip,0.0,
                     0.0,mp2noecut,kappa,regp)
    else:
        pyscf_driver_mp2no(out_prefix,True,diag_type,appnos,
                     nread_from_h5,h5_add_orbs,nskip,eigcut,
                     mp2noecut,nextracut,kappa,regp)

def qe_driver_hamil(qe_info,out_prefix='pwscf',
                    nread_from_h5=0,h5_add_orbs='',ndet=1,eigcut=1e-3,
                    nextracut=1e-6,thresh=1e-5,ncholmax=15,get_hf=True,
                    get_mp2=True,update_qe_bands=False):

# subroutine pyscf_driver_hamil(out_prefix_, nread_from_h5, h5_add_orbs_, &
#       ndet, eigcut, nextracut, thresh, ncholmax, &
#       get_hf, get_mp2, update_qe_bands, ehf, emp2)
    ehf, emp2 = pyscf_driver_hamil(out_prefix, nread_from_h5, h5_add_orbs, 
             ndet, eigcut, nextracut, thresh, ncholmax, 
             get_hf, get_mp2, update_qe_bands)
    return ehf, emp2
    



