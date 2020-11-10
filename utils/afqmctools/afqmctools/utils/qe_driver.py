from h5py import File
import numpy
from pyscf.pbc import gto, tools
from pyscf.pbc.dft import numint
from pyscf import gto as molgto
import os
import sys
import numpy
from mpi4py import MPI
from afqmctools.utils.gto_basis_utils import extend_gto_id
try:
    from pyscf_driver import (pyscf_driver_init, pyscf_driver_get_info, pyscf_driver_end,
                    pyscf_driver_mp2,pyscf_driver_hamil,pyscf_driver_mp2no)
except ImportError:
    print("Warning: module pyscf_driver not found. AFQMC QE converter required to "
          "use this module.")

def make_cell(latt,sp_label,atid,atpos,basis_='gthdzvp',pseudo_='gthpbe',mesh=None,prec=1e-8):
    """Generates a PySCF gto.Cell object.

    Parameters
    ----------
    latt: (3,3) floating point array.
      Lattice vectors in Bohr. Used to define cell.a
    sp_label: array of strings.
      Array containing species labels/symbols. 
    atid: integer array. 
      Array that contains the mapping between the atoms in the unit cell and their
      species label (as defined by sp_label).
    basis_: string. Default: 'gthdzvp'
      Basis set string. Used to define cell.basis.
    pseudo_: string. Default: 'gthpbe'
      Pseudopotential string. Used to define cell.pseudo.
    mesh: 3-d array. Default: None
      Used to define cell.mesh.
    prec: floating point. Default: 1e-8
      Used to define cell.precision     

    Returns
    -------
    cell: PySCF get.Cell object.
    """
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

def write_esh5_orbitals(cell, name, kpts = numpy.zeros((1,3),dtype=numpy.float64)):
    """Writes periodic AO basis to hdf5 file.  

    Parameters
    ----------
    cell: PySCF get.Cell object
      PySCF cell object which contains information of the system, including 
      AO basis set, FFT mesh, unit cell information, etc.
    name: string 
      Name of hdf5 file.
    kpts: array. Default: numpy.zeros((1,3)
      K-point array of dimension (nkpts, 3)
    dtype: datatype. Default: numpy.float64
      Datatype of orbitals in file.   

    """

    def to_qmcpack_complex(array):
        shape = array.shape
        return array.view(numpy.float64).reshape(shape+(2,))
    nao = cell.nao_nr()

    fh5 = File(name,'w')
    coords = cell.gen_uniform_grids(cell.mesh)
    
    kpts = numpy.asarray(kpts)
    nkpts = len(kpts)
    norbs = numpy.zeros((nkpts,),dtype=int)
    norbs[:] = nao

    grp = fh5.create_group("OrbsG")
    dset = grp.create_dataset("reciprocal_vectors", data=cell.reciprocal_vectors())
    dset = grp.create_dataset("number_of_kpoints", data=len(kpts))
    dset = grp.create_dataset("kpoints", data=kpts)
    dset = grp.create_dataset("number_of_orbitals", data=norbs)
    dset = grp.create_dataset("fft_grid", data=cell.mesh)
    dset = grp.create_dataset("grid_type", data=int(0))
    nnr = cell.mesh[0]*cell.mesh[1]*cell.mesh[2]
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
    """Splits a communicator into image communicators, consistent with QE partitioning.
       nimage consecutive ranks in comm belong to the same image communicator.
       The number of distinct image communcators is comm.size/nimage. 

    Parameters
    ----------
    nimage: integer
      Number of image communicators, must divide comm.size.
    comm: mpi4py MPI comunicator. Default: MPI.COMM_WORLD
      A valid mpi4py communicator. 

    Returns
    -------
    intra_image: mpi4py MPI comunicator.
      Communicator between mpi tasks within an image.
    inter_image:mpi4py MPI comunicator.
      Communicator between mpi tasks on different images, but having the same rank in the image.
    """
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
    """Initializes the QE driver. Must be called before any routine that calls
       the QE driver is executed. Requires a pre-existing QE successful run.

    Parameters
    ----------
    norb: integer
        Number of orbitals read from QE calculation.
    qe_prefix: string
        prefix parameter from QE run.
    qe_outdir: string
        outdir parameter from QE run. (location of QE files).
    atm_labels: array of strings
        Array containing the species labels.
    intra_image: mpi4py communicator. Default: MPI.COMM_WORLD
        Intra image communicator.
    inter_image: mpi4py communicator. Default: None
        Inter image communicator.
    npools: integer. Default: 1
        Number of QE pools used in the driver.
    outdir: string. Default: ./qedrv
        Output directory of the driver. Does not need to be the same as the QE parameter.
    set_soft_links: Bool. Default: True
        If true, soft links to QE files/foulders from qe_outdir will be placed in outdir.
    verbose: Bool. Default: True
        Sets verbosity in driver.
    add_image_tag: Bool. Default True.
        If True, outdir is modified by adding a tag that identifies the image.
        This is needed if running with multiple images simultaneously,
        otherwise the files from different images might conflict with each other.
    Returns
    -------
    qe_info: Python Dictionary
      Dictonary containing all stored information about the QE driver.
      Contents:
        'species' : string array    # array with species labels.
        'nsp' : integer,            # number of species
        'nat' : integer,            # number of atoms
        'at_id' :  integer array    # array with the ids of atoms in the unit cell.  
        'at_pos' : (nat,3) fp array # array with atom positions
        'nkpts' : integer           # number of kpoints
        'kpts' : (nkpts,3) fp array # k-points
        'latt' : (3,3) fp array     # lattice vectors
        'npwx' : integer            # npwx parameter from QE.
        'mesh' : 3D integer array   # FFT mesh 
        'ngm' : integer             # ngm parameter from QE. 
        'outdir' : string           # Location of folder with driver files. 
    """
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
    if intra_rank == 0:
#        if remove_dir:
#            os.system('rm -rf '+fname+'\n')
        if set_soft_links:
            os.system('mkdir '+fname)
            os.system('ln -s ./'+qe_outdir+'/'+qe_prefix+'.xml '+fname+'/'+qe_prefix+'.xml')
            os.system('ln -s ./'+qe_outdir+'/'+qe_prefix+'.save/ '+fname+'/'+qe_prefix+'.save')
    MPI.COMM_WORLD.barrier()

    # initialize driver:
    nkpts, nat, nsp, npwx, ngm, mesh = pyscf_driver_init(inter_size, npools, intra_size/npools,
                                                         norb, qe_prefix, fname, verbose)
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
    """Finishes and performs clean-up on the QE driver.  
       After a call to this routine, further calls to the driver are undefined.
    """
    if(MPI.COMM_WORLD.rank==0):
        print(" Closing QE driver.")
    pyscf_driver_end()

def gen_qe_gto(qe_info,bset,x=[],
               fname='pyscf.orbitals.h5',prec=1e-12):
    """ Writes periodic AO basis set in real space to hdf5 file.
        This routine constructs a new gaussian basis set from a OptimizableBasisSet
        object and an array of optimizable parameters. 
        With the resulting basis set, a new gto.Cell object is constructed consistent 
        with the QE calculation and used to generate the periodic AO basis set.

    Parameters
    ----------
    qe_info: Python Dictionary.
        Dictionary with information from QE calculation, generated by qe_driver_init.
    bset: Object of type OptimizableBasisSet.
        Contains information about a (possibly dynamic) basis set.
    x: fp array. Default: [] (no variable parameters in bset)
        Array with variational parameters in the bset object.
    fname: string. Default: 'pyscf.orbitals.h5'
        Name of hdf5 file.
    prec: floating point number. Default: 1e-12
        Precision used to generate AO orbitals in real space. Controls sum over periodic images. 

    Returns
    -------
    nao: integer
        Number of atomic orbitals generates.
    """
    assert(len(x) == bset.number_of_params)
    basis = {}
    for I, atm in enumerate(qe_info['species']):
        basis.update({atm: molgto.parse( bset.basis_str(atm,x) )})

    cell = make_cell(qe_info['latt'],
                 qe_info['species'],
                 qe_info['at_id'],
                 qe_info['at_pos'],
                 basis_=basis,
                 mesh=qe_info['mesh'],prec=prec)
    nao = cell.nao_nr()
    write_esh5_orbitals(cell, qe_info['outdir']+fname, kpts=qe_info['kpts'])
    return nao

def qe_driver_MP2(qe_info,out_prefix='pyscf_drv',
                        diag_type='keep_occ',
                        nread_from_h5=0,h5_add_orbs='',
                        eigcut=1e-3,nextracut=1e-6,kappa=0.0,regp=0):
    """ Calls the MP2 routine in the driver. 

    Parameters
    ----------
    qe_info: Python Dictionary.
        Dictionary with information from QE calculation, generated by qe_driver_init.
    out_prefix: string. Default: 'pyscf_drv'
        Prefix used in all the files generated by the driver.
    diag_type: string. Default: 'keep_occ'
        Defines the type of HF diagonalization performed before the MP2 calculation.
        Options:
            'keep_occ': Only the virtual orbitals/eigenvalues are calculated.
                        Occupied orbitals/eigenvalues are kept from the QE calculation.
            'full': All orbitals/eigenvalues are recalculated.
            'fullpw': A basis set is generated that contains all the plane waves 
                      below the QE wfn cutoff. The HF eigenvalues/orbitals and MP2NO 
                      are calculated in this basis.
    nread_from_h5: integer. Default: 0
        Number of orbitals to read from h5_add_orbs.  
    h5_add_orbs: string. Default: ''
        Name of hdf5 file with additional orbitals to add to the basis set.
    eigcut: fp number. Default: 1e-3
        Cutoff used during the generation of the spin independent basis in UHF/GHF
        calculations. Only the eigenvalues of the overlap matrix (alpha/beta) 
        above this cutoff are kept in the calculation. In order to reproduce
        the UHF/GHF energy accurately, this number must be set to a small value (e.g. 1e-8). 
    nextracut: fp number. Default: 1e-6
        Cutoff used when adding states from h5_add_orbs to the basis set.
        When a new state from the file is being added to the orbital set, 
        the component along all current orbitals in the set is removed.
        The resulting (orthogonal) state is added only if the norm of the unnormalized
        orbital is larger than nextracut (state is afterwards normalized).
        This is used as a way to remove linear dependencies from the basis set. 
    """
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
    """ Calls the MP2NO routine in the driver. 

    Parameters
    ----------
    qe_info: Python Dictionary.
        Dictionary with information from QE calculation, generated by qe_driver_init.
    out_prefix: string. Default: 'pyscf_drv'
        Prefix used in all the files generated by the driver.
    appnos: Bool. Default: False.
        If True, generates approximate natural orbitals.
    diag_type: string. Default: 'keep_occ'
        Defines the type of HF diagonalization performed before the MP2 calculation.
        Options:
            'keep_occ': Only the virtual orbitals/eigenvalues are calculated.
                        Occupied orbitals/eigenvalues are kept from the QE calculation.
            'full': All orbitals/eigenvalues are recalculated.
            'fullpw': A basis set is generated that contains all the plane waves 
                      below the QE wfn cutoff. The HF eigenvalues/orbitals and MP2NO 
                      are calculated in this basis.
    nread_from_h5: integer. Default: 0
        Number of orbitals to read from h5_add_orbs.  
    h5_add_orbs: string. Default: ''
        Name of hdf5 file with additional orbitals to add to the basis set.
    nskip: integer. Default: 0
        Number of states above the HOMO state of the solid to skip
        during the calculation of MP2 NOs. This can be used to avoid divergencies
        in metals. The assumption being that these states will be included in the 
        orbital set directly.
    eigcut: fp number. Default: 1e-3
        Cutoff used during the generation of the spin independent basis in UHF/GHF
        calculations. Only the eigenvalues of the overlap matrix (alpha/beta) 
        above this cutoff are kept in the calculation. In order to reproduce
        the UHF/GHF energy accurately, this number must be set to a small value (e.g. 1e-8). 
    nextracut: fp number. Default: 1e-6
        Cutoff used when adding states from h5_add_orbs to the basis set.
        When a new state from the file is being added to the orbital set, 
        the component along all current orbitals in the set is removed.
        The resulting (orthogonal) state is added only if the norm of the unnormalized
        orbital is larger than nextracut (state is afterwards normalized).
        This is used as a way to remove linear dependencies from the basis set. 
    mp2noecut: fp number. Default: 1e-6
        Cutoff used when adding natural orbitals from the MP2 RDM,
        only states with eigenvalue > mp2noecut will be kept.
        If this number is < 0.0, then a specific number of states is kept and is 
        given by nint(-mp2noecut). 
 
    """
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
    """ Calls the MP2 routine in the driver. 

    Parameters
    ----------
    qe_info: Python Dictionary.
        Dictionary with information from QE calculation, generated by qe_driver_init.
    out_prefix: string. Default: 'pyscf_drv'
        Prefix used in all the files generated by the driver.
    nread_from_h5: integer. Default: 0
        Number of orbitals to read from h5_add_orbs.  
    h5_add_orbs: string. Default: ''
        Name of hdf5 file with additional orbitals to add to the basis set.
    ndet: integer. Default: 1
        Maximum number of determinants allowed in the trial wave-function.
    eigcut: fp number. Default: 1e-3
        Cutoff used during the generation of the spin independent basis in UHF/GHF
        calculations. Only the eigenvalues of the overlap matrix (alpha/beta) 
        above this cutoff are kept in the calculation. In order to reproduce
        the UHF/GHF energy accurately, this number must be set to a small value (e.g. 1e-8). 
    nextracut: fp number. Default: 1e-6
        Cutoff used when adding states from h5_add_orbs to the basis set.
        When a new state from the file is being added to the orbital set, 
        the component along all current orbitals in the set is removed.
        The resulting (orthogonal) state is added only if the norm of the unnormalized
        orbital is larger than nextracut (state is afterwards normalized).
        This is used as a way to remove linear dependencies from the basis set. 
    thresh: floating point. Detault: 1e-5
        Value used to stop the iterative calculation of Cholesky vectors. The iterations
        stop when the error on a diagonal element falls below this value. 
    ncholmax: integer. Default: 15
        Maximum number of Cholesky vectors allowed (in units of the number of orbitals).
        If the iterative calculation has not converged when this number of Cholesky vectors
        is found, the calculation stops.
    get_hf: Bool. Default: True
        If True, calculate the HF eigenvalues/eigenvectors.
    get_mp2: Bool. Default: True
        If True, calculate the MP2 energy. (If True, get_hf will be set to true.)
    update_qe_bands: Bool. Default: False
        If True, the orbitals in the QE restart file will overwritten with the
        basis set generated by the driver. Orbitals beyond norb (set in qe_driver_init)
        will be left unmodified. 
    Returns
    -------
    ehf: floting point
        The HF energy on this basis. (return 0.0 if not requested)
    emp2: floting point
        The MP2 energy on this basis. (return 0.0 if not requested)
    """
    if get_mp2:
        get_hf = True
    ehf, emp2 = pyscf_driver_hamil(out_prefix, nread_from_h5, h5_add_orbs,
             ndet, eigcut, nextracut, thresh, ncholmax,
             get_hf, get_mp2, update_qe_bands)
    return ehf, emp2
