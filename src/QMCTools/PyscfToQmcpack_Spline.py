######################################################################################
## This file is distributed under the University of Illinois/NCSA Open Source License.
## See LICENSE file in top directory for details.
##
## Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
##
## File developed by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, University of Illinois Urbana-Champaign
##                    Thomas Applencourt, applencourt@anl.gov,  Argonne National Laboratory
##                    Hyeondeok Shin, hshin@anl.gov, Argonne National Laboratory          
##                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
##                  
## File created by: Thomas Applencourt, applencourt@anl.gov, Argonne National Laboratory
#######################################################################################

import numpy as np
import h5py
import os
try:
  from lxml import etree
except:
  import sys
  sys.exit("Error: lxml python module is needed for the generation of the XML file by PyscfToQmcpack_Spline.py. Install as other packages, either directly or with a package manager.")

try:
  import pandas as pd
except:
  import sys
  sys.exit("Error: Pandas python module is needed for the save_eigensystem and eigensystem functions by PyscfToQmcpack_Spline.py. Install as other packages, either directly or with a package manager.")


def pyscf2qmcpackspline(cell,mf,title="Default", kpts=[], kmesh=[],  sp_twist=[]):
  import sys, re

  Restricted=True
  PBC=False
  Gamma=False

  # Python version check
  if sys.version_info < (3, 2):
       sys.exit("Python < 3.2 not supported")

  # FFT mesh check
  if np.any(cell.mesh % 2 == 0):
       sys.exit("Even number of FFT mesh not supported")	
    
  #Twists generation not yet implemented
  if len(sp_twist)== 0:
       sp_twist=[0.0,0.0,0.0]

  val=str(mf)
  ComputeMode= re.split('[. ]',val)

  SizeMode=len(ComputeMode)
  for i in range(SizeMode):
     if ComputeMode[i] in ("UHF","KUHF","UKS"):
           Restricted=False
           sys.exit("Unrestricted calculations not supported")

     if ComputeMode[i]=="pbc":
           PBC = True


  if not PBC:
        sys.exit("Open boundary condition without lattice vectors not supported")
     
  if PBC and len(kpts) == 0:
        Gamma=True

  if len(kpts)!= 0:
     sys.exit("K-point scf not supported")
  else:
     loc_cell=cell
 
  h5_fname = "{}.h5".format(title)
  xml_fname = "{}.xml".format(title)

  tile = [1, 1, 1]
  tilematrix = [0]*9
  for dim, val in enumerate(tile):
    tilematrix[4*dim] = val
   

  tilematrix_str = " ".join(map(str,tilematrix))

  gvecs, eig_df = save_eigensystem(mf, save=False)

  # generate wave function file
  # ================================================
  generate_pwscf_h5(loc_cell,gvecs,eig_df,h5_fname)  

  # generate QMCPACK input file
  # ================================================ 
  h5_handle = h5py.File(h5_fname,'r')
  inp = InputXml()
  # build <simulationcell>
  sc_node = inp.simulationcell_from_cell(loc_cell)
  # build <particleset>
  pset_node = inp.particleset_from_hdf5(h5_handle)
  psedinit_node = inp.particleset_initialposition_from_hdf5(h5_handle)
  wavefunction_node = inp.wavefunction(h5_handle,h5_fname, tilematrix_str)
  hamiltonien_node = inp.hamiltonian(h5_handle, loc_cell)

  # assemble <qmcsystem>
  sys_node = etree.Element('qmcsystem')

  sys_children = [sc_node]
  for child in sys_children:
    sys_node.append(child)


  # write input
  root = etree.Element('simulation')
  doc = etree.ElementTree(root)
  #Change name!!!
  root.append(etree.fromstring('<project id="{}" series="0"/>'.format(title)))
  root.append(sys_node)

  for node in [pset_node,psedinit_node,wavefunction_node,hamiltonien_node]:
      root.append(node)

  a = inp.vmc_dmc(h5_handle)
  root.extend(a)

  doc.write(xml_fname,pretty_print=True)
  print ('Wavefunction successfully saved to QMCPACK HDF5 spline format')



def get_supercell(cell,kmesh=[]):
  latt_vec = cell.lattice_vectors()
  if len(kmesh)==0:
      # Guess kmesh
      scaled_k = cell.get_scaled_kpts(kpts).round(8)
      kmesh = (len(np.unique(scaled_k[:,0])),
               len(np.unique(scaled_k[:,1])),
               len(np.unique(scaled_k[:,2])))

  R_rel_a = np.arange(kmesh[0])
  R_rel_b = np.arange(kmesh[1])
  R_rel_c = np.arange(kmesh[2])
  R_vec_rel = lib.cartesian_prod((R_rel_a, R_rel_b, R_rel_c))
  R_vec_abs = np.einsum('nu, uv -> nv', R_vec_rel, latt_vec)


  # R_rel_mesh has to be construct exactly same to the Ts in super_cell function
  scell = tools.super_cell(cell, kmesh)
  return scell, kmesh
#end def get_supercell


def save_eigensystem(mf,gvec_fname = 'gvectors.dat'
                     ,eigsys_fname = 'eigensystem.json',save=True):
  import os
  if os.path.isfile(eigsys_fname) and os.path.isfile(gvec_fname):
    gvecs = np.loadtxt(gvec_fname)
    eig_df = pd.read_json(eigsys_fname).set_index(
      ['ikpt','ispin','istate'],drop=True).sort_index()
  else:
    data = []
    ikpt  = 0 # gamma-point calculation
    ispin = 0 # restricted (same orbitals for up and down electrons)
    # get MOs in plane-wave basis
    aoR = ao_on_grid(mf.cell)
    gvecs,psig = mo_coeff_to_psig(mf.mo_coeff,aoR,mf.cell.gs,mf.cell.vol)
    nstate,npw,ncomp = psig.shape
    for istate in range(nstate):
      entry = {'ikpt':ikpt,'ispin':ispin,'istate':istate,
        'reduced_k':mf.kpt,'evalue':mf.mo_energy[istate],'evector':psig[istate,:,:]}
      data.append(entry)
    # end for istate
    eig_df = pd.DataFrame(data).set_index(
      ['ikpt','ispin','istate'],drop=True).sort_index()
  # end if
  if save:
    eig_df.reset_index().to_json(eigsys_fname)
    np.savetxt(gvec_fname,gvecs)
  # end if
  return gvecs,eig_df
# end def save_eigensystem

def ao_on_grid(cell):
  from pyscf.pbc.dft import gen_grid,numint
  coords = gen_grid.gen_uniform_grids(cell)
  aoR    = numint.eval_ao(cell,coords)
  return aoR
# end def ao_on_grid

def mo_coeff_to_psig(mo_coeff,aoR,cell_gs,cell_vol,int_gvecs=None):
  """
   Inputs:
     mo_coeff: molecular orbital in AO basis, each column is an MO, shape (nao,nmo)
     aoR: atomic orbitals on a real-space grid, each column is an AO, shape (ngrid,nao)
     cell_gs: 2*cell_gs+1 should be the shape of real-space grid (e.g. (5,5,5))
     cell_vol: cell volume, used for FFT normalization
     int_gvecs: specify the order of plane-waves using reciprocal lattice points
   Outputs:
       3. plane-wave coefficients representing the MOs, shape (ngrid,nmo)
  """
  import sys

  # provide the order of reciprocal lattice vectors to skip
  if int_gvecs is None: # use internal order
    nx,ny,nz = cell_gs
    from itertools import product
    int_gvecs = np.array([gvec for gvec in product(
      range(-nx,nx+1),range(-ny,ny+1),range(-nz,nz+1))],dtype=int)
  else:
    assert (int_gvecs.dtype is int)
  # end if
  npw = len(int_gvecs) # number of plane waves 

  # put molecular orbitals on real-space grid
  moR = np.dot(aoR,mo_coeff)
  nao,nmo = moR.shape
  rgrid_shape = 2*np.array(cell_gs)+1
  assert nao == np.prod(rgrid_shape)

  # for each MO, FFT to get psig
  psig = np.zeros([nmo,npw,2]) # store real & complex
  for istate in range(nmo):
    # fill real-space FFT grid
    rgrid = moR[:,istate].reshape(rgrid_shape)
    # get plane-wave coefficients (on reciprocal-space FFT grid)
    moG   = np.fft.fftn(rgrid)/np.prod(rgrid_shape)*np.sqrt(cell_vol)
    orb_norm = np.sum(moG*np.conj(moG)).real

    if abs(1.-orb_norm) > 1.e-6:
      print('Orbital normalization failed in state:'+str(istate)+' with norm:'+str(orb_norm))
      sys.exit(0)

    # transfer plane-wave coefficients to psig in specified order
    for igvec in range(npw):
      comp_val = moG[tuple(int_gvecs[igvec])]
      psig[istate,igvec,:] = comp_val.real,comp_val.imag
    # end for igvec
  # end for istate
  return int_gvecs,psig
# end def mo_coeff_to_psig

def generate_pwscf_h5(cell,gvecs,eig_df,h5_fname):

  # if eigensystem was saved to disk, use the following to read
  #import numpy as np
  #import pandas as pd
  #gvecs = np.loadtxt('../1_eigsys/gvectors.dat')
  #eig_df= pd.read_json('../1_eigsys/eigensystem.json').set_index(
  #  ['ikpt','ispin','istate'],drop=True).sort_index()

  new = h5py.File(h5_fname,'w')
  ref = PwscfH5()
  nelecs = ref.system_from_cell(new,cell)
  ref.create_electrons_group(new,gvecs,eig_df,nelecs)

  # transfer version info. !!!! hard code for now
  new.create_dataset('application/code',data=[np.string_('PySCF')])
  new.create_dataset('application/version',data=[np.string_('1.7.5')])
  new.create_dataset('format',data=[np.string_('ES-HDF')])
  new.create_dataset('version',data=[2,1,0])
  new.close()
# end def generate_pwscf_h5

# =======================================================================
# Class for bspline h5 generator
# =======================================================================
class PwscfH5:
  def __init__(self):
    self.locations = {
      'gvectors':'electrons/kpoint_0/gvectors',
       'nkpt':'electrons/number_of_kpoints',
       'nspin':'electrons/number_of_spins',
       'nstate':'electrons/kpoint_0/spin_0/number_of_states', # !!!! same number of states per kpt
       'axes':'supercell/primitive_vectors'
        }
    self.dtypes = {
        'nkpt':int,
        'nspin':int,
        'nstate':int
        }
    self.fp = None # h5py.File object (like a file pointer)
  def __del__(self):
    if self.fp is not None:
      self.fp.close()

  # =======================================================================
  # Basic Read Methods i.e. basic read/write and path access
  # =======================================================================
  def read(self,fname,force=False):
    """ open 'fname' for reading and save handle in this class """
    if not os.path.isfile(fname):
      raise RuntimeError('%s not found' % fname)
    if (self.fp is None) or force:
      self.fp = h5py.File(fname)
    else:
      raise RuntimeError('already tracking a file %s'%str(self.fp))

  def val(self,loc):
    """ get value array of an arbitrary entry at location 'loc' """
    return self.fp[loc][()]

  def get(self,name):
    """ get value array of a known entry """
    loc   = self.locations[name]
    return self.fp[loc][()]

  # =======================================================================
  # Advance Read Methods i.e. more specific to QMCPACK 3.0.0
  # =======================================================================

  # construct typical paths
  #  e.g. electrons/kpoint_0/spin_0/state_0
  @staticmethod
  def kpoint_path(ikpt):
    path = 'electrons/kpoint_%d' % (ikpt)
    return path
  @staticmethod
  def spin_path(ikpt,ispin):
    path = 'electrons/kpoint_%d/spin_%d' % (ikpt,ispin)
    return path
  @staticmethod
  def state_path(ikpt,ispin,istate):
    path = 'electrons/kpoint_%d/spin_%d/state_%d/' % (ikpt,ispin,istate)
    return path

  # access specific eigenvalue or eigenvector
  def psig(self,ikpt=0,ispin=0,istate=0):
    psig_loc = self.state_path(ikpt,ispin,istate)+'psi_g'
    return self.fp[psig_loc][()]

  def psir(self,ikpt=0,ispin=0,istate=0):
    psir_loc = self.state_path(ikpt,ispin,istate)+'psi_r'
    return self.fp[psir_loc][()]

  def eigenvalues(self):
    """ return all eigenvalues, shape=(nkpt,nspin,nstate) """
    nkpt   = self.get('nkpt')[0]
    nspin  = self.get('nspin')[0]
    nstate = self.get('nstate')[0] # !!!! same number of states per kpt

    evals  = np.zeros([nkpt,nspin,nstate])
    for ikpt in range(nkpt):
      for ispin in range(nspin):
        path = self.spin_path(ikpt,ispin)
        evals[ikpt,ispin,:] = self.val(
        os.path.join(path,'eigenvalues')
        )
    return evals

  @classmethod
  def psig_to_psir(self,gvecs,psig,rgrid_shape,vol):
    """ contruct orbital given in planewave basis
    Inputs: 
     gvecs: gvectors in reciprocal lattice units i.e. integers
     psig: planewave coefficients, should have the same length as gvecs
     vol: simulation cell volume, used to normalized fft
    Output:
     rgrid: orbital on a real-space grid """
    assert len(gvecs) == len(psig)

    kgrid = np.zeros(rgrid_shape,dtype=complex)
    for igvec in range(len(gvecs)):
      kgrid[tuple(gvecs[igvec])] = psig[igvec]
    # end for
    rgrid = np.fft.ifftn(kgrid) * np.prod(rgrid_shape)/vol
    return rgrid
  # end def psig_to_psir

  def get_psir_from_psig(self,ikpt,ispin,istate,rgrid_shape=None,mesh_factor=1.0):
    """ FFT psig to psir at the given (kpoint,spin,state) """
    # get lattice which defines the FFT grid
    axes = self.get('axes')
    vol  = np.dot(np.cross(axes[0],axes[1]),axes[2])
    # get MO in plane-wave basis
    gvecs = self.get('gvectors').astype(int)
    psig_arr = self.psig(ikpt=ikpt,ispin=ispin,istate=istate)
    psig = psig_arr[:,0] + 1j*psig_arr[:,1]
    # determine real-space grid size (QMCPACK 3.0.0 convention)
    #  ref: QMCWaveFunctions/Experimental/EinsplineSetBuilder.cpp::ReadGvectors_ESHDF()
    if rgrid_shape is not None: # !!!! override grid size
      pass
    else:
      rgrid_shape = map(int, np.ceil(gvecs.max(axis=0)*4*mesh_factor) )
    # end if
    psir = self.psig_to_psir(gvecs,psig,rgrid_shape,vol)
    return psir
  # end def get_psir_from_psig

  # build entire eigensystem as a dataframe
  def eigensystem(self):
    """ construct dataframe containing eigenvalues and eigenvectors
    labeled by (kpoint,spin,state) indices """
    import pandas as pd

    data = []
    nkpt = self.get('nkpt')
    nspin= self.get('nspin')
    for ikpt in range(nkpt):
      k_grp = self.fp[self.kpoint_path(ikpt)]
      rkvec = k_grp['reduced_k'][()]
      for ispin in range(nspin):
        spin_loc = self.spin_path(ikpt,ispin)
        sp_grp   = self.fp[spin_loc]
        nstate   = sp_grp['number_of_states'][(0)]
        evals    = sp_grp['eigenvalues'][()]
        for istate in range(nstate):
          st_loc = self.state_path(ikpt,ispin,istate)
          st_grp = self.fp[st_loc]
          evector= st_grp['psi_g'][()] # shape (ngvec,2) (real,complex)
          entry  = {'ikpt':ikpt,'ispin':ispin,'istate':istate,
            'reduced_k':rkvec,'evalue':evals[istate],'evector':evector}
          data.append(entry)
        # end for istate
      # end for ispin
    # end for ikpt
    df = pd.DataFrame(data).set_index(['ikpt','ispin','istate'],drop=True)
    return df
  # end def eigensystem

  # =======================================================================
  # Advance Write Methods, some specialized for pyscf
  # =======================================================================
  @staticmethod
  def create_electrons_group(h5_handle,gvec,df,nelec):
    """ create and fill the /electrons group in hdf5 handle
    Inputs:
      h5_handle: hdf5 handle generated by h5py.File
      gvec: 2D numpy array of reciprocal space vectors (npw,ndim)
      df: dataframe containing the eigensystem,
         indexed by (kpt,spin,state), contains (evalue,evector,reduced_k)
      nelec: a list of the number of electrons per atom (if no pseudopotential, then 'species_id' returned by system_from_cell should do)
    Output:
      None
    Effect:
      fill /electrons group in 'h5_handle' """
    flat_df = df.reset_index()
    kpoints = flat_df['ikpt'].unique()
    spins   = flat_df['ispin'].unique()
    nkpt,nspin = len(kpoints),len(spins)
    # transfer orbitals (electrons group)
    for ikpt in range(nkpt):
      # !!!! assume no symmetry was used to generate the kpoints
      kpt_path = 'electrons/kpoint_%d'%ikpt
      kgrp = h5_handle.create_group(kpt_path)
      kgrp.create_dataset('num_sym',data=[1])
      kgrp.create_dataset('symgroup',data=[1])
      kgrp.create_dataset('weight',data=[1])

      rkvec = df.loc[ikpt,'reduced_k'].values[0]
      kgrp.create_dataset('reduced_k',data=rkvec)
      if ikpt == 0: # store gvectors in kpoint_0
        kgrp.create_dataset('gvectors',data=gvec)
        kgrp.create_dataset('number_of_gvectors',data=[len(gvec)])
      # end if 

      for ispin in range(nspin): # assume ispin==0
        nstate = len(df.loc[(ikpt,ispin)])
        spin_path = os.path.join(kpt_path,'spin_%d'%ispin)
        spgrp     = h5_handle.create_group(spin_path)
        spgrp.create_dataset('number_of_states',data=[nstate])

        evals = np.zeros(nstate) # fill eigenvalues during eigenvector read
        for istate in range(nstate):
          state_path = os.path.join(spin_path,'state_%d'%istate)
          psig = df.loc[(ikpt,ispin,istate),'evector']
          psig_path = os.path.join(state_path,'psi_g')
          h5_handle.create_dataset(psig_path,data=psig)
          evals[istate] = df.loc[(ikpt,ispin,istate),'evalue']
        # end for istate
        spgrp.create_dataset('eigenvalues',data=evals)
      # end for ispin
    # end for ikpt
    # transfer orbital info
    h5_handle.create_dataset('electrons/number_of_electrons',data=nelec)
    h5_handle.create_dataset('electrons/number_of_kpoints',data=[nkpt])
    # !!!! hard-code restricted orbitals
    h5_handle.create_dataset('electrons/number_of_spins',data=[1])
  # end def create_electrons_group

  @staticmethod
  def system_from_cell(h5_handle,cell):
    """ create and fill the /supercell and /atoms groups
     Inputs:
       h5_handle: hdf5 handle generated by h5py.File
       cell: pyscf.pbc.gto.Cell class
     Outputs:
       species_id: a list of atomic numbers for each atom
     Effect:
       fill /supercell and /atoms group in 'h5_handle'
    """

    # write lattice
    axes = cell.lattice_vectors() # always in bohr
    h5_handle.create_dataset('supercell/primitive_vectors',data=axes)


    # write atoms
    pos  = cell.atom_coords() # always in bohr
    elem = [cell.atom_symbol(i) for i in range(cell.natm)]
    assert len(pos) == len(elem)
    h5_handle.create_dataset('atoms/number_of_atoms',data=[len(elem)])
    h5_handle.create_dataset('atoms/positions',data=pos)

    # write species info
    species, indices  = np.unique(elem, return_index=True)
    h5_handle.create_dataset('atoms/number_of_species',data=[len(species)])
    atomic_number = {}
    number_of_electrons = {}
    species_map = {}
    nelecs = cell.nelec

    for ispec,name in enumerate(species):
      species_map[name] = ispec
      atomic_number[name] = cell.atom_nelec_core(indices[ispec]) + cell.atom_charge(indices[ispec])
      number_of_electrons[name] = cell.atom_charge(indices[ispec])
      spec_grp = h5_handle.create_group('atoms/species_%d'%ispec)

      # write name
      if name not in species_map.keys():
        raise NotImplementedError('unknown element %s' % name)
      # end if
      spec_grp.create_dataset('name',data=[np.string_(name)])

      # write atomic number and valence
      Zn   = atomic_number[name]
      spec_grp.create_dataset('atomic_number',data=[Zn])
      Zps = number_of_electrons[name]
      spec_grp.create_dataset('valence_charge',data=[Zps])
    # end for ispec 
    species_ids = [species_map[name] for name in elem]
    h5_handle.create_dataset('atoms/species_ids',data=species_ids)

    return nelecs
  # end def system_from_cell

# end class PwscfH5


# =======================================================================
# Class for xml generator
# =======================================================================
class InputXml:
  def __init__(self):
    pass
  # end def

  # =======================================================================
  # Basic Methods (applicable to all xml files)
  # =======================================================================
  def read(self,fname):
    self.fname = fname
    parser = etree.XMLParser(remove_blank_text=True)
    self.root = etree.parse(fname,parser)
  # end def

  def write(self,fname=None,pretty_print=True):
    if fname is None:
      self.root.write(self.fname,pretty_print=pretty_print)
    else:
      self.root.write(fname,pretty_print=pretty_print)
    # end if
  # end def

  def show(self,node):
    """ print text representation of an xml node """
    print (etree.tostring(node,pretty_print=True))

  # pass along xpath expression e.g. './/particleset'
  def find(self,xpath):
    return self.root.find(xpath)
  # end def

  def find_all(self,xpath):
    return self.root.findall(xpath)
  # end def


  @classmethod
  def arr2text(self,arr):
    """ format convert a numpy array into a text string """
    text = ''
    if len(arr.shape) == 1: # vector
      text = " ".join(arr.astype(str))
    elif len(arr.shape) == 2: # matrix
      mat  = [self.arr2text(line) for line in arr]
      text = "\n      " + "\n      ".join(mat) + "\n"
    else:
      raise RuntimeError('arr2text can only convert vector or matrix.')
    # end if
    return text
  # end def

  @classmethod
  def text2arr(self,text,dtype=float,flatten=False):
    tlist = text.strip(' ').strip('\n').split('\n')
    if len(tlist) == 1:
      return np.array(tlist,dtype=dtype)
    else:
      if flatten:
        mytext = '\n'.join(['\n'.join(line.split()) for line in tlist])
        myarr = self.text2arr(mytext)
        return myarr.flatten()
      else:
        return np.array([line.split() for line in tlist],dtype=dtype)
      # end if
    # end if
  # end def

  @classmethod
  def node2dict(self,node):
    entry = dict(node.attrib)
    if node.text:
      entry.update({'text':node.text})
    # end if
    return entry
  # end def node2dict

  # =======================================================================
  # Simple Methods Specific to QMCPACK
  # =======================================================================
  def find_pset(self,name='e'):
    """ return xml node specifying the particle set with given name
     by default return the quantum particle set 'e' """
    return self.find('.//particleset[@name="%s"]'%name)
  # end find_pset

  # =======================================================================
  # Advance Methods i.e. specific to pyscf or QMCPACK 3.0
  # =======================================================================

  # ----------------
  # simulationcell
  def simulationcell_from_cell(self,cell,bconds='p p p',lr_cut=15.0):
    """ construct the <simulationcell> xml element from pyscf.pbc.gto.Cell class
     Inputs:
       cell: pyscf.pbc.gto.Cell class, should have lattice_vectors() and unit
       bconds: boundary conditions in each of the x,y,z directions, p for periodic, n for non-periodic, default to 'p p p ' 
       lr_cut: long-range cutoff parameter rc*kc, default to 15
     Output: 
       etree.Element representing <simulationcell>
     Effect:
       none
    """

    # write primitive lattice vectors
    axes = cell.lattice_vectors() # rely on pyscf to return a.u.
    lat_node = etree.Element('parameter'
      ,attrib={'name':'lattice','units':'bohr'})
    lat_node.text = self.arr2text(axes) + "      "

    # write boundary conditions
    bconds_node = etree.Element('parameter',{'name':'bconds'})
    bconds_node.text = bconds

    # write long-range cutoff parameter
    lr_node = etree.Element('parameter',{'name':'LR_dim_cutoff'})
    lr_node.text = str(lr_cut)

    # build <simulationcell>
    sc_node = etree.Element('simulationcell')
    sc_node.append(lat_node)
    sc_node.append(bconds_node)
    sc_node.append(lr_node)
    return sc_node
  # end def simulationcell_from_cell
  # ----------------

  def particleset_from_hdf5(self,h5_handle):
    atom_grp = h5_handle.get('atoms')
    nspec = atom_grp.get('number_of_species')[(0)]
    species_ids = atom_grp.get('species_ids')[()]
    positions = atom_grp.get('positions')[()]

    natom_total = atom_grp.get('number_of_atoms')[(0)]
    # Get the name of the atoms
    groups = []
    for ispec in range(nspec):
      # turn h5 group into dictionary (i.e. h5ls -d)
      sp_grp         = atom_grp.get('species_%d'%ispec)
      name           = sp_grp.get('name')[(0)]
      valence_charge = sp_grp.get('valence_charge')[(0)]
      atomic_number  = sp_grp.get('atomic_number')[(0)]

      # locate particles of this species
      atom_idx = np.where(species_ids==ispec)
      pos_arr  = positions[atom_idx]
      natom    = len(pos_arr)
      # build xml node
      charge_node = etree.Element('parameter',{'name':'charge'})
      charge_node.text = str(valence_charge)
      valence_node = etree.Element('parameter',{'name':'valence'})
      valence_node.text = str(valence_charge)
      atomic_number_node = etree.Element('parameter',{'name':'atomicnumber'})
      atomic_number_node.text = str(atomic_number)

      grp_children = [charge_node,valence_node,atomic_number_node]
      grp_node = etree.Element('group',{'name':name})
      for child in grp_children:
        grp_node.append(child)
      groups.append(grp_node)


    pos_node = etree.Element('attrib',{'name':'position','datatype':'posArray'})
    pos_node.text = self.arr2text(positions) + "    "
    groups.append(pos_node)

    ionid_node = etree.Element('attrib',{'name':'ionid','datatype':'stringArray'})
    ionid_node.text = self.arr2text(np.array([atom_grp.get('species_{}/name'.format(id_))[(0)]  for id_ in species_ids]))

    groups.append(ionid_node)

    # build <particleset>
    pset_node = etree.Element('particleset',{'name':'ion0','size':str(natom_total)})
    for group in groups:
      pset_node.append(group)

    return pset_node

  def particleset_initialposition_from_hdf5(self,h5_handle):
    pset_node = etree.Element('particleset',{'name':'e','random':"yes", 'randomsrc':'ion0'})

    # size = number of electron up and down
    elec_alpha_beta = h5_handle.get('electrons/number_of_electrons')[()]
    l_name = ("u","d")
    for name, electron in zip(l_name,elec_alpha_beta):
      groupe_node = etree.Element('group',{'name':"{}".format(name),'size':"{}".format(electron)})
      param_node = etree.Element('parameter',{'name':'charge'})
      param_node.text=str(-1)
      groupe_node.append(param_node)
      pset_node.append(groupe_node)
    return pset_node

  def wavefunction(self,h5_handle,h5_path,tilematrix):
    wf_node = etree.Element('wavefunction', {'name':"psi0", 'target':"e"})

    determinantset_node = etree.Element('determinantset', {"type":"einspline",
                                                           "href":"{}".format(h5_path),
                                                           "source":"ion0",
                                                           "tilematrix":"{}".format(tilematrix),
                                                           "twistnum":"0",
                                                           "meshfactor":"1.0"})
    basisset_node = etree.Element('basisset')
    determinantset_node.append(basisset_node)

    atom_grp = h5_handle.get('electrons')
    alpha, beta = atom_grp.get('number_of_electrons')[()]

    # Slaterdet,imamt
    slaterdet_node = etree.fromstring('''
    <slaterdeterminant>
        <determinant id="updet" size="{}" ref="updet">
          <occupation mode="ground" spindataset="0">
          </occupation>
        </determinant>
        <determinant id="downdet" size="{}" ref="downdet">
          <occupation mode="ground" spindataset="0">
          </occupation>
        </determinant>
      </slaterdeterminant>
      '''.format(alpha, beta))

    determinantset_node.append(slaterdet_node)

    wf_node.append(determinantset_node)
    # Jastrow
    jastrow_node = etree.fromstring('''
    <jastrow name="J2" type="Two-Body" function="Bspline" print="yes">
      <correlation speciesA="u" speciesB="u" size="10">
        <coefficients id="uu" type="Array"> 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0</coefficients>
      </correlation>
      <correlation speciesA="u" speciesB="d" size="10">
        <coefficients id="ud" type="Array"> 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0</coefficients>
      </correlation>
    </jastrow>''')
    wf_node.append(jastrow_node)

    # Species_id
    atom_grp = h5_handle.get('atoms')
    species_ids = atom_grp.get('species_ids')[()]
    list_atom = sorted(set(atom_grp.get('species_{}/name'.format(id_))[(0)].decode()  for id_ in species_ids))

    jastrow_node =   etree.Element('jastrow', {'name':"J1", 'type':"One-Body", "function":"Bspline", "print":"yes", "source":"ion0"})

    for atom in list_atom:
      jastrow_node.append(etree.fromstring('''
      <correlation elementType="{}" cusp="0.0" size="10">
        <coefficients id="{}" type="Array"> 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0</coefficients>
      </correlation>'''.format(atom, atom)))

    wf_node.append(jastrow_node)
    return wf_node

  def hamiltonian(self,h5_handle, cell):

    atom_grp = h5_handle.get('atoms')
    species_ids = atom_grp.get('species_ids')[()]
    list_atom = sorted(set(atom_grp.get('species_{}/name'.format(id_))[(0)].decode()  for id_ in species_ids))

    if cell.has_ecp():
      hamiltonian_node = etree.fromstring('''
      <hamiltonian name="h0" type="generic" target="e">
        <pairpot name="ElecElec" type="coulomb" source="e" target="e" physical="true"/>
        <pairpot name="IonIon" type="coulomb" source="ion0" target="ion0"/>
        <pairpot name="PseudoPot" type="pseudo" source="ion0" wavefunction="psi0" format="xml">
        </pairpot>
  </hamiltonian>''')
      # Add the list of ecp file for each atom.
      for element in list_atom:
        pset_node = etree.Element('pseudo',{'elementType':element,'href':"{}.qmcpp.xml".format(element)})
        hamiltonian_node[-1].append(pset_node)
    else:
      hamiltonian_node = etree.fromstring('''
      <hamiltonian name="h0" type="generic" target="e">
        <pairpot name="ElecElec" type="coulomb" source="e" target="e" physical="true"/>
        <pairpot name="IonIon" type="coulomb" source="ion0" target="ion0"/>
        <pairpot type="coulomb" name="ElecIon" source="ion0" target="e"/>
      </hamiltonian>''')

    return hamiltonian_node

  def vmc_dmc(self, h5_handle):

    vmc_init_comment = '''
    Example initial VMC to measure initial energy and variance 
    '''
    vmc_init = '''
    <qmc method="vmc" move="pbyp" checkpoint="-1">
      <estimator name="LocalEnergy" hdf5="no"/>
      <parameter name="warmupSteps">100</parameter>
      <parameter name="blocks">20</parameter>
      <parameter name="steps">50</parameter>
      <parameter name="substeps">8</parameter>
      <parameter name="timestep">0.5</parameter>
      <parameter name="usedrift">no</parameter>
  </qmc>
    '''
    loop_comment = '''
    Example initial VMC optimization 
 
    Number of steps required will be computed from total requested sample 
    count and total number of walkers 
    '''
    loop= '''
    <loop max="4">
      <qmc method="linear" move="pbyp" checkpoint="-1">
        <estimator name="LocalEnergy" hdf5="no"/>
        <parameter name="warmupSteps">100</parameter>
        <parameter name="blocks">20</parameter>
        <parameter name="timestep">0.5</parameter>
        <parameter name="walkers">1</parameter>
        <parameter name="samples">16000</parameter>
        <parameter name="substeps">4</parameter>
        <parameter name="usedrift">no</parameter>
        <parameter name="MinMethod">OneShiftOnly</parameter>
        <parameter name="minwalkers">0.003</parameter>
      </qmc>
  </loop>
    '''
    loop_followup_comment='''
    Example follow-up VMC optimization using more samples for greater accuracy
    '''

    loop_followup = '''
    <loop max="10">
      <qmc method="linear" move="pbyp" checkpoint="-1">
        <estimator name="LocalEnergy" hdf5="no"/>
        <parameter name="warmupSteps">100</parameter>
        <parameter name="blocks">20</parameter>
        <parameter name="timestep">0.5</parameter>
        <parameter name="walkers">1</parameter>
        <parameter name="samples">64000</parameter>
        <parameter name="substeps">4</parameter>
        <parameter name="usedrift">no</parameter>
        <parameter name="MinMethod">OneShiftOnly</parameter>
        <parameter name="minwalkers">0.3</parameter>
      </qmc>
  </loop>
    '''
    # Generate production VMC and DMC
    vmc_comment = '''Production VMC and DMC
    Examine the results of the optimization before running these blocks.
    e.g. Choose the best optimized jastrow from all obtained, put in 
    wavefunction file, do not reoptimize.'''

    vmc = ''' 
    <qmc method="vmc" move="pbyp" checkpoint="-1">
      <estimator name="LocalEnergy" hdf5="no"/>
      <parameter name="warmupSteps">100</parameter>
      <parameter name="blocks">200</parameter>
      <parameter name="steps">50</parameter>
      <parameter name="substeps">8</parameter>
      <parameter name="timestep">0.5</parameter>
      <parameter name="usedrift">no</parameter>
      <!--Sample count should match targetwalker count for DMC. Will be obtained from all nodes.-->
      <parameter name="samples">16000</parameter>
  </qmc>'''

    dmc_comment = ""
    dmc ='''
    <qmc method="dmc" move="pbyp" checkpoint="20">
      <estimator name="LocalEnergy" hdf5="no"/>
      <parameter name="targetwalkers">16000</parameter>
      <parameter name="reconfiguration">no</parameter>
      <parameter name="warmupSteps">100</parameter>
      <parameter name="timestep">0.005</parameter>
      <parameter name="steps">100</parameter>
      <parameter name="blocks">100</parameter>
      <parameter name="nonlocalmoves">yes</parameter>
  </qmc>
    '''
    a = []

    l_xml = (vmc_init, loop, loop_followup, vmc, dmc)
    l_xml_comment =  (vmc_init_comment, loop_comment, loop_followup_comment, vmc_comment, dmc_comment)
    for comment, qmc in zip(l_xml_comment,l_xml):
      a.append(etree.Comment(comment))
      a.append(etree.fromstring(qmc))
    return a

  # ----------------
  # numerics
  # grid
  def radial_function(self,node):
    assert node.tag=='radfunc'

    # read grid definitions ( e.g. np.linspace(ri,rf,npts) for linear grid )
    gnode = node.find('.//grid')      # expected attributes: 
    grid_defs = self.node2dict(gnode) #   type,ri,rf,npts,units
    ri = float(grid_defs['ri'])
    rf = float(grid_defs['rf'])
    npts = int(grid_defs['npts'])
    gtype = grid_defs['type']
    units = grid_defs['units']

    # read 
    dnode = node.find('.//data')
    rval  = self.text2arr(dnode.text,flatten=True) # read as 1D vector
    assert len(rval) == npts
    entry = {'type':gtype,'units':units,'ri':ri,'rf':rf,'npts':npts,'rval':rval}
    return entry
  # end def radial_function
  # ----------------

# end class
