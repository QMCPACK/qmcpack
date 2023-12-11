######################################################################################
## This file is distributed under the University of Illinois/NCSA Open Source License.
## See LICENSE file in top directory for details.
##
## Copyright (c) 2021 QMCPACK developers.
##
## File developed by: Juha Tiihonen, tiihonen@iki.fi, University of Jyvaskyla
##
## File created by: Juha Tiihonen, tiihonen@iki.fi, University of Jyvaskyla 
#######################################################################################

from numpy import array

# This file defines python classes to load/write data needed by the ESHDF file format.
#
# The main logic of the class structure is to inherit from Base -> Representation -> Application
#
# currently implemented:
#
# Base:
#   EshdfFileBase
# Representation:
#   EshdfFilePw
# Application:
#   EshdfFilePwGpaw  <- found in Eshdf_gpaw.py


# define dtypes for the HDF5
int_dtype     = '<i4'   # default integer dtype
str_dtype     = 'S'     # default string dtype
float_dtype   = 'f'     # default float dtype
eshdf_version = [2,1,0] # current version


# Class for atomic species
class AtomicSpecies():
    id      = None
    mass    = None
    Z       = None
    name    = None
    valence = None
    def __init__(
            self,
            id,
            name,
            mass,
            Z,
            valence
            ):
        id           = id
        self.mass    = mass
        self.Z       = Z
        self.name    = namei
        self.valence = valence
    #end def
#end class


# Base class for an ESHDF file
class EshdfFileBase():
    
    # class flags and defaults
    loaded      = False
    outfile     = 'eshdf.h5'

    # init null or default ESHDF data

    # /application must be overridden
    code    = ''      # string
    version = [0,0,1] # 3-int list
    # /atoms must be overridden
    num_atoms   = 0
    num_species = 0
    positions   = [[]] # num_atoms x 3 array
    species     = []   # list of AtomicSpecies (num_species)
    species_ids = []   # list of integers (num_atoms)
    # /electrons must be overridden
    num_elec = 0
    num_spin = 0
    # /format
    fmt     = 'ES-HDF'
    # /supercell must be overridden
    supercell = [[1,0,0],[0,1,0],[0,0,1]] # 3x3 matrix
    # /version
    version = eshdf_version

    # standard init: infile, outfile, **kwargs
    def __init__(
        self,
        infile    = None,
        outfile   = None,
        **kwargs
        ):
        self.infile  = infile
        self.outfile = outfile
        if infile is not None:
            self.load(infile=infile,**kwargs)
        #end if
    #end def

    # standard load: infile, **kwargs
    #   override this in implementations to read data
    def load(
        self,
        infile  = None,
        **kwargs,
        ):
        print('Loading of data not implemented.')
        self.loaded = True
    #end def

    # check integrity of the data
    def check_integrity(self):
        #TODO: various tests
        return True
    #end def

    # standard write: outfile
    #   default: write to HDF5
    def write(self, outfile = None):
        if not self.loaded:
            raise(DataNotLoaded('Must load data first'))
        #end if
        if outfile is None:
            outfile = self.outfile
        #end if
        if self.check_integrity():
            self._write_h5(outfile)
        #end if
    #end def

    # write to HDF5
    def _write_h5(self, outfile):
        import h5py as h5
        with h5.File(outfile,mode='w') as f:
            self._write_h5_application(f)
            self._write_h5_atoms(f)
            self._write_h5_electrons(f)
            self._write_h5_format(f)
            self._write_h5_supercell(f)
            self._write_h5_version(f)
        #end with
    #end def

    def _write_h5_application(self,f):
        g = f.create_group('application')
        g.create_dataset('code',shape=1,data=array(self.application,dtype=str_dtype))
        g.create_dataset('version',shape=3,data=array(self.app_version,dtype=int_dtype))
    #end def

    def _write_h5_atoms(self,f):
        g  = f.create_group('atoms')
        g.create_dataset('number_of_atoms',data=array([self.num_atoms],dtype=int_dtype))
        g.create_dataset('number_of_species',data=array([self.num_species],dtype=int_dtype))
        g.create_dataset('positions',data=self.positions)
        for s,sid in enumerate(self.species):
            gs = g.create_group('species_'+str(s))
            gs.create_dataset('atomic_number',data=array([self.species_Z[s]]),dtype=int_dtype)
            gs.create_dataset('mass',data=array([self.species_mass[s]],dtype=float_dtype))
            gs.create_dataset('name',data=array([self.species_name[s]],dtype=str_dtype))
            gs.create_dataset('valence_charge',data=array([self.species_valence[s]],dtype=float_dtype))
        #end for
        g.create_dataset('species_ids',data=array(self.species_ids),dtype=int_dtype)
    #end def

    # writing of the orbital information must be overridden based on the representation
    #   Pw, real-space, etc
    def _write_h5_electrons(self,f):
        g = f.create_group('electrons')
        g.create_dataset('number_of_electrons',data=self.num_elec)
        g.create_dataset('number_of_spins',data=array([self.num_spin],dtype=int_dtype))
        print('To write kpoint/orbital information, must not use the base class EshdfFileBase\nbut an appropriate representation, such as EshdfFilePw!')
    #end def

    def _write_h5_format(self,f):
        f.create_dataset('format',data=array(self.fmt,dtype=str_dtype))
    #end def

    def _write_h5_supercell(self,f):
        g = f.create_group('supercell')
        g.create_dataset('primitive_vectors',data=self.supercell)
    #end def
    
    def _write_h5_version(self,f):
        f.create_dataset('version',data=array(self.version,dtype=int_dtype))
    #end def

#end class


# Class for an ESHDF file in the PW representation
class EshdfFilePw(EshdfFileBase):

    # /electrons
    #  WIP
    density   = False
    rho_gv    = []
    rho_g     = []
    kpts      = []
    gvectors  = []
    numsym    = []
    eps       = []
    num_bands = [[]]
    psi_g     = [[[]]]
    symgroup  = []
    W_kpt     = []

    def _write_h5_electrons(self,f):
        g = f.create_group('electrons')
        if self.density:
            gd = g.create_group('density')
            gd.create_dataset('gvectors',data=self.rho_gv)
            gd.create_dataset('number_of_gvectors',data=array([self.rho_gv.shape[0]],dtype=int_dtype))
            for s in range(self.num_spin):
                gd0 = gd.create_group('spin_{}'.format(s))
                gd0.create_dataset('density_g',data=self.rho_g[s])
            #end for
        #end if
        for k,kpt in enumerate(self.kpts):
            kp0 = g.create_group('kpoint_{}'.format(k))
            if k==0:
                kp0.create_dataset('gvectors',data=self.gvectors)
                kp0.create_dataset('number_of_gvectors',data=array([self.gvectors.shape[0]],dtype=int_dtype))
            #end if
            kp0.create_dataset('numsym',data=array([self.numsym[k]],dtype=int_dtype))
            kp0.create_dataset('reduced_k',data=kpt)
            for s in range(self.num_spin):
                kp00 = kp0.create_group('spin_{}'.format(s))
                kp00.create_dataset('eigenvalues',data=self.eps[k][s])
                kp00.create_dataset('number_of_states',data=array([self.num_bands[k][s]],dtype=int_dtype))
                for n in range(self.num_bands[k][s]):
                    kp00s = kp00.create_group('state_{}'.format(n))
                    kp00s.create_dataset('psi_g',data=self.psi_g[k][s][n])
                #end for
            #end for
            kp0.create_dataset('symgroup',data=array([self.symgroup[k]],dtype=int_dtype))
            kp0.create_dataset('weight',data=array([self.W_kpt[k]],dtype=float_dtype))
        #end for
        g.create_dataset('number_of_electrons',data=self.num_elec)
        g.create_dataset('number_of_kpoints',data=array([len(self.kpts)],dtype=int_dtype))
        g.create_dataset('number_of_spins',data=array([self.num_spin],dtype=int_dtype))
    #end def

#end class

