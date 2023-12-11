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

# This file implements Eshdf bindings (Eshdf.py) for reading GPAW orbitals

from numpy import array,pi,dot,real,imag,concatenate,unique,zeros
from Eshdf import EshdfFilePw

# Calculate vector norm of 2-float representation of a complex vector
def calculate_norm(v):
    return ((v**2).sum())**0.5
#end def

# Helper function (obsolete) to match gvector g in an array of gvectors gset
# and return the value of the associated density coefficients from dset
def match_gv(g,gset,dset):
    g0 = gset[:,0]==g[0]
    g1 = gset[:,1]==g[1]
    g2 = gset[:,2]==g[2]
    return dset[logical_and(logical_and(g0,g1),g2),:][0]
#end def

# Add the negative half of gvectors removed by symmetry (full array expected by QMCPACK)
def desymmetrize_z(gv,psi):
    gv2  = -gv[1:].copy()
    psi2 = psi[1:].copy()
    psi2[:,1] *= -1
    gv2  = concatenate((gv, gv2), axis=0)
    psi2 = concatenate((psi,psi2),axis=0)
    return gv2,psi2
#end def

# Eshdf File class from GPAW
class EshdfFilePwGpaw(EshdfFilePw):

    # if infile and outfile provided, load automatically
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

    # initiate loading of data from the provided objects
    #   Compute and add charge density, if density=True
    def load(
        self,
        infile  = None,
        calc    = None,
        atoms   = None,
        density = False,
        **kwargs,
        ):
        self.density = density
        if infile is not None:
            if '.gpw' in infile:
                self.load_gpaw_restart(infile,**kwargs)
            else:
                print('Unsupported file {} not loaded.'.format(infile))
            #end if
        elif calc is not None and atoms is not None:
            self.load_gpaw_calc(calc,atoms)
        else:
            print('Could not load data for ESHDF.')
        #end if
    #end def

    # load data from GPAW restart file
    def load_gpaw_restart(
        self,
        infile,
        **kwargs,
        ):
        from gpaw import restart
        atoms,calc = restart(infile)
        try: # check if the wavefunction exists
            wfs = calc.wfs
        except:
            print('ERROR: No wavefunctions found in the GPAW restart file {}'.format(filename))
            print("Make sure to write the GPAW file with mode='all'!")
            raise NoWavefunctionFound
        #end try
        if self.density: 
            # run SCF to re-populate density
            print('Running SCF to re-populate density')
            calc.scf.converged = False
            calc.results.pop('energy')
            E = calc.get_potential_energy()
        #end def
        self.load_gpaw_calc(calc,atoms)
    #end def

    # GPAW: load data from calc and atoms objects
    #   Currently only supports PW mode
    def load_gpaw_calc(
        self,
        calc,
        atoms,
        ):
        mode = calc.wfs.mode
        if mode=='pw':
            self._load_application()
            self._load_atoms(calc,atoms)
            self._load_electrons(calc,atoms)
        else:
            print('Could not load data from GPAW Calculator in mode={}. Only pw mode is currently supported.'.format(mode))
        #end if
        self.loaded = True
    #end def

    # GPAW: load "application" data
    def _load_application(self):
        from gpaw import __version__ as gpaw_version
        from re import sub
        self.application = 'gpaw'
        self.app_version = sub('\D',' ',gpaw_version).split()[:3]
    #end def

    # GPAW: load "atoms" data
    def _load_atoms(self,calc,atoms):
        from ase.units import Bohr
        num_atoms = atoms.get_global_number_of_atoms()
        supercell = calc.wfs.gd.cell_cv
        syms      = atoms.get_chemical_symbols()
        masses    = atoms.get_masses()
        positions = atoms.get_positions()/Bohr # convert from A to Bohr
        Zs        = atoms.get_atomic_numbers()
        # figure out species data
        sym_set         = set(syms)
        num_species     = len(sym_set)
        species_ids     = []
        species_mass    = []
        species_Z       = []
        species_name    = []
        species_valence = []
        for s,name in enumerate(sym_set):
            i = syms.index(name)
            species_mass.append(masses[i])
            species_Z.append(Zs[i])
            species_valence.append(Zs[i]) # fixme for PP calculations
            species_name.append(name)
            for sym in syms:
                if sym==name:
                    species_ids.append(s)
                #end if
            #end for
        #end for
        self.species         = sym_set
        self.num_species     = num_species
        self.num_atoms       = num_atoms
        self.positions       = positions
        self.species_ids     = species_ids
        self.species_mass    = species_mass
        self.species_valence = species_valence
        self.species_Z       = species_Z
        self.species_name    = species_name
        self.supercell       = supercell
    #end def

    # GPAW: transform reciprocal vectors to Miller indices
    def _get_pw_indices(self,pd,k):
        return array(dot(pd.get_reciprocal_vectors(k,add_q=False),pd.gd.cell_cv/2/pi).round(0),dtype=int)
    #end def

    # GPAW: load "electrons" data
    def _load_electrons(self,calc,atoms):
        mesh     = calc.wfs.gd.N_c
        num_spin = calc.wfs.nspins
        N_tot    = calc.wfs.nvalence
        if num_spin==1:
            N_up = N_tot // 2
            N_dn = N_tot - N_up
        else: # spin-polarized # FIXME
            N_up = N_tot // 2
            N_dn = N_tot - N_up
        #end if
        # kpoints
        num_kpt  = calc.wfs.kd.get_count()
        W_kpt    = calc.wfs.kd.weight_k
        try:
            kpts    = calc.wfs.kd.get_bz_q_points()
        except:
            kpts    = calc.wfs.kd.bzk_kc
        #end try

        # Combine the GPAW plane wave data for each separate kpoint into one big array, gvall
        gvall = self._get_pw_indices(calc.wfs.pd,0)
        imin  = 0
        imax  = len(gvall)
        iis   = [[imin,imax]] # indices to mark the beginning and end of each kpoint
        for k in range(1,num_kpt):
            gv    = self._get_pw_indices(calc.wfs.pd,k)
            gvall = concatenate((gvall,gv),axis=0)
            imin  = imax
            imax += len(gv)
            iis.append([imin,imax])
        #end for
        # Gv is a new gvector basis with only unique vectors
        # indices help to associate the pw correct coefficients later on
        Gv,indices = unique(gvall,axis=0,return_inverse=True)

        data_numsym     = []
        data_symgroup   = []
        data_reduced_k  = []
        data_occupation = []
        data_num_bands  = []
        data_eps        = []
        data_psi_g      = []
        for k,kpt in enumerate(kpts):
            numsym = 1   # FIXME
            symgroup = 1 # FIXME
            data_numsym.append(numsym)
            data_symgroup.append(symgroup)
            ii = iis[k]
            # TODO: checks for spin polarized calculations?
            datas_eps       = []
            datas_psi_g     = []
            datas_num_bands = []
            for s in range(num_spin):
                occ       = calc.wfs.kpt_u[k].f_n
                eps       = calc.wfs.kpt_u[k].eps_n
                num_bands = len(eps)
                u         = k*num_spin+s
                datasn_psi_g = []
                for n in range(num_bands):
                    # get wavefunction in PW
                    wf_c  = calc.wfs.kpt_u[u].psit_nG[n]
                    wf    = array([real(wf_c),imag(wf_c)])
                    psi   = zeros((len(Gv),2))
                    psi[indices[ii[0]:ii[1]]] = wf.T/calculate_norm(wf) # normalize the crude way
                    datasn_psi_g.append(psi)
                #end for
                datas_num_bands.append(num_bands)
                datas_eps.append(eps)
                datas_psi_g.append(datasn_psi_g)
            #end for
            data_eps.append(datas_eps)
            data_num_bands.append(datas_num_bands)
            data_psi_g.append(datas_psi_g)
        #end for
        self.gvectors   = Gv
        self.num_elec   = [N_up,N_dn]
        self.kpts       = kpts
        self.numsym     = data_numsym
        self.symgroup   = data_symgroup
        self.reduced_k  = data_reduced_k
        self.W_kpt      = W_kpt
        self.num_bands  = data_num_bands
        self.num_spin   = num_spin
        self.eps        = data_eps
        self.psi_g      = data_psi_g
        if self.density:
            # FIXME: spin polarized
            # read data
            rho_gv   = self._get_pw_indices(calc.density.pd2,0)
            rho_psic = calc.density.nt_Q
            # process
            rho_psi  = array([real(rho_psic),imag(rho_psic)]).T
            gv,rho   = desymmetrize_z(rho_gv,rho_psi)
            n        = N_tot/rho[0,0]/calc.density.gd.volume # normalize at the gamma point
            # save
            self.rho_gv   = gv
            self.rho_g    = rho*n
        #end if
    #end def

    # file writing functions are defined in the base classes

#end class
