##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  physical_system.py                                                #
#    Representations of particles collected together in complete     #
#    systems.                                                        #
#                                                                    #
#  Content summary:                                                  #
#    PhysicalSystem                                                  #
#      Class representing electrons+ions for a simulation.           #
#                                                                    #
#    generate_physical_system                                        #
#      User function to create arbitrary physical systems.           #
#                                                                    #
#====================================================================#

import os
from pathlib import Path
from copy import deepcopy
import numpy as np
from .developer import DevBase, obj, error, warn
from .periodic_table import Elements
from .structure import Structure, generate_structure, read_structure


class PhysicalSystem(DevBase):

    ghost_aliases = ["Xx"]

    def __init__(self,structure=None,net_charge=0,net_spin=0,**valency):
        self.pseudized = False
        if structure is None:
            self.structure = Structure()
        else:
            self.structure = structure
        #end if

        self.folded_system = None
        if self.structure.has_folded():
            if self.structure.is_tiled():
                vratio = structure.volume()/structure.folded_structure.volume()
                ncells = int(round(vratio))
                if abs(vratio-ncells)>1e-4:
                    self.error('volume of system does not divide evenly into folded system')
                #end if
                if net_charge%ncells!=0:
                    self.error('net charge of system does not divide evenly into folded system')
                #end if
                if isinstance(net_spin,str):
                    net_spin_fold = net_spin
                elif net_spin%ncells!=0:
                    self.error('net_spin of system does not divide evenly into folded system')
                else:
                    net_spin_fold = net_spin//ncells 
                #end if
                net_charge_fold = net_charge//ncells
            elif not self.structure.has_axes(): # folded molecule
                # net charge/spin are not physically meaningful
                # for a point group folded molecule
                # set them to safe values; they will not be used later
                net_charge_fold = 0
                net_spin_fold   = 'low'
            else:
                self.error('folded structure is not correctly integrated with full structure\nfolded physical system cannot be constructed')
            #end if
                
            self.folded_system = PhysicalSystem(
                structure  = structure.folded_structure,
                net_charge = net_charge_fold,
                net_spin   = net_spin_fold,
                **valency
                )
        #end if
        if valency is not None and len(valency) > 0:
            self.pseudize(**valency)
        else:
            self.valency = None
        self.net_charge = net_charge
        self.net_spin   = net_spin

        self.check_folded_system()
    #end def __init__


    def pseudize(self,**valency):
        for ion in valency.keys():
            if ion not in self.ion_labels:
                self.error(ion+' is not in the physical system',exit=False)

        self.valency = obj(**valency)
        self.pseudized = True
    #end def pseudize

        
    def check_folded_system(self,*,exit=True,message=False):
        msg = ''
        sys_folded    = self.folded_system is not None
        struct_folded = self.structure.folded_structure is not None
        if sys_folded!=struct_folded:
            msg+='folding of physical system and structure is not consistent\nsystem folded: {0}\nstructure folded: {1}\n'.format(sys_folded,struct_folded)
        #end if
        if sys_folded and id(self.structure.folded_structure)!=id(self.folded_system.structure):
            msg+='structure of folded system and folded structure are distinct\nthis is not allowed and may be a developer error'
        #end if
        success = len(msg)==0
        if not success and exit:
            self.error(msg)
        #end if
        if not message:
            return success
        else:
            return success,msg
        #end if
    #end def check_folded_system


    def check_consistent(self,tol=1e-8,*,exit=True,message=False):
        fs,fm = self.check_folded_system(exit=False,message=True)
        cs,cm = self.structure.check_consistent(tol,exit=False,message=True)
        msg = ''
        if not fs:
            msg += fm+'\n'
        #end if
        if not cs:
            msg += cm+'\n'
        #end if
        consistent = len(msg)==0
        if not consistent and exit:
            self.error(msg)
        #end if
        if not message:
            return consistent
        else:
            return consistent,msg
        #end if
    #end def check_consistent


    def is_valid(self):
        return self.check_consistent(exit=False)
    #end def is_valid


    def change_units(self,units):
        self.structure.change_units(units,folded=False)
        if self.folded_system is not None:
            self.folded_system.change_units(units)
        #end if
    #end def change_units


    def group_atoms(self):
        self.structure.group_atoms(folded=False)
        if self.folded_system is not None:
            self.folded_system.group_atoms()
        #end if
    #end def group_atoms


    def rename(self,*,folded=True,**name_pairs):
        self.structure.rename(folded=False,**name_pairs)
        if self.pseudized:
            for old,new in name_pairs.items():
                if old in self.valency:
                    if new not in self.valency:
                        self.valency[new] = self.valency[old]
                    #end if
                    del self.valency[old]
                #end if
            #end for
        #end if
        if self.folded_system is not None and folded:
            self.folded_system.rename(folded=folded,**name_pairs)
        #end if
    #end def rename


    def copy(self):
        cp = deepcopy(self)
        if self.folded_system is not None and self.structure.folded_structure is not None:
            del cp.folded_system.structure
            cp.folded_system.structure = cp.structure.folded_structure
        #end if
        return cp
    #end def copy


    def load(self,filepath):
        DevBase.load(self,filepath)
        if self.folded_system is not None and self.structure.folded_structure is not None:
            del self.folded_system.structure
            self.folded_system.structure = self.structure.folded_structure
        #end if
    #end def load


    def tile(self,*td,**kwargs):
        extensive = True
        net_spin  = None
        if 'extensive' in kwargs:
            extensive = kwargs['extensive']
        #end if
        if 'net_spin' in kwargs:
            net_spin = kwargs['net_spin']
        #end if
        supercell = self.structure.tile(*td)
        supercell.remove_folded()
        if extensive:
            ncells = int(round(supercell.volume()/self.structure.volume()))
            net_charge = ncells*self.net_charge
            if net_spin is None:
                net_spin = ncells*self.net_spin
            #end if
        else:
            net_charge = self.net_charge
            if net_spin is None:
                net_spin   = self.net_spin
            #end if
        #end if
        system = deepcopy(self)
        supersystem = PhysicalSystem(
            structure  = supercell,
            net_charge = net_charge,
            net_spin   = net_spin,
            **self.valency
            )
        supersystem.folded_system = system
        supersystem.structure.set_folded(system.structure)
        return supersystem
    #end def tile


    def has_folded(self):
        return self.folded_system is not None
    #end def has_folded


    def remove_folded_system(self):
        self.folded_system = None
        self.structure.remove_folded_structure()
    #end def remove_folded_system


    def remove_folded(self):
        self.remove_folded_system()
    #end def remove_folded


    def get_smallest(self):
        if self.has_folded():
            return self.folded_system
        else:
            return self
        #end if
    #end def get_smallest

    
    def is_magnetic(self):
        return self.net_spin!=0 or self.structure.is_magnetic()
    #end def is_magnetic

    
    def spin_polarized_orbitals(self):
        return self.is_magnetic()
    #end def spin_polarized_orbitals


    # test needed
    def large_Zeff_elem(self,Zmin):
        elem = []
        for atom,Zeff in self.valency.items():
            if Zeff>Zmin:
                elem.append(atom)
            #end if
        #end for
        return elem
    #end def large_Zeff_elem


    # test needed
    def ae_pp_species(self):
        species = set(self.structure.elem)
        if self.pseudized:
            pp_species = set(self.valency.keys())
            ae_species = species-pp_species
        else:
            pp_species = set()
            ae_species = species
        #end if
        return ae_species,pp_species
    #end def ae_pp_species


    def kf_rpa(self):
      nelecs = (self.n_up, self.n_down)
      volume = self.structure.volume()
      kvol1 = (2*np.pi)**3/volume  # k-space volume per particle
      kfs = [(3*nelec*kvol1/(4*np.pi))**(1./3) for nelec in nelecs]
      return np.array(kfs)
    #end def kf_rpa


    @property
    def n_elec(self):
        ions = self.structure.elem.tolist()
        tot_charge = 0
        for ion in ions:
            if self.valency is not None:
                if ion in self.valency:
                    tot_charge += self.valency[ion]
                else:
                    _, element = Elements.is_element(ion, return_element=True)
                    tot_charge += element.atomic_number
            else:
                _, element = Elements.is_element(ion, return_element=True)
                tot_charge += element.atomic_number

        return tot_charge - self.net_charge

    @property
    def n_up(self):
        return (self.n_elec + self.net_spin) // 2

    @property
    def n_down(self):
        return (self.n_elec - self.net_spin) // 2

    @property
    def n_species(self):
        return len(set(self.structure.elem))

    @property
    def n_ions(self):
        return len(self.structure.elem)

    @property
    def ion_labels(self):
        return set(self.structure.elem)

    @property
    def Zeff(self):
        if self.valency is not None:
            return self.valency

        Zeff = {}
        for ion in self.ion_labels:
            _, element = Elements.is_element(ion, return_element=True)
            Zeff[ion] = element.atomic_number

        return Zeff
#end class PhysicalSystem


ps_defaults = dict(
    type='crystal',
    kshift = (0,0,0),
    net_charge=0,
    net_spin=0,
    pretile=None,
    tiling=None,
    tiled_spin=None,
    extensive=True
    )
def generate_physical_system(**kwargs):
    for var,val in ps_defaults.items():
        if var not in kwargs:
            kwargs[var] = val
        #end if
    #end for
    type = kwargs['type']
    if type=='atom' or type=='dimer' or type=='trimer':
        del kwargs['kshift']
        del kwargs['tiling']
        #if not 'units' in kwargs:
        #    kwargs['units'] = 'B'
        ##end if
        tiling = None
    else:
        tiling = kwargs['tiling']
    #end if

    if 'structure' in kwargs:
        s = kwargs['structure']
        is_str = isinstance(s,str)
        if is_str or isinstance(s, Path):
            if os.path.exists(s):
                if 'elem' in kwargs:
                    s = read_structure(s,elem=kwargs['elem'])
                else:
                    s = read_structure(s)
                #end if
                if 'axes' in kwargs:
                    s.reset_axes(kwargs['axes'])
                #end if
                kwargs['structure'] = s
            else:
                slow = s.lower()
                format = None
                if '.' in slow:
                    format = slow.rsplit('.')[1]
                elif 'poscar' in slow:
                    format = 'poscar'
                #end if
                is_path = '/' in s
                is_file = format in set('xyz xsf poscar cif fhi-aims'.split())
                if is_path or is_file:
                    error('user provided structure file does not exist\nstructure file path: '+s,'generate_physical_system')
                #end if
            #end if
        #end if
    #end if

    generation_info = obj(**deepcopy(kwargs))

    net_charge = kwargs['net_charge']
    net_spin   = kwargs['net_spin']
    tiled_spin = kwargs['tiled_spin']
    extensive  = kwargs['extensive']
    del kwargs['net_spin']
    del kwargs['net_charge']
    del kwargs['tiled_spin']
    del kwargs['extensive']
    if 'particles' in kwargs:
        warn("The keyword `particles` is no longer valid. Please remove from your scripts!")
        del kwargs['particles']

    pretile = kwargs['pretile']
    del kwargs['pretile']
    valency = dict()
    remove = []
    for var in kwargs:
        #if var in Matter.elements:
        if Elements.is_element(var):
            valency[var] = kwargs[var]
            remove.append(var)
        #end if
    #end if
    generation_info.valency = deepcopy(valency)
    for var in remove:
        del kwargs[var]
    #end for

    if pretile is None:
        structure = generate_structure(**kwargs)
    else:
        for d in range(len(pretile)):
            if tiling[d]%pretile[d]!=0:
                error('pretile does not divide evenly into tiling\n  tiling provided: {0}\n  pretile provided: {1}'.format(tiling,pretile),'generate_physical_system')
            #end if
        #end for
        tiling = tuple(np.array(tiling)//np.array(pretile))
        kwargs['tiling'] = pretile
        pre = generate_structure(**kwargs)
        pre.remove_folded_structure()
        structure = pre.tile(tiling)
    #end if
    if isinstance(tiling, tuple):
        tiling_mat = np.diag(tiling)
    elif tiling is None:
        tiling_mat = np.eye(3)
    else:
        tiling_mat = tiling

    if not np.array_equal(tiling_mat, np.eye(3)) and structure.has_folded():
        # Has some supercell tiling
        fps = PhysicalSystem(
            structure  = structure.folded_structure,
            net_charge = net_charge,
            net_spin   = net_spin,
            **valency
            )
        structure.remove_folded()
        folded_structure = fps.structure
        if extensive:
            ncells = int(round(structure.volume()/folded_structure.volume()))
            net_charge = ncells*net_charge
            if not isinstance(net_spin,str):
                net_spin   = ncells*net_spin
            #end if
        #end if
        if tiled_spin is not None:
            net_spin = tiled_spin
        #end if
        ps = PhysicalSystem(
            structure  = structure,
            net_charge = net_charge,
            net_spin   = net_spin,
            **valency
            )
        structure.set_folded(folded_structure)
        ps.folded_system = fps
    else:
        # No supercell tiling
        ps = PhysicalSystem(
            structure  = structure,
            net_charge = net_charge,
            net_spin   = net_spin,
            **valency
            )
    #end if
    
    ps.generation_info = generation_info

    return ps
#end def generate_physical_system



# test needed
def ghost_atoms(*particles):
    for particle in particles:
        PhysicalSystem.ghost_aliases.append(particle)
#end def ghost_atoms
