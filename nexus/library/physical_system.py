##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  physical_system.py                                                #
#    Representations of matter, particles, and particles collected   #
#    together in complete systems.                                   #
#                                                                    #
#  Content summary:                                                  #
#    PhysicalSystem                                                  #
#      Class representing electrons+ions for a simulation.           #
#                                                                    #
#    generate_physical_system                                        #
#      User function to create arbitrary physical systems.           #
#                                                                    #
#    Matter                                                          #
#      Base class for all forms of matter.                           #
#      Class contains a list of all matter known to Nexus.           #
#                                                                    #
#    Particle                                                        #
#      Class representing a particular particle species.             #
#                                                                    #
#    Ion, PseudoIon                                                  #
#      Specialized Particle classes for ion species and ions         #
#      represented by pseudopotentials.                              #
#                                                                    #
#    Particles                                                       #
#      A collection of particles.                                    #
#      PhysicalSystem objects contain a Particles instance.          #
#                                                                    #
#====================================================================#


from numpy import dot,array
from numpy.linalg import inv
from generic import obj
from developer import DevBase
from unit_converter import convert
from periodic_table import is_element
from structure import Structure
from debug import *


class Matter(DevBase):
    particle_collection = None

    @classmethod
    def set_elements(cls,elements):
        cls.elements = set(elements)
    #end def set_elements

    @classmethod
    def set_particle_collection(cls,pc):
        cls.particle_collection = pc
    #end def set_particle_collection

    @classmethod
    def new_particles(cls,*particles,**named_particles):
        cls.particle_collection.add_particles(*particles,**named_particles)
    #end def new_particles

    def is_element(self,name,symbol=False):
        return is_element(name,symbol=symbol)
    #end def is_element
#end class Matter


class Particle(Matter):
    def __init__(self,name=None,mass=None,charge=None,spin=None):
        self.name   = name  
        self.mass   = mass  
        self.charge = charge
        self.spin   = spin  
    #end def __init__

    def set_count(self,count):
        self.count  = count 
    #end def set_count
#end class Particle


class Ion(Particle):
    def __init__(self,name=None,mass=None,charge=None,spin=None,
                 protons=None,neutrons=None):
        Particle.__init__(self,name,mass,charge,spin)
        self.protons  = protons
        self.neutrons = neutrons
    #end def __init__

    def pseudize(self,valence):
        ps = PseudoIon()
        ps.transfer_from(self)
        ps.charge = valence
        ps.core_electrons    = ps.protons - valence
        return ps
    #end def pseudize
#end class


class PseudoIon(Ion):
    def __init__(self,name=None,mass=None,charge=None,spin=None,
                 protons=None,neutrons=None,core_electrons=None):
        Ion.__init__(self,name,mass,charge,spin,protons,neutrons)
        self.core_electrons    = core_electrons    
    #end def __init__
#end class PseudoIon


class Particles(Matter):
    def __init__(self,*particles,**named_particles):
        self.add_particles(*particles,**named_particles)
    #end def __init__

    def add_particles(self,*particles,**named_particles):
        if len(particles)==1 and isinstance(particles[0],list):
            particles = particles[0]
        #end if
        for particle in particles:
            self[particle.name] = particle
        #end for
        for name,particle in named_particles.iteritems():
            self[name] = particle
        #end for
    #end def add_particles

    def get_particle(self,name):
        p = None
        if name in self:
            p = self[name]
        else:
            iselem,symbol = self.is_element(name,symbol=True)
            if iselem and symbol in self:
                p = self[symbol].copy()
                p.name = name
                self[name] = p
            #end if
        #end if
        return p
    #end def get_particle

    def get(self,quantity):
        q = obj()
        for name,particle in self.iteritems():
            q[name] = particle[quantity]
        #end for
        return q
    #end def get

    def rename(self,**name_pairs):
        for old,new in name_pairs.iteritems():
            if old in self:
                o = self[old]
                o.name = new
                del self[old]
                if new in self:
                    self[new].count += o.count
                else:
                    self[new] = o
                #end if
            #end if
        #end for
    #end def rename

    def get_ions(self):
        ions = obj()
        for name,particle in self.iteritems():
            if self.is_element(name):
                ions[name] = particle
            #end if
        #end for
        return ions
    #end def get_ions
        
    def count_ions(self,species=False):
        nions = 0
        nspecies = 0
        for name,particle in self.iteritems():
            if self.is_element(name):
                nspecies += 1
                nions += particle.count
            #end if
        #end for
        if species:
            return nions,nspecies
        else:
            return nions
        #end if
    #end def count_ions

    def get_electrons(self):
        electrons = obj()
        for electron in ('up_electron','down_electron'):
            if electron in self:
                electrons[electron] = self[electron]
            #end if
        #end for
        return electrons
    #end def get_electrons

    def count_electrons(self):
        nelectrons = 0
        for electron in ('up_electron','down_electron'):
            if electron in self:
                nelectrons += self[electron].count
            #end if
        #end for
        return nelectrons
    #end def count_electrons

    def electron_counts(self):
        counts = []
        for electron in ('up_electron','down_electron'):
            if electron in self:
                counts.append(self[electron].count)
            else:
                counts.append(0)
            #end if
        #end for
        return counts
    #end def electron_counts
#end class Particles

me_amu = convert(1.,'me','amu')

plist = [
    Particle('up_electron'  ,1.0,-1, 1),
    Particle('down_electron',1.0,-1,-1),
    ]
from periodic_table import ptable
for name,a in ptable.elements.iteritems():
    spin = 0 # don't have this data
    protons  = a.atomic_number
    neutrons = int(round(a.atomic_weight['amu']-a.atomic_number))
    p = Ion(a.symbol,a.atomic_weight['me'],a.atomic_number,spin,protons,neutrons)
    plist.append(p)
#end for
for name,iso in ptable.isotopes.iteritems():
    for mass_number,a in iso.iteritems():
        spin = 0 # don't have this data
        protons  = a.atomic_number
        neutrons = int(round(a.atomic_weight['amu']-a.atomic_number))
        p = Ion(a.symbol+'_'+str(mass_number),a.atomic_weight['me'],a.atomic_number,spin,protons,neutrons)
        plist.append(p)
    #end for
#end for

Matter.set_elements(ptable.elements.keys())
Matter.set_particle_collection(Particles(plist))

del plist



class PhysicalSystem(Matter):

    def __init__(self,structure=None,net_charge=0,net_spin=0,particles=None,**valency):
        self.pseudized = False

        if structure is None:
            self.structure = Structure()
        else:
            self.structure = structure
        #end if
        if particles is None:
            self.particles = Particles()
        else:
            self.particles = particles.copy()
        #end if

        self.folded_system = None
        if self.structure.folded_structure!=None:
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
                net_spin_fold = net_spin/ncells 
            #end if
            self.folded_system = PhysicalSystem(
                structure  = structure.folded_structure,
                net_charge = net_charge/ncells,
                net_spin   = net_spin_fold,
                particles  = particles,
                **valency
                )
        #end if

        self.valency_in = obj(**valency)
        self.net_charge_in = net_charge
        self.net_spin_in   = net_spin

        self.update_particles(clear=False)

        self.check_folded_system()
    #end def __init__


    def update_particles(self,clear=True):
        #add ions
        pc = dict()
        elem = list(self.structure.elem)
        for ion in set(elem):
            pc[ion] = elem.count(ion)
        #end for
        missing = set(pc.keys())-set(self.particles.keys())
        if len(missing)>0 or len(elem)==0:
            if clear:
                self.particles.clear()
            #end if
            self.add_particles(**pc)

            #pseudize
            if len(self.valency_in)>0:
                self.pseudize(**self.valency_in)
            #end if

            #add electrons
            self.generate_electrons(self.net_charge_in,self.net_spin_in)
        #end if
    #end def update_particles


    def update(self):
        self.net_charge = self.structure.background_charge
        self.net_spin   = 0
        for p in self.particles:
            self.net_charge += p.count*p.charge
            self.net_spin   += p.count*p.spin
        #end for
        self.net_charge = int(round(float(self.net_charge)))
        self.net_spin   = int(round(float(self.net_spin)))
    #end def update


    def add_particles(self,**particle_counts):
        pc = self.particle_collection # all known particles
        plist = []
        for name,count in particle_counts.iteritems():
            particle = pc.get_particle(name)
            if particle is None:
                self.error('particle {0} is unknown'.format(name))
            else:
                particle = particle.copy()
            #end if
            particle.set_count(count)
            plist.append(particle)
        #end for
        self.particles.add_particles(plist)
        self.update()
    #end def add_particles


    def generate_electrons(self,net_charge=0,net_spin=0):
        nelectrons = -net_charge + self.net_charge
        if net_spin=='low':
            net_spin = nelectrons%2
        #end if
        nup   = float(nelectrons + net_spin - self.net_spin)/2
        ndown = float(nelectrons - net_spin + self.net_spin)/2        
        if abs(nup-int(nup))>1e-3:
            self.error('requested spin state {0} incompatible with {1} electrons'.format(net_spin,nelectrons))
        #end if
        nup   = int(nup)
        ndown = int(ndown)
        self.add_particles(up_electron=nup,down_electron=ndown)
    #end def generate_electrons


    def pseudize(self,**valency):
        errors = False
        for ion,valence_charge in valency.iteritems():
            if ion in self.particles:
                ionp = self.particles[ion]
                if isinstance(ionp,Ion):
                    self.particles[ion] = ionp.pseudize(valence_charge)
                    self.pseudized = True
                else:
                    self.error(ion+' cannot be pseudized',exit=False)
                #end if
            else:
                self.error(ion+' is not in the physical system',exit=False)
                errors = True
            #end if
        #end for
        if errors:
            self.error('system cannot be generated')
        #end if
        self.valency = obj(**valency)
        self.update()
    #end def pseudize

        
    def check_folded_system(self,exit=True,message=False):
        msg = ''
        sys_folded    = self.folded_system!=None
        struct_folded = self.structure.folded_structure!=None
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


    def check_consistent(self,tol=1e-8,exit=True,message=False):
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


    def change_units(self,units):
        self.structure.change_units(units,folded=False)
        if self.folded_system!=None:
            self.folded_system.change_units(units)
        #end if
    #end def change_units


    def group_atoms(self):
        self.structure.group_atoms(folded=False)
        if self.folded_system!=None:
            self.folded_system.group_atoms()
        #end if
    #end def group_atoms


    def rename(self,folded=True,**name_pairs):
        self.particles.rename(**name_pairs)
        self.structure.rename(folded=False,**name_pairs)
        if self.pseudized:
            for old,new in name_pairs.iteritems():
                if old in self.valency:
                    if new not in self.valency:
                        self.valency[new] = self.valency[old]
                    #end if
                    del self.valency[old]
                #end if
            #end for
            self.valency_in = self.valency
        #end if
        if self.folded_system!=None and folded:
            self.folded_system.rename(folded=folded,**name_pairs)
        #end if
    #end def rename


    def copy(self):
        cp = DevBase.copy(self)
        if self.folded_system!=None and self.structure.folded_structure!=None:
            del cp.folded_system.structure
            cp.folded_system.structure = cp.structure.folded_structure
        #end if
        return cp
    #end def copy


    def load(self,filepath):
        DevBase.load(self,filepath)
        if self.folded_system!=None and self.structure.folded_structure!=None:
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
                net_spin   = ncells*self.net_spin
            #end if
        else:
            net_charge = self.net_charge
            if net_spin is None:
                net_spin   = self.net_spin
            #end if
        #end if
        system = self.copy()
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


    def remove_folded_system(self):
        self.folded_system = None
        self.structure.remove_folded_structure()
    #end def remove_folded_system


    def remove_folded(self):
        self.remove_folded_system()
    #end def remove_folded


    def get_primitive(self):
        if self.folded_system is None:
            fs = self
        else:
            fs = self.folded_system
            while fs.folded_system!=None:
                fs = fs.folded_system
            #end while
        #end if
        return fs
    #end def get_primitive


    def folded_representation(self,arg0,arg1=None):
        self.error('folded_representation needs a developers attention to make it equivalent with tile')
        if isinstance(arg0,PhysicalSystem):
            folded_system = arg0
        elif isinstance(arg0,str):
            shape = arg0
            tiling    = arg1
            if tiling is None:
                tiling = (1,1,1)
            #end if
            if not 'generation_info' in self:
                self.error('system was not formed with generate_physical_system, cannot form folded representation')
            #end if
            structure,element,scale,units,net_charge,net_spin,particles,valency = \
                self.generation_info.tuple('structure','element','scale','units', \
                                               'net_charge','net_spin','particles','valency')
            folded_system = generate_physical_system(
                structure  = structure,
                shape      = shape,
                element    = element,
                tiling     = tiling,
                scale      = scale,
                units      = units,
                net_charge = net_charge,
                net_spin   = net_spin,
                particles  = particles,
                **valency
                )
        else:
            self.error('unrecognized inputs in folded_representation')
        #end if
        tilematrix,kmap = self.structure.fold(folded_system.structure,'tilematrix','kmap')
        self.set(
            folded_system = folded_system,
            tilematrix    = tilematrix,
            kmap          = kmap
            )
        return folded_system
    #end def folded_representation


    def large_Zeff_elem(self,Zmin):
        elem = []
        for atom,Zeff in self.valency.iteritems():
            if Zeff>Zmin:
                elem.append(atom)
            #end if
        #end for
        return elem
    #end def large_Zeff_elem


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
#end class PhysicalSystem


import os
from structure import generate_structure,read_structure
from copy import deepcopy
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
    for var,val in ps_defaults.iteritems():
        if not var in kwargs:
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
        if is_str:
            if os.path.exists(s):
                if 'elem' in kwargs:
                    kwargs['structure'] = read_structure(s,elem=kwargs['elem'])
                else:
                    kwargs['structure'] = read_structure(s)
                #end if
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
                    PhysicalSystem.class_error('user provided structure file does not exist\nstructure file path: '+s,'generate_physical_system')
                #end if
            #end if
        #end if
    #end if

    generation_info = obj()
    generation_info.transfer_from(deepcopy(kwargs))

    net_charge = kwargs['net_charge']
    net_spin   = kwargs['net_spin']
    tiled_spin = kwargs['tiled_spin']
    extensive  = kwargs['extensive']
    del kwargs['net_spin']
    del kwargs['net_charge']
    del kwargs['tiled_spin']
    del kwargs['extensive']
    if 'particles' in kwargs:
        particles = kwargs['particles']
        del kwargs['particles']
    else:
        generation_info.particles = None
    #end if
    pretile = kwargs['pretile']
    del kwargs['pretile']
    valency = dict()
    remove = []
    for var in kwargs:
        #if var in Matter.elements:
        if is_element(var):     
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
                PhysicalSystem.class_error('pretile does not divide evenly into tiling\n  tiling provided: {0}\n  pretile provided: {1}'.format(tiling,pretile),'generate_physical_system')
            #end if
        #end for
        tiling = tuple(array(tiling)/array(pretile))
        kwargs['tiling'] = pretile
        pre = generate_structure(**kwargs)
        pre.remove_folded_structure()
        structure = pre.tile(tiling)
    #end if

    if tiling!=None and tiling!=(1,1,1) and structure.has_folded():
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
        if tiled_spin!=None:
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




if __name__=='__main__':
    
    from structure import generate_structure

    ps = PhysicalSystem(
        structure = generate_structure('diamond','sc','Ge',(2,2,2),scale=5.639,units='A'),
        net_charge = 1,
        net_spin   = 1,
        Ge = 4
        )

    print 'net_charge',ps.net_charge
    print 'net_spin  ',ps.net_spin
    print ps.particles

#end if
