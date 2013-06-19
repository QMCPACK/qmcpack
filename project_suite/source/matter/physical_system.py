
from numpy import dot
from numpy.linalg import inv
from generic import obj
from developer import DevBase
from unit_converter import convert
from structure import Structure

class Matter(DevBase):
    particle_collection = None

    @classmethod
    def new_particles(*particles,**named_particles):
        self.particle_collection.add_particles(*particles,**named_particles)
    #end def new_particles
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

    def get(self,quantity):
        q = obj()
        for name,particle in self.iteritems():
            q[name] = particle[quantity]
        #end for
        return q
    #end def get

    def get_ions(self):
        ions = obj()
        for name,particle in self.iteritems():
            if name in self.elements:
                ions[name] = particle
            #end if
        #end for
        return ions
    #end def get_ions
        
    def count_ions(self):
        nions = 0
        nspecies = 0
        for name,particle in self.iteritems():
            if name in self.elements:
                nspecies += 1
                nions += particle.count
            #end if
        #end for
        return nions,nspecies
    #end def count_ions

    def get_electrons(self):
        electrons = obj()
        for electron in ['up_electron','down_electron']:
            if electron in self:
                electrons[electron] = self[electron]
            #end if
        #end for
        return electrons
    #end def get_electrons

    def count_electrons(self):
        nelectrons = 0
        for electron in ['up_electron','down_electron']:
            if electron in self:
                nelectrons += self[electron].count
            #end if
        #end for
        return nelectrons
    #end def count_electrons
#end class Particles

me_amu = convert(1.,'me','amu')

plist = [
    Particle('up_electron'  ,1.0,-1, 1),
    Particle('down_electron',1.0,-1,-1),
    ]
from periodic_table import PeriodicTable
pt = PeriodicTable()
for name,a in pt.elements.iteritems():
    spin = 0 # don't have this data
    protons  = a.atomic_number
    neutrons = int(round(a.atomic_weight['amu']-a.atomic_number))
    p = Ion(a.symbol,a.atomic_weight['me'],a.atomic_number,spin,protons,neutrons)
    plist.append(p)
#end for
for name,iso in pt.isotopes.iteritems():
    for mass_number,a in iso.iteritems():
        spin = 0 # don't have this data
        protons  = a.atomic_number
        neutrons = int(round(a.atomic_weight['amu']-a.atomic_number))
        p = Ion(a.symbol+'_'+str(mass_number),a.atomic_weight['me'],a.atomic_number,spin,protons,neutrons)
        plist.append(p)
    #end for
#end for

Matter.elements = pt.elements.keys()
Matter.particle_collection = Particles(plist)

del plist
del pt



class PhysicalSystem(Matter):

    def __init__(self,structure=None,net_charge=0,net_spin=0,particles=None,**valency):

        self.structure = structure
        self.particles = particles
        if structure==None:
            self.structure = Structure()
        #end if
        if particles==None:
            self.particles = Particles()
        #end if

        self.folded_system = None
        if self.structure.folded_structure!=None:
            self.folded_system = PhysicalSystem(
                structure  = structure.folded_structure,
                net_charge = net_charge,
                net_spin   = net_spin,
                particles  = particles,
                **valency
                )
        #end if

        #add ions
        pc = dict()
        elem = list(self.structure.elem)
        for ion in set(elem):
            pc[ion] = elem.count(ion)
        #end for
        self.add_particles(**pc)

        #pseudize
        if len(valency)>0:
            self.pseudize(**valency)
        #end if

        #add electrons
        self.generate_electrons(net_charge,net_spin)

        self.check_folded_system()
    #end def __init__


    def update(self):
        self.net_charge = 0
        self.net_spin   = 0
        for p in self.particles:
            self.net_charge += p.count*p.charge
            self.net_spin   += p.count*p.spin
        #end for
        self.net_charge = int(round(float(self.net_charge)))
        self.net_spin = int(round(float(self.net_spin)))
    #end def update


    def add_particles(self,**particle_counts):
        plist = []
        for name,count in particle_counts.iteritems():
            particle = self.particle_collection[name].copy()
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

        
    def check_folded_system(self):
        sys_folded    = self.folded_system!=None
        struct_folded = self.structure.folded_structure!=None
        if sys_folded!=struct_folded:
            self.error('folding of physical system and structure is not consistent\n  system folded: {0}\n  structure folded: {1}'.format(sys_folded,struct_folded))
        #end if
        if sys_folded and id(self.structure.folded_structure)!=id(self.folded_system.structure):
            self.error('structure of folded system and folded structure are distinct\n  this is not allowed and may be a developer error')
        #end if
    #end def check_folded_system


    def change_units(self,units):
        self.structure.change_units(units)
        if self.folded_system!=None:
            self.folded_system.change_units(units)
        #end if
    #end def change_units


    def group_atoms(self):
        self.structure.group_atoms()
        if self.folded_system!=None:
            self.folded_system.group_atoms()
        #end if
    #end def group_atoms


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
        if 'extensive' in kwargs:
            extensive = kwargs['extensive']
        #end if
        supercell = self.structure.tile(*td)
        if extensive:
            ncells = int(round(supercell.volume()/self.structure.volume()))
            net_charge = ncells*self.net_charge
            net_spin   = ncells*self.net_spin
        else:
            net_charge = self.net_charge
            net_spin   = self.net_spin
        #end if
        supersystem = PhysicalSystem(
            structure  = supercell,
            net_charge = net_charge,
            net_spin   = net_spin,
            **self.valency
            )
        return supersystem
    #end def tile


    def remove_folded_system(self):
        print 'removing folded system',self.folded_system.__class__.__name__
        self.folded_system = None
        self.structure.remove_folded_structure()
    #end def remove_folded_system


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
#end class PhysicalSystem



from structure import generate_structure
from copy import deepcopy
ps_defaults = dict(
    type='crystal',
    kshift = (0,0,0),
    net_charge=0,
    net_spin=0
    )
def generate_physical_system(**kwargs):

    for var,val in ps_defaults.iteritems():
        if not var in kwargs:
            kwargs[var] = val
        #end if
    #end for
    if 'type' in kwargs and kwargs['type']=='atom':
        del kwargs['kshift']
        if not 'units' in kwargs:
            kwargs['units'] = 'B'
        #end if
    #end if

    generation_info = obj()
    generation_info.transfer_from(deepcopy(kwargs))


    net_charge = kwargs['net_charge']
    net_spin = kwargs['net_spin']
    del kwargs['net_spin']
    del kwargs['net_charge']
    if 'particles' in kwargs:
        particles = kwargs['particles']
        del kwargs['particles']
    else:
        generation_info.particles = None
    #end if
    valency = dict()
    remove = []
    for var in kwargs:
        if var in Matter.elements:
            valency[var] = kwargs[var]
            remove.append(var)
        #end if
    #end if
    generation_info.valency = deepcopy(valency)
    for var in remove:
        del kwargs[var]
    #end for

    ps = PhysicalSystem(
        structure  = generate_structure(**kwargs),
        net_charge = net_charge,
        net_spin   = net_spin,
        **valency
        )
    
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
