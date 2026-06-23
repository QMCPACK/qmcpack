##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################

"""Representations of electrons, ions, and complete physical systems."""

from __future__ import annotations
import os
from pathlib import Path
from abc import ABC, abstractmethod
from copy import deepcopy
from typing import Self
import numpy as np
from .developer import DevBase, obj
from .periodic_table import Elements, ElementLike
from .structure import Structure, generate_structure, read_structure


class ElectronsPositronsBase(ABC):
    """Base class for ``Electrons`` and ``Positrons``.

    Note that this class does not make guarantees about having an
    integer amount of particles, but provides ``is_fractional`` to check
    for non-integer numbers of particles. If the class is created with
    integer-value floats for ``count`` and ``spin``, they will be
    converted into ints.

    Attributes
    ----------
    count : int or float, property
        The total number of particles
    spin : int or float, property
        The total spin of the system.

        An up-spin particle has a spin of +1/2, a down-spin particle has
        a spin of -1/2.
    n_up : int or float, read-only property
        The number of up-spin particles.
        Not defined for spin-orbit systems.
    n_down : int or float, read-only property
        The number of down-spin particles.
        Not defined for spin-orbit systems.
    total_charge : int or float, read-only property
        The total charge of the particles, equal to
        ``self.unit_charge * self.count``.
    multiplicity : int or float, read-only property
        The spin multiplicity of the particles, equal to :math:`2S+1`
        where :math:`S` is the spin.
    """

    @property
    @abstractmethod
    def unit_charge(self) -> int:
        ...

    def __init__(
        self,
        count     : int | float,
        spin      : int | float,
        spin_orbit: bool = False,
        ):
        self.count      = count
        self.spin       = spin
        self.spin_orbit = spin_orbit

    @property
    def count(self) -> int | float:
        return self._count

    @count.setter
    def count(self, new_count: int | float) -> None:
        if abs((int_count := int(new_count)) - new_count) < 1e-8:
            self._count = int_count
        else:
            self._count = new_count

    @property
    def spin(self) -> int | float:
        return self._spin

    @spin.setter
    def spin(self, new_spin: int | float) -> None:
        if abs((int_spin := int(new_spin)) - new_spin) < 1e-8:
            self._spin = int_spin
        else:
            self._spin = new_spin

    def n_up_down(self) -> tuple[int, int] | tuple[float, float]:
        """Return a tuple representing the number of up- and down-spin electrons.

        Examples
        --------
        >>> Electrons(count=16, spin=1).n_up_down()
        (9, 7) # (up, down)
        >>> Electrons(count=15, spin=1).n_up_down()
        (8.5, 6.5)
        >>> Electrons(count=16, spin=1/2).n_up_down()
        (8.5, 7.5)
        >>> Electrons(count=15, spin=1/2).n_up_down()
        (8, 7)
        >>> Electrons(count=15, spin=-1/2).n_up_down()
        (7, 8)
        """
        if self.spin_orbit:
            # Use self.__class__.__name__ to get name of subclass, not base class
            raise RuntimeError(
                f"{self.__class__.__name__} can not be split into up- and down-spin with a spin-orbit system!"
                )
        if isinstance(self.count, int):
            if isinstance(self.spin, int):
                if self.count % 2 == 0:
                    n_up   = (self.count // 2) + self.spin
                    n_down = (self.count // 2) - self.spin
                else:
                    n_up   = (self.count / 2) + self.spin
                    n_down = (self.count / 2) - self.spin
            else:
                if self.count % 2 == 0:
                    n_up   = (self.count / 2) + self.spin
                    n_down = (self.count / 2) - self.spin
                else:
                    if self.spin > 0:
                        n_down = self.count // 2
                        n_up   = self.count - n_down
                    else:
                        n_up   = self.count // 2
                        n_down = self.count - n_up
        else:
            n_up   = (self.count / 2) + self.spin
            n_down = (self.count / 2) - self.spin

        return n_up, n_down

    @property
    def n_up(self) -> int | float:
        if self.spin_orbit:
            raise RuntimeError(
                f"Up-spin {self.__class__.__name__} are not defined with a spin-orbit system!"
                )
        return self.n_up_down()[0]

    @property
    def n_down(self) -> int | float:
        if self.spin_orbit:
            raise RuntimeError(
                f"Down-spin {self.__class__.__name__} are not defined with a spin-orbit system!"
                )
        return self.n_up_down()[1]

    def is_fractional(self) -> bool:
        """Returns ``True`` if the count of particles is not an int."""
        return isinstance(self.count, float)

    @property
    def total_charge(self) -> int | float:
        return self.unit_charge * self.count

    @property
    def multiplicity(self) -> int | float:
        """Defined as :math:`2S+1` where :math:`S` is ``self.spin``.

        Undefined for spin-orbit systems and fractional counts.
        """
        if self.spin_orbit:
            raise RuntimeError("Multiplicity is undefined for spin-orbit systems!")
        elif self.is_fractional():
            raise RuntimeError("Multiplicity is undefined for fractional counts!")
        else:
            return (2 * abs(self.spin)) + 1

    def __eq__(self, other: Self) -> bool:
        return (
            self.unit_charge    == other.unit_charge
            and self.count      == other.count
            and self.spin       == other.spin
            and self.spin_orbit is other.spin_orbit
            )

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"count={self.count}, "
            f"spin={self.spin}, "
            f"spin_orbit={self.spin_orbit})"
            )
#end class ElectronsPositronsBase


class Electrons(ElectronsPositronsBase):
    """Class representing a collection of electrons."""
    @property
    def unit_charge(self) -> int:
        return -1
#end class Electrons


class Positrons(ElectronsPositronsBase):
    """Class representing a collection of positrons."""
    @property
    def unit_charge(self) -> int:
        return 1
#end class Positrons


class IonSpecies:
    """Class representing a collection of ions of the same type.

    Attributes
    ----------
    element : Elements
        The element for this ion collection.
    count : int or float
        The number of ions in this collection.
    label : str
        The label for the ion collection.
    unit_charge : int or float
        The charge associated with a single one of the ions.
    unit_spin : int or float
        The spin of a single one of the ions.
    Zeff : int or None
        The effective nuclear charge of one of the ion.
    symbol : str, read-only property
        The atomic symbol of the element.
    total_charge : float, read-only property
        Unit charge multiplied by count.
    total_spin : float, read-only property
        Unit spin multiplied by count.
    total_mass : float, read-only property
        Atomic weight multiplied by count.

    Parameters
    ----------
    element : ElementLike
        A member of the ``Elements`` enum, atomic symbol, or atomic
        number.
    count : int or float
        The number of ions in this collection.
    label : str, optional
        The label for the ion. If not given, defaults to
        ``element.symbol``.
    unit_charge : int or float, default=0
        The charge associated with a single one of the ions.
    unit_spin : int or float, default=0
        The spin of a single one of the ions.
    Zeff : int, optional
        The effective nuclear charge of the ion.
    """

    def __init__(
        self,
        element    : ElementLike,
        count      : int | float,
        label      : str | None  = None,
        unit_charge: int | float = 0,
        unit_spin  : int | float = 0,
        Zeff       : int | None  = None,
        ):
        self.element     = Elements(element)
        self.count       = count
        self.label       = label if label is not None else self.element.symbol
        self.unit_charge = unit_charge
        self.unit_spin   = unit_spin
        self.Zeff        = Zeff

    def is_pseudo(self) -> bool:
        """Check if this ion is pseudized."""
        return not (self.Zeff is None or self.Zeff == 0)

    def pseudize(self, Zeff: int):
        """Equivalent to setting ``self.Zeff = Zeff``."""
        self.Zeff = Zeff

    def is_ghost(self) -> bool:
        """Check if this collection of ions represents ghost atoms."""
        return self.element is Elements.Unknown

    @property
    def symbol(self) -> str:
        """Atomic symbol of the element."""
        return self.element.symbol

    @property
    def total_mass(self) -> float:
        """Total mass of all ions in the collection."""
        return self.element.atomic_weight * self.count

    @property
    def total_charge(self) -> int | float:
        """Total charge of all ions in the collection."""
        return self.unit_charge * self.count

    @property
    def total_spin(self) -> int | float:
        """Total spin of all ions in the collection."""
        return self.unit_spin * self.count

    def __eq__(self, other: Self) -> bool:
        return (
            self.element is other.element
            and self.count       == other.count
            and self.label       == other.label
            and self.unit_charge == other.unit_charge
            and self.unit_spin   == other.unit_spin
            and self.Zeff        == other.Zeff
            )

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"element={self.element}, "
            f"count={self.count}, "
            f"label={self.label}, "
            f"unit_charge={self.unit_charge}, "
            f"unit_spin={self.unit_spin}, "
            f"Zeff={self.Zeff})"
            )
#end class IonSpecies


class PhysicalSystem:
    """A system of electrons and ions with a structure.

    The ``PhysicalSystem`` is used to create inputs for all simulations.
    The difference between a ``PhysicalSystem`` and a ``Structure`` is
    that a ``Structure`` contains no information about charge, spin,
    pseudopotentials, elements, or electrons. The ``PhysicalSystem``
    class contains a ``Structure``, but it also includes the additional
    data required to define a real system that would be used in a
    calculation.

    Attributes
    ----------
    structure : Structure
        The structure of the ions in the system.
    ions : list of Ions
        The unique ions of the system. These are unique if and only if
        they share the same element, count, label, charge, spin, mass
        number, and effective nuclear charge. See ``Ions`` for more
        information.
    electrons : Electrons
        The up-spin and down-spin electrons in the system.
    particles : Particles
        A subclass of ``Particles`` that is defined by a simulation
        class.

    Parameters
    ----------
    structure : Structure
        The structure of the ions in the system.
    total_charge : int or float, default=0
        The total charge of the system.
    total_spin : int or float, default=0
        The total spin of the system.
    Zeffs : dict[ElementLike: int or float], default=dict()
        A dictionary mapping the effective nuclear charges of the system
        to the ions of the system.
    """

    def __init__(
        self,
        structure: Structure,
        ions     : dict[str: IonSpecies] | None = None,
        electrons: Electrons | None = None,
        positrons: Positrons | None = None,
        ):
        if structure is None:
            self.structure = Structure()
        else:
            self.structure = structure

        self.folded_system = None
        if self.structure.has_folded():
            if self.structure.is_tiled():
                vratio = structure.volume() / structure.folded_structure.volume()
                ncells = round(vratio)
                net_charge_fold = total_charge // ncells

                if abs(vratio - ncells) > 1e-4:
                    self.error('volume of system does not divide evenly into folded system')

                if total_charge % ncells != 0:
                    self.error('net charge of system does not divide evenly into folded system')

                if isinstance(total_spin, str):
                    net_spin_fold = total_spin
                elif total_spin % ncells != 0:
                    self.error('net_spin of system does not divide evenly into folded system')
                else:
                    net_spin_fold = total_spin // ncells

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
                total_charge = net_charge_fold,
                total_spin   = net_spin_fold,
                Zeffs        = Zeffs
                )
        #end if

        elems = list(self.structure.elem)
        ions = []
        for element_label in set(elems):
            Zeff = Zeffs.get(element_label)
            is_elem, element = Elements.is_element(element_label, return_element=True)
            if not is_elem:
                self.error(
                    f"Tried to initialize a physical system with unknown element {element}."
                    )
            ion = IonSpecies(
                element = element,
                count   = elems.count(element_label),
                label   = element_label,
                Zeff    = Zeff,
            )
            ions.append(ion)
        self.ions = ions

        self.valency_in = obj(Zeffs)

        self.check_folded_system()
    #end def __init__

    @property
    def ion_charge(self) -> int | float:
        ion_charge = 0
        for ion in self.ions:
            ion_charge += ion.total_charge

        return ion_charge

    @property
    def ion_spin(self) -> int | float:
        ion_spin = 0
        for ion in self.ions:
            ion_spin += ion.total_spin

        return ion_spin

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
        for name,count in particle_counts.items():
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
        for ion,valence_charge in valency.items():
            if ion in self.particles:
                ionp = self.particles[ion]
                if isinstance(ionp,IonSpecies):
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


    def rename(self,folded=True,**name_pairs):
        self.particles.rename(**name_pairs)
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
            self.valency_in = self.valency
        #end if
        if self.folded_system is not None and folded:
            self.folded_system.rename(folded=folded,**name_pairs)
        #end if
    #end def rename


    def copy(self):
        cp = DevBase.copy(self)
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
        system = self.copy()
        supersystem = PhysicalSystem(
            structure  = supercell,
            total_charge = net_charge,
            total_spin   = net_spin,
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
      nelecs = self.particles.electron_counts()
      volume = self.structure.volume()
      kvol1 = (2*np.pi)**3/volume  # k-space volume per particle
      kfs = [(3*nelec*kvol1/(4*np.pi))**(1./3) for nelec in nelecs]
      return np.array(kfs)
    #end def kf_rpa
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
        # particles = kwargs['particles']
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
                PhysicalSystem.class_error('pretile does not divide evenly into tiling\n  tiling provided: {0}\n  pretile provided: {1}'.format(tiling,pretile),'generate_physical_system')
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
            total_charge = net_charge,
            total_spin   = net_spin,
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
            total_charge = net_charge,
            total_spin   = net_spin,
            **valency
            )
        structure.set_folded(folded_structure)
        ps.folded_system = fps
    else:
        # No supercell tiling
        ps = PhysicalSystem(
            structure  = structure,
            total_charge = net_charge,
            total_spin   = net_spin,
            **valency
            )
    #end if

    ps.generation_info = generation_info

    return ps
#end def generate_physical_system



# test needed
def ghost_atoms(*particles):
    for particle in particles:
        ...
    #end for
#end def ghost_atoms
