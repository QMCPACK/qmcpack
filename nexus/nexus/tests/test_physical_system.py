import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.PHYSICAL_SYSTEM)

from ..generic import generic_settings
generic_settings.raise_error = True

import numpy as np
from ..testing import value_eq, object_eq
from ..physical_system import IonSpecies, PhysicalSystem, Electrons, Positrons
from ..periodic_table import Elements
from ..structure import Structure

from .test_structure import structure_same


def system_same(s1,s2,pseudized=True,tiled=False):
    same = True
    keys = ('net_charge','net_spin','pseudized','particles')
    o1 = s1.obj(keys)
    o2 = s2.obj(keys)
    qsame = object_eq(o1,o2)
    vsame = True
    if pseudized:
        vsame = s1.valency==s2.valency
    #end if
    ssame = structure_same(s1.structure,s2.structure)
    fsame = True
    if tiled:
        fsame = system_same(s1.folded_system,s2.folded_system)
    #end if
    same = qsame and vsame and ssame and fsame
    return same
#end def system_same


def test_electrons():
    ref_charge       = -16
    ref_multiplicity = 3
    ref_n_up         = 9
    ref_n_down       = 7

    electrons = Electrons(count=16, spin=1, spin_orbit=False)

    assert(electrons.total_charge    == ref_charge)
    assert(electrons.multiplicity    == ref_multiplicity)
    assert(electrons.n_up            == ref_n_up)
    assert(electrons.n_down          == ref_n_down)
    assert(not electrons.is_fractional())

    ref_charge       = -15
    ref_multiplicity = 3
    ref_n_up         = 8.5
    ref_n_down       = 6.5

    electrons = Electrons(count=15, spin=1, spin_orbit=False)

    assert(electrons.total_charge == ref_charge)
    assert(not electrons.is_fractional())
    assert(electrons.multiplicity == ref_multiplicity)
    assert(electrons.n_up         == ref_n_up)
    assert(electrons.n_down       == ref_n_down)

    ref_charge       = -16
    ref_multiplicity = 2
    ref_n_up         = 8.5
    ref_n_down       = 7.5

    electrons = Electrons(count=16, spin=0.5, spin_orbit=False)

    assert(electrons.total_charge == ref_charge)
    assert(not electrons.is_fractional())
    assert(electrons.multiplicity == ref_multiplicity)
    assert(electrons.n_up         == ref_n_up)
    assert(electrons.n_down       == ref_n_down)

    ref_charge       = -15
    ref_multiplicity = 2
    ref_n_up         = 8
    ref_n_down       = 7

    electrons = Electrons(count=15, spin=0.5, spin_orbit=False)

    assert(electrons.total_charge == ref_charge)
    assert(not electrons.is_fractional())
    assert(electrons.multiplicity == ref_multiplicity)
    assert(electrons.n_up         == ref_n_up)
    assert(electrons.n_down       == ref_n_down)

    ref_charge       = -15
    ref_multiplicity = 2
    ref_n_up         = 7
    ref_n_down       = 8

    electrons = Electrons(count=15, spin=-0.5, spin_orbit=False)

    assert(electrons.total_charge == ref_charge)
    assert(not electrons.is_fractional())
    assert(electrons.multiplicity == ref_multiplicity)
    assert(electrons.n_up         == ref_n_up)
    assert(electrons.n_down       == ref_n_down)


def test_electrons_eq():
    ref_charge       = -16
    ref_multiplicity = 3
    ref_n_up         = 9
    ref_n_down       = 7

    electrons1 = Electrons(count=16, spin=1, spin_orbit=False)

    assert(electrons1.total_charge    == ref_charge)
    assert(electrons1.multiplicity    == ref_multiplicity)
    assert(electrons1.n_up            == ref_n_up)
    assert(electrons1.n_down          == ref_n_down)
    assert(not electrons1.is_fractional())

    electrons2 = Electrons(count=16, spin=1, spin_orbit=False)

    assert(electrons2.total_charge    == ref_charge)
    assert(electrons2.multiplicity    == ref_multiplicity)
    assert(electrons2.n_up            == ref_n_up)
    assert(electrons2.n_down          == ref_n_down)
    assert(not electrons2.is_fractional())

    assert(electrons1 == electrons2)


def test_positrons():
    """More minimal test since the code is the same as for ``Electrons``."""

    ref_charge       = 16
    ref_multiplicity = 3
    ref_n_up         = 9
    ref_n_down       = 7

    positrons1 = Positrons(count=16, spin=1, spin_orbit=False)

    assert(positrons1.total_charge    == ref_charge)
    assert(positrons1.multiplicity    == ref_multiplicity)
    assert(positrons1.n_up            == ref_n_up)
    assert(positrons1.n_down          == ref_n_down)
    assert(not positrons1.is_fractional())

    positrons2 = Positrons(count=16, spin=1, spin_orbit=False)

    assert(positrons2.total_charge    == ref_charge)
    assert(positrons2.multiplicity    == ref_multiplicity)
    assert(positrons2.n_up            == ref_n_up)
    assert(positrons2.n_down          == ref_n_down)
    assert(not positrons2.is_fractional())

    assert(positrons1 == positrons2)


def test_electron_positron_neq():

    ref_positron_charge       = 16
    ref_positron_multiplicity = 3
    ref_positron_n_up         = 9
    ref_positron_n_down       = 7

    positrons = Positrons(count=16, spin=1, spin_orbit=False)

    assert(positrons.total_charge    == ref_positron_charge)
    assert(positrons.multiplicity    == ref_positron_multiplicity)
    assert(positrons.n_up            == ref_positron_n_up)
    assert(positrons.n_down          == ref_positron_n_down)
    assert(not positrons.is_fractional())

    ref_electron_charge       = -16
    ref_electron_multiplicity = 3
    ref_electron_n_up         = 9
    ref_electron_n_down       = 7

    electrons = Electrons(count=16, spin=1, spin_orbit=False)

    assert(electrons.total_charge    == ref_electron_charge)
    assert(electrons.multiplicity    == ref_electron_multiplicity)
    assert(electrons.n_up            == ref_electron_n_up)
    assert(electrons.n_down          == ref_electron_n_down)
    assert(not electrons.is_fractional())

    assert(electrons != positrons)


def test_custom_ion_species():
    ion = IonSpecies(
        element     = Elements.Iron,
        count       = 12,
        label       = "Fe1",
        unit_charge = 2,
        unit_spin   = 1,
        Zeff        = 16,
        )

    assert(ion.element      is Elements.Iron)
    assert(ion.count        == 12)
    assert(ion.label        == "Fe1")
    assert(ion.unit_charge  == 2)
    assert(ion.unit_spin    == 1)
    assert(ion.Zeff         == 16)
    assert(ion.is_pseudo()  is True)
    assert(ion.is_ghost()   is False)
    assert(ion.symbol       == "Fe")
    assert(ion.total_mass   == 670.14)
    assert(ion.total_spin   == 12)
    assert(ion.total_charge == 24)


def test_minimal_ion_species():
    """Test to make sure the defaults are populated correctly."""
    ion = IonSpecies(element="Fe", count=12)

    assert(ion.element      is Elements.Iron)
    assert(ion.label        == "Fe")
    assert(ion.count        == 12)
    assert(ion.unit_charge  == 0)
    assert(ion.unit_spin    == 0)
    assert(ion.Zeff         is None)
    assert(ion.is_pseudo()  is False)
    assert(ion.is_ghost()   is False)
    assert(ion.symbol       == Elements.Iron.symbol)
    assert(ion.total_mass   == 670.14)
    assert(ion.total_spin   == 0)
    assert(ion.total_charge == 0)


def test_ion_species_eq():
    ref_element      = Elements.Iron
    ref_label        = "Fe1"
    ref_count        = 12
    ref_unit_charge  = 0
    ref_unit_spin    = 0
    ref_Zeff         = None
    ref_is_pseudo    = False
    ref_is_ghost     = False
    ref_symbol       = Elements.Iron.symbol
    ref_total_mass   = 670.14
    ref_total_charge = 0
    ref_total_spin   = 0

    ion1 = IonSpecies(element="Fe", label="Fe1", count=12)

    assert(ion1.element      is ref_element)
    assert(ion1.label        == ref_label)
    assert(ion1.count        == ref_count)
    assert(ion1.unit_charge  == ref_unit_charge)
    assert(ion1.unit_spin    == ref_unit_spin)
    assert(ion1.Zeff         is ref_Zeff)
    assert(ion1.is_pseudo()  is ref_is_pseudo)
    assert(ion1.is_ghost()   is ref_is_ghost)
    assert(ion1.symbol       == ref_symbol)
    assert(ion1.total_mass   == ref_total_mass)
    assert(ion1.total_charge == ref_total_charge)
    assert(ion1.total_spin   == ref_total_spin)

    ion2 = IonSpecies(element=Elements.Iron, label="Fe1", count=12)

    assert(ion2.element      is ref_element)
    assert(ion2.label        == ref_label)
    assert(ion2.count        == ref_count)
    assert(ion2.unit_charge  == ref_unit_charge)
    assert(ion2.unit_spin    == ref_unit_spin)
    assert(ion2.Zeff         is ref_Zeff)
    assert(ion2.is_pseudo()  is ref_is_pseudo)
    assert(ion2.is_ghost()   is ref_is_ghost)
    assert(ion2.symbol       == ref_symbol)
    assert(ion2.total_mass   == ref_total_mass)
    assert(ion2.total_charge == ref_total_charge)
    assert(ion2.total_spin   == ref_total_spin)

    assert(ion1 == ion2)


def test_ion_species_repr():
    ref_repr = "IonSpecies(element=Fe, count=12, label=Fe, unit_charge=0, unit_spin=0, Zeff=None)"
    ion = IonSpecies(element="Fe", count=12)
    assert(repr(ion) == ref_repr)


def test_ion_species_hash():
    ion1 = IonSpecies(
        element     = Elements.Carbon,
        count       = 1,
        label       = "C2",
        unit_charge = 0,
        unit_spin   = 0,
        Zeff        = 6,
        )
    ion2 = IonSpecies(
        element     = "C", # Use a slightly different element specifier to make sure it resolves.
        count       = 1,
        label       = "C2",
        unit_charge = 0,
        unit_spin   = 0,
        Zeff        = 6,
        )
    assert(ion1 == ion2) # Make sure they're actually equal.
    assert(hash(ion1) == hash(ion2))


def test_ion_species_from_structure():
    ref_ions = {
        "C1": IonSpecies(element=Elements.Carbon,   count=1, label="C1", unit_charge=0, unit_spin=0, Zeff=6),
        "C2": IonSpecies(element=Elements.Carbon,   count=1, label="C2", unit_charge=0, unit_spin=0, Zeff=6),
        "H" : IonSpecies(element=Elements.Hydrogen, count=5, label="H",  unit_charge=0, unit_spin=0, Zeff=1),
        "N" : IonSpecies(element=Elements.Nitrogen, count=1, label="N",  unit_charge=0, unit_spin=0, Zeff=7),
        "O" : IonSpecies(element=Elements.Oxygen,   count=2, label="O",  unit_charge=0, unit_spin=0, Zeff=8),
        }

    structure = Structure(
        elem = ["N", "C1", "C2", "O", "O", "H", "H", "H", "H", "H"],
        pos = np.zeros((10,3)),
        )
    ions = IonSpecies.from_structure(
        structure=structure,
        )

    assert(ions == ref_ions)

    ref_ions = {
        "C1": IonSpecies(element=Elements.Carbon,   count=1, label="C1", unit_charge=2, unit_spin=0,   Zeff=4),
        "C2": IonSpecies(element=Elements.Carbon,   count=1, label="C2", unit_charge=4, unit_spin=0.5, Zeff=6),
        "H" : IonSpecies(element=Elements.Hydrogen, count=5, label="H",  unit_charge=1, unit_spin=0.5, Zeff=1),
        "N" : IonSpecies(element=Elements.Nitrogen, count=1, label="N",  unit_charge=3, unit_spin=1,   Zeff=5),
        "O" : IonSpecies(element=Elements.Oxygen,   count=2, label="O",  unit_charge=2, unit_spin=0,   Zeff=6),
        }

    structure = Structure(
        elem = ["N", "C1", "C2", "O", "O", "H", "H", "H", "H", "H"],
        pos = np.zeros((10,3)),
        )
    ions = IonSpecies.from_structure(
        structure   = structure,
        elem_charge = dict(N=3, C1=2, C2=4,   O=2, H=1),
        elem_spin   = dict(N=1, C1=0, C2=0.5, O=0, H=0.5),
        elem_Zeff   = dict(N=5, C1=4, C2=6,   O=6, H=1),
        )

    assert(ions == ref_ions)


def test_physical_system_initialization(tmp_path):
    from ..developer import obj
    from ..structure import generate_structure
    from ..physical_system import generate_physical_system

    d2 = generate_structure(
        structure = 'diamond',
        cell      = 'prim',
        )
    d2_path = tmp_path / 'diamond2.xsf'
    d2.write(d2_path)

    d8 = generate_structure(
        structure = 'diamond',
        cell      = 'conv',
        )
    d8_path = tmp_path / 'diamond8.xsf'
    d8.write(d8_path)


    d8_tile = d2.tile([[ 1, -1,  1],
                       [ 1,  1, -1],
                       [-1,  1,  1]])

    d8_tile_pos_ref = np.array([
        [0.  , 0.  , 0.  ],
        [0.25, 0.25, 0.25],
        [0.5 , 0.5 , 0.  ],
        [0.75, 0.75, 0.25],
        [0.  , 0.5 , 0.5 ],
        [0.25, 0.75, 0.75],
        [0.5 , 0.  , 0.5 ],
        [0.75, 0.25, 0.75]])

    assert(value_eq(d8_tile.pos_unit(),d8_tile_pos_ref,atol=1e-8))


    direct_notile = generate_physical_system(
        units = 'A',
        axes  = [[3.57, 0.00, 0.00],
                 [0.00, 3.57, 0.00],
                 [0.00, 0.00, 3.57]],
        elem  = 8*['C'],
        posu  = [[0.00, 0.00, 0.00],
                 [0.25, 0.25, 0.25],
                 [0.00, 0.50, 0.50],
                 [0.25, 0.75, 0.75],
                 [0.50, 0.00, 0.50],
                 [0.75, 0.25, 0.75],
                 [0.50, 0.50, 0.00],
                 [0.75, 0.75, 0.25]],
        C     = 4,
        )

    direct_tile = generate_physical_system(
        units  = 'A',
        axes   = [[1.785, 1.785, 0.   ],
                  [0.   , 1.785, 1.785],
                  [1.785, 0.   , 1.785]],
        elem   = 2*['C'],
        posu   = [[0.00, 0.00, 0.00],
                  [0.25, 0.25, 0.25]],
        tiling = [[ 1, -1,  1],
                  [ 1,  1, -1],
                  [-1,  1,  1]],
        C      = 4,
        )

    struct_notile = generate_physical_system(
        structure = d8,
        C         = 4,
        )

    struct_tile = generate_physical_system(
        structure = d2,
        tiling    = [[ 1, -1,  1],
                     [ 1,  1, -1],
                     [-1,  1,  1]],
        C         = 4,
        )

    read_notile = generate_physical_system(
        structure = d8_path,
        C         = 4,
        )

    read_tile = generate_physical_system(
        structure = d2_path,
        tiling    = [[ 1, -1,  1],
                     [ 1,  1, -1],
                     [-1,  1,  1]],
        C         = 4,
        )

    gen_notile = generate_physical_system(
        lattice   = 'cubic',        # cubic tetragonal orthorhombic rhombohedral
                                    # hexagonal triclinic monoclinic
        cell      = 'conventional', # primitive, conventional
        centering = 'F',            # P A B C I R F
        constants = 3.57,           # a,b,c,alpha,beta,gamma
        units     = 'A',            # A or B
        atoms     = 'C',            # species in primitive cell
        basis     = [[0,0,0],       # basis vectors (optional)
                     [.25,.25,.25]],
        C         = 4,
        )

    gen_tile = generate_physical_system(
        lattice   = 'cubic',        # cubic tetragonal orthorhombic rhombohedral
                                    # hexagonal triclinic monoclinic
        cell      = 'primitive',    # primitive, conventional
        centering = 'F',            # P A B C I R F
        constants = 3.57,           # a,b,c,alpha,beta,gamma
        units     = 'A',            # A or B
        atoms     = 'C',            # species in primitive cell
        basis     = [[0,0,0],       # basis vectors (optional)
                     [.25,.25,.25]],
        tiling    = [[ 1, -1,  1],
                     [ 1,  1, -1],
                     [-1,  1,  1]],
        C         = 4,
        )

    lookup_notile = generate_physical_system(
        structure = 'diamond',
        cell      = 'conv',
        C         = 4,
        )

    lookup_tile = generate_physical_system(
        structure = 'diamond',
        cell      = 'prim',
        tiling    = [[ 1, -1,  1],
                     [ 1,  1, -1],
                     [-1,  1,  1]],
        C         = 4,
        )

    pref = obj(
        C = obj(
            charge          = 4,
            core_electrons  = 2,
            count           = 8,
            mass            = 21894.7135906,
            name            = 'C',
            neutrons        = 6,
            protons         = 6,
            spin            = 0,
            ),
        down_electron = obj(
            charge          = -1,
            count           = 16,
            mass            = 1.0,
            name            = 'down_electron',
            spin            = -1,
            ),
        up_electron = obj(
            charge          = -1,
            count           = 16,
            mass            = 1.0,
            name            = 'up_electron',
            spin            = 1,
            ),
        )

    # check direct system w/o tiling
    ref = direct_notile
    sref = ref.structure
    assert(ref.net_charge==0)
    assert(ref.net_spin==0)
    assert(ref.pseudized)
    assert(object_eq(ref.valency,obj(C=4)))
    assert(object_eq(ref.particles.to_obj(),pref))
    assert(structure_same(sref,d8))
    assert(value_eq(sref.axes,3.57*np.eye(3)))
    assert(tuple(sref.bconds)==tuple('ppp'))
    assert(list(sref.elem)==8*['C'])
    assert(value_eq(tuple(sref.pos[-1]),(2.6775,2.6775,0.8925)))
    assert(sref.units=='A')
    assert(object_eq(ref.particles.get_ions().to_obj(),obj(C=pref.C)))
    assert(object_eq(ref.particles.get_electrons().to_obj(),obj(down_electron=pref.down_electron,up_electron=pref.up_electron)))

    # check direct system w/ tiling
    ref = direct_tile
    sref = ref.structure
    assert(ref.net_charge==0)
    assert(ref.net_spin==0)
    assert(ref.pseudized)
    assert(object_eq(ref.valency,obj(C=4)))
    assert(object_eq(ref.particles.to_obj(),pref))
    assert(structure_same(sref,d8_tile))
    assert(value_eq(sref.axes,3.57*np.eye(3)))
    assert(tuple(sref.bconds)==tuple('ppp'))
    assert(list(sref.elem)==8*['C'])
    assert(value_eq(tuple(sref.pos[-1]),(2.6775,0.8925,2.6775)))
    assert(sref.units=='A')
    ref = direct_tile.folded_system
    sref = ref.structure
    pref.C.count = 2
    pref.down_electron.count = 4
    pref.up_electron.count = 4
    assert(ref.net_charge==0)
    assert(ref.net_spin==0)
    assert(ref.pseudized)
    assert(object_eq(ref.valency,obj(C=4)))
    assert(object_eq(ref.particles.to_obj(),pref))
    assert(structure_same(sref,d2))
    assert(value_eq(sref.axes,1.785*np.array([[1.,1,0],[0,1,1],[1,0,1]])))
    assert(tuple(sref.bconds)==tuple('ppp'))
    assert(list(sref.elem)==2*['C'])
    assert(value_eq(tuple(sref.pos[-1]),(0.8925,0.8925,0.8925)))
    assert(sref.units=='A')


    ref_notile = direct_notile
    ref_tile   = direct_tile

    assert(system_same(struct_notile,ref_notile))
    assert(system_same(read_notile  ,ref_notile))
    assert(system_same(gen_notile   ,ref_notile))
    assert(system_same(lookup_notile,ref_notile))

    assert(system_same(struct_tile,ref_tile,tiled=True))
    assert(system_same(read_tile  ,ref_tile,tiled=True))
    assert(system_same(gen_tile   ,ref_tile,tiled=True))
    assert(system_same(lookup_tile,ref_tile,tiled=True))

    systems_notile = [
        direct_notile,
        struct_notile,
        read_notile  ,
        gen_notile   ,
        lookup_notile,
        ]
    systems_tile = [
        direct_tile,
        struct_tile,
        read_tile  ,
        gen_tile   ,
        lookup_tile,
        ]
    systems = systems_notile+systems_tile
    for sys in systems:
        assert(sys.is_valid())
    #end for

    # test has_folded
    for sys in systems_notile:
        assert(not sys.has_folded())
    #end for
    for sys in systems_tile:
        assert(sys.has_folded())
    #end for

    # test copy
    for sys in systems:
        c = sys.copy()
        assert(id(c)!=id(sys))
        assert(c.is_valid())
        assert(system_same(c,sys,tiled=sys.has_folded()))
    #end for

    # test load
    for i,sys in enumerate(systems):
        path = tmp_path / 'system_{}'.format(i)
        sys.save(path)
        sys2 = PhysicalSystem()
        sys2.load(path)
        assert(sys2.is_valid())
        assert(system_same(sys2,sys,tiled=sys.has_folded()))
    #end for

    # test particle counts
    p = direct_notile.particles
    assert(p.count_ions()==8)
    assert(p.count_ions(species=True)==(8,1))
    assert(p.count_electrons()==32)
    assert(p.electron_counts()==[16,16])
#end def test_physical_system_initialization



def test_change_units():
    from ..physical_system import generate_physical_system

    sys = generate_physical_system(
        units = 'A',
        axes  = [[3.57, 0.00, 0.00],
                 [0.00, 3.57, 0.00],
                 [0.00, 0.00, 3.57]],
        elem  = 8*['C'],
        posu  = [[0.00, 0.00, 0.00],
                 [0.25, 0.25, 0.25],
                 [0.00, 0.50, 0.50],
                 [0.25, 0.75, 0.75],
                 [0.50, 0.00, 0.50],
                 [0.75, 0.25, 0.75],
                 [0.50, 0.50, 0.00],
                 [0.75, 0.75, 0.25]],
        C     = 4,
        )

    s = sys.structure

    assert(value_eq(s.pos[-1],np.array([2.6775,2.6775,0.8925])))
    sys.change_units('B')
    assert(value_eq(s.pos[-1],np.array([5.05974172,5.05974172,1.68658057])))
#end def test_change_units



def test_rename():
    from ..developer import obj
    from ..physical_system import generate_physical_system

    sys = generate_physical_system(
        units  = 'A',
        axes   = [[1.785, 1.785, 0.   ],
                  [0.   , 1.785, 1.785],
                  [1.785, 0.   , 1.785]],
        elem   = ['C1','C2'],
        posu   = [[0.00, 0.00, 0.00],
                  [0.25, 0.25, 0.25]],
        tiling = [[ 1, -1,  1],
                  [ 1,  1, -1],
                  [-1,  1,  1]],
        C1     = 4,
        C2     = 4,
        )

    ref = sys
    assert(object_eq(ref.valency,obj(C1=4,C2=4)))
    assert(list(ref.structure.elem)==4*['C1','C2'])
    assert(ref.particles.count_ions()==8)
    assert(ref.particles.count_ions(species=True)==(8,2))
    ref = sys.folded_system
    assert(object_eq(ref.valency,obj(C1=4,C2=4)))
    assert(list(ref.structure.elem)==['C1','C2'])
    assert(ref.particles.count_ions()==2)
    assert(ref.particles.count_ions(species=True)==(2,2))

    sys.rename(C1='C',C2='C')

    ref = sys
    assert(object_eq(ref.valency,obj(C=4)))
    assert(list(ref.structure.elem)==8*['C'])
    assert(ref.particles.count_ions()==8)
    assert(ref.particles.count_ions(species=True)==(8,1))
    ref = sys.folded_system
    assert(object_eq(ref.valency,obj(C=4)))
    assert(list(ref.structure.elem)==2*['C'])
    assert(ref.particles.count_ions()==2)
    assert(ref.particles.count_ions(species=True)==(2,1))

#end def test_rename



def test_tile():
    from ..physical_system import generate_physical_system

    d2_ref = generate_physical_system(
        units  = 'A',
        axes   = [[1.785, 1.785, 0.   ],
                  [0.   , 1.785, 1.785],
                  [1.785, 0.   , 1.785]],
        elem   = 2*['C'],
        posu   = [[0.00, 0.00, 0.00],
                  [0.25, 0.25, 0.25]],
        C      = 4,
        )

    d8_ref = generate_physical_system(
        units  = 'A',
        axes   = [[1.785, 1.785, 0.   ],
                  [0.   , 1.785, 1.785],
                  [1.785, 0.   , 1.785]],
        elem   = 2*['C'],
        posu   = [[0.00, 0.00, 0.00],
                  [0.25, 0.25, 0.25]],
        tiling = [[ 1, -1,  1],
                  [ 1,  1, -1],
                  [-1,  1,  1]],
        C      = 4,
        )

    d8 = d2_ref.tile([[ 1, -1,  1],
                      [ 1,  1, -1],
                      [-1,  1,  1]])

    assert(system_same(d8,d8_ref,tiled=True))
#end def test_tile



def test_kf_rpa():
    from .test_structure import example_structure_h4
    from ..physical_system import generate_physical_system
    s1 = example_structure_h4()
    ps = generate_physical_system(
        structure = s1,
        net_charge = 1,
        net_spin = 1,
        H = 1
        )
    kfs = ps.kf_rpa()
    assert np.isclose(kfs[0], 1.465, atol=1e-3)
    assert np.isclose(kfs[1], 1.465/2**(1./3), atol=1e-3)
#end def test_kf_rpa


@pytest.mark.xfail(
    raises=RuntimeError,
    reason="Can not split into up- and down-spin with spin_orbit=True",
)
def test_spin_orbit_fail_updown():
    electrons = Electrons(count=16, spin=1, spin_orbit=True)
    electrons.n_up_down()
#end def test_spin_orbit_fail_updown


@pytest.mark.xfail(
    raises=RuntimeError,
    reason="Can not split into up-spin with spin_orbit=True",
)
def test_spin_orbit_fail_up():
    electrons = Electrons(count=16, spin=1, spin_orbit=True)
    electrons.n_up
#end def test_spin_orbit_fail_up


@pytest.mark.xfail(
    raises=RuntimeError,
    reason="Can not split into down-spin with spin_orbit=True",
)
def test_spin_orbit_fail_down():
    electrons = Electrons(count=16, spin=1, spin_orbit=True)
    electrons.n_down
#end def test_spin_orbit_fail_down


@pytest.mark.xfail(
    raises=RuntimeError,
    reason="Multiplicity undefined with spin_orbit=True",
)
def test_spin_orbit_fail_multiplicity():
    electrons = Electrons(count=16, spin=1, spin_orbit=True)
    electrons.multiplicity
#end def test_spin_orbit_fail_multiplicity


@pytest.mark.xfail(
    raises=RuntimeError,
    reason="Multiplicity undefined with fractional count",
)
def test_fractional_count_fail_multiplicity():
    electrons = Electrons(count=16.5, spin=1, spin_orbit=False)
    electrons.multiplicity
#end def test_fractional_count_multiplicity
