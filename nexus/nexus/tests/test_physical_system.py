import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.PHYSICAL_SYSTEM)

from ..generic import generic_settings
generic_settings.raise_error = True

import numpy as np
from ..testing import value_eq, object_eq
from ..physical_system import Ion, PhysicalSystem#, Electrons
from ..periodic_table import Elements

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


# def test_electrons():
#     ref_charge = -10
#     ref_multiplicity = 3
#     ref_n_up = 6
#     ref_n_down = 4

#     electrons = Electrons(count=10, spin=1)

#     assert(electrons.charge       == ref_charge)
#     assert(electrons.multiplicity == ref_multiplicity)
#     assert(electrons.n_up         == ref_n_up)
#     assert(electrons.n_down       == ref_n_down)

#     ref_charge = -35
#     ref_multiplicity = 2
#     ref_n_up = 18
#     ref_n_down = 17

#     electrons = Electrons(count=35, spin=0.5)

#     assert(electrons.charge       == ref_charge)
#     assert(electrons.multiplicity == ref_multiplicity)
#     assert(electrons.n_up         == ref_n_up)
#     assert(electrons.n_down       == ref_n_down)


def test_custom_ion():
    ion = Ion(
        element     = Elements.Iron,
        label       = "Fe1",
        charge      = 2,
        spin        = 1,
        mass_number = 54,
        Zeff        = 16,
    )

    assert(ion.element       is Elements.Iron)
    assert(ion.label         == "Fe1")
    assert(ion.charge        == 2)
    assert(ion.spin          == 1)
    assert(ion.mass_number   == 54)
    assert(ion.Zeff          == 16)
    assert(ion.is_pseudo()   is True)
    assert(ion.is_ghost()    is False)
    assert(ion.name          == "Iron")
    assert(ion.symbol        == "Fe")
    assert(ion.atomic_weight == 55.845)
    assert(ion.mass          == 55.845)
    assert(ion.atomic_number == 26)
    assert(ion.protons       == 26)
    assert(ion.neutrons      == 28)


def test_minimal_ion():
    """Test to make sure the defaults are populated correctly."""
    ion = Ion("Fe")

    assert(ion.element       is Elements.Iron)
    assert(ion.label         == "Fe")
    assert(ion.charge        == 0)
    assert(ion.spin          == 0)
    assert(ion.mass_number   == 56)
    assert(ion.Zeff          is None)
    assert(ion.is_pseudo()   is False)
    assert(ion.is_ghost()    is False)
    assert(ion.name          == Elements.Iron.name)
    assert(ion.symbol        == Elements.Iron.symbol)
    assert(ion.atomic_weight == Elements.Iron.atomic_weight)
    assert(ion.mass          == Elements.Iron.atomic_weight)
    assert(ion.atomic_number == Elements.Iron.atomic_number)
    assert(ion.protons       == Elements.Iron.atomic_number)
    assert(ion.neutrons      == Elements.Iron.neutrons())


def test_ion_setters():
    """Test to make sure the setters for the ``Ion`` class work.

    This test also checks to make sure there are no side-effects from the
    setters. This means we check every property after each setter is used.
    """
    ion = Ion(
        element     = Elements.Iron,
        mass_number = 54,
    )
    # Mass number is set to 54 so we can use Chromium's isotopes without a warning

    assert(ion.element       is Elements.Iron)
    assert(ion.label         == "Fe")
    assert(ion.charge        == 0)
    assert(ion.spin          == 0)
    assert(ion.mass_number   == 54)
    assert(ion.Zeff          is None)
    assert(ion.is_pseudo()   is False)
    assert(ion.is_ghost()    is False)
    assert(ion.name          == Elements.Iron.name)
    assert(ion.symbol        == Elements.Iron.symbol)
    assert(ion.atomic_weight == Elements.Iron.atomic_weight)
    assert(ion.atomic_number == Elements.Iron.atomic_number)
    assert(ion.neutrons      == Elements.Iron.neutrons(mass_number=54))

    ref_element = Elements.Chromium
    ion.atomic_number = 24

    assert(ion.element       is ref_element)
    assert(ion.label         == "Fe") # We don't want to change custom labels
    assert(ion.charge        == 0)
    assert(ion.spin          == 0)
    assert(ion.mass_number   == 54) # The mass number also shouldn't get changed
    assert(ion.Zeff          is None)
    assert(ion.is_pseudo()   is False)
    assert(ion.is_ghost()    is False)
    assert(ion.name          == ref_element.name)
    assert(ion.symbol        == ref_element.symbol)
    assert(ion.atomic_weight == ref_element.atomic_weight)
    assert(ion.atomic_number == ref_element.atomic_number)
    assert(ion.neutrons      == ref_element.neutrons(mass_number=54))

    new_mass_number = 52
    ion.neutrons = 28

    assert(ion.element       is ref_element)
    assert(ion.label         == "Fe")
    assert(ion.charge        == 0)
    assert(ion.spin          == 0)
    assert(ion.mass_number   == new_mass_number)
    assert(ion.Zeff          is None)
    assert(ion.is_pseudo()   is False)
    assert(ion.is_ghost()    is False)
    assert(ion.name          == ref_element.name)
    assert(ion.symbol        == ref_element.symbol)
    assert(ion.atomic_weight == ref_element.atomic_weight)
    assert(ion.atomic_number == ref_element.atomic_number)
    assert(ion.neutrons      == ref_element.neutrons(mass_number=new_mass_number))



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
