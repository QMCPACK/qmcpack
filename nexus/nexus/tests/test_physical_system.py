import pytest
from copy import deepcopy

from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.PHYSICAL_SYSTEM)

import numpy as np
from ..structure import generate_structure
from ..testing import value_eq,object_eq
from ..physical_system import generate_physical_system
from ..pseudoset import PseudoSet

from .test_structure import structure_same


def sub_obj(s,keys):
    from ..developer import obj
    return obj({k:s[k] for k in keys})


def system_same(s1,s2,*,pseudized=True,tiled=False):
    same = True
    keys = ('net_charge','net_spin','pseudized')
    o1 = sub_obj(s1,keys)
    o2 = sub_obj(s2,keys)
    qsame = object_eq(o1,o2)
    vsame = True
    if pseudized:
        vsame = dict(**s1.valency)==dict(**s2.valency)
    #end if
    ssame = structure_same(s1.structure,s2.structure)
    fsame = True
    if tiled:
        fsame = system_same(s1.folded_system,s2.folded_system)
    #end if
    same = qsame and vsame and ssame and fsame
    return same
#end def system_same


def test_physical_system_initialization(tmp_path):
    from ..developer import obj
    from ..structure import generate_structure
    from ..physical_system import generate_physical_system
    from ..physical_system import PhysicalSystem

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

    # check direct system w/o tiling
    ref = direct_notile
    sref = ref.structure
    assert(ref.net_charge==0)
    assert(ref.net_spin==0)
    assert(ref.pseudized)
    assert(object_eq(ref.valency,obj(C=4)))
    assert(structure_same(sref,d8))
    assert(value_eq(sref.axes,3.57*np.eye(3)))
    assert(tuple(sref.bconds)==tuple('ppp'))
    assert(list(sref.elem)==8*['C'])
    assert(value_eq(tuple(sref.pos[-1]),(2.6775,2.6775,0.8925)))
    assert(sref.units=='A')

    # check direct system w/ tiling
    ref = direct_tile
    sref = ref.structure
    assert(ref.net_charge==0)
    assert(ref.net_spin==0)
    assert(ref.pseudized)
    assert(object_eq(ref.valency,obj(C=4)))
    assert(structure_same(sref,d8_tile))
    assert(value_eq(sref.axes,3.57*np.eye(3)))
    assert(tuple(sref.bconds)==tuple('ppp'))
    assert(list(sref.elem)==8*['C'])
    assert(value_eq(tuple(sref.pos[-1]),(2.6775,0.8925,2.6775)))
    assert(sref.units=='A')
    ref = direct_tile.folded_system
    sref = ref.structure
    assert(ref.net_charge==0)
    assert(ref.net_spin==0)
    assert(ref.pseudized)
    assert(object_eq(ref.valency,obj(C=4)))
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
        c = deepcopy(sys)
        assert(id(c)!=id(sys))
        assert(c.is_valid())
        assert(system_same(c,sys,tiled=sys.has_folded()))
    #end for

    # test load
    for i,sys in enumerate(systems):
        path = tmp_path / f'system_{i}'
        sys.save(path)
        sys2 = PhysicalSystem()
        sys2.load(path)
        assert(sys2.is_valid())
        assert(system_same(sys2,sys,tiled=sys.has_folded()))
    #end for

    # test particle counts
    assert(direct_notile.n_ions    ==  8)
    assert(direct_notile.n_species ==  1)
    assert(direct_notile.n_elec    == 32)
    assert(direct_notile.n_up      == 16)
    assert(direct_notile.n_down    == 16)
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
    assert(ref.n_ions==8)
    assert(ref.n_species==2)
    ref = sys.folded_system
    assert(object_eq(ref.valency,obj(C1=4,C2=4)))
    assert(list(ref.structure.elem)==['C1','C2'])
    assert(ref.n_ions==2)
    assert(ref.n_species==2)

    sys.rename(C1='C',C2='C')

    ref = sys
    assert(object_eq(ref.valency,obj(C=4)))
    assert(list(ref.structure.elem)==8*['C'])
    assert(ref.n_ions==8)
    assert(ref.n_species==1)
    ref = sys.folded_system
    assert(object_eq(ref.valency,obj(C=4)))
    assert(list(ref.structure.elem)==2*['C'])
    assert(ref.n_ions==2)
    assert(ref.n_species==1)

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


def test_update_pseudos(tmp_path):
    pseudo_names = (
        "C.BFD.xml",
        "H.BFD.xml",
        "O.BFD.xml",
        "C.BFD.upf",
        "H.BFD.upf",
        "O.BFD.upf",
        "Fe.BFD.upf",
        "C/POTCAR",
        "H/POTCAR",
        "O/POTCAR",
        )

    psp_dir = tmp_path / "pseudos"
    psp_dir.mkdir()
    assert psp_dir.exists(), "Failed to create pseudo directory!"

    pseudo_list = []
    for psp in pseudo_names:
        if "POTCAR" in psp:
            potcar_dir = psp_dir / psp.split("/")[0]
            potcar_dir.mkdir()
            assert potcar_dir.exists(), "Failed to create POTCAR directory!"

        pseudo = (psp_dir / psp).resolve()
        pseudo.touch()
        assert pseudo.exists(), "Failed to create pseudo file!"
        pseudo_list.append(pseudo)

    pseudoset = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        extensions = {"rmg": ".xml"},
        code_Zeff_map = {
            "qmcpack":  {"C": 4, "H": 1, "O": 6},
            "espresso": {"C": 4, "H": 1, "O": 6, "Fe": 16},
            "rmg":      {"C": 4, "H": 1, "O": 6},
            "vasp":     {"C": 4, "H": 1, "O": 6},
        }
    )

    structure = generate_structure(
        type       = "trimer",
        trimer     = ["O", "H", "H"],
        units      = "A",
        separation = [0.984, 0.984],
        angle      = 102.25,
    )

    mol = generate_physical_system(
        structure  = structure,
        net_charge = 0,
        pseudos    = pseudoset,
    )

    assert(mol.valency.__dict__ == {"H": 1, "O": 6})

    mol = generate_physical_system(
        structure  = structure,
        net_charge = 0,
        H = 2,
        O = 5
    )

    with pytest.raises(
        RuntimeError,
        match=(
            "The Z-valence for 'H' does not match the value in the provided pseudos!\n"
            "    Existing: 2\n"
            "    Provided: 1\n"
            "The Z-valence for 'O' does not match the value in the provided pseudos!\n"
            "    Existing: 5\n"
            "    Provided: 6\n"
        )
    ):
        mol.update_pseudos(pseudoset)
#end def test_update_pseudos

