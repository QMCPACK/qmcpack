#!/env/bin/python

import numpy as np
from testing import value_eq,object_eq,object_diff


reference_inputs     = dict()
reference_structures = dict()
generated_structures = dict()
crystal_structures   = dict()


def structure_same(s1,s2):
    keys = ('units','elem','pos','axes','kpoints','kweights','kaxes')
    o1 = s1.obj(keys)
    o2 = s2.obj(keys)
    return object_eq(o1,o2)
#end def structure_same


def structure_diff(s1,s2):
    keys = ('units','elem','pos','axes','kpoints','kweights','kaxes')
    o1 = s1.obj(keys)
    o2 = s2.obj(keys)
    return object_diff(o1,o2,full=True)
#end def structure_diff


def get_reference_inputs():
    from generic import obj
    if len(reference_inputs)==0:
        ref_in = obj()
        ref_in.diamond_prim = obj(
            units = 'A',
            axes  = [[1.785, 1.785, 0.   ],
                     [0.   , 1.785, 1.785],
                     [1.785, 0.   , 1.785]],
            elem  = ['C','C'],
            pos   = [[0.    , 0.    , 0.    ],
                     [0.8925, 0.8925, 0.8925]],
            )
        ref_in.diamond_conv = obj(
            units = 'A',
            axes  = [[3.57, 0   , 0.  ],
                     [0.  , 3.57, 0.  ],
                     [0   , 0.  , 3.57]],
            elem  = ['C','C','C','C','C','C','C','C'],
            pos   = [[0.0000, 0.0000, 0.0000],
                     [0.8925, 0.8925, 0.8925],
                     [0.0000, 1.7850, 1.7850],
                     [0.8925, 2.6775, 2.6775],
                     [1.7850, 0.0000, 1.7850],
                     [2.6775, 0.8925, 2.6775],
                     [1.7850, 1.7850, 0.0000],
                     [2.6775, 2.6775, 0.8925]],
            )
        ref_in.wurtzite_prim = obj(
            units = 'A',
            axes  = [[ 3.350, 0.00000000, 0.00],
                     [-1.675, 2.90118510, 0.00],
                     [ 0.000, 0.00000000, 5.22]],
            elem  = ['Zn','O','Zn','O'],
            pos   = [[0.00000000e+00, 0.00000000e+00, 3.26250000e+00],
                     [0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                     [1.67500000e+00, 9.67061701e-01, 6.52500000e-01],
                     [1.67500000e+00, 9.67061701e-01, 2.61000000e+00]],
            )
        ref_in.oxygen_prim = obj(
            units = 'A',
            axes  = [[ 2.70150000e+00, -1.71450000e+00,  0.00000000e+00],
                     [ 2.70150000e+00,  1.71450000e+00,  0.00000000e+00],
                     [-3.43801471e+00,  0.00000000e+00,  3.74799291e+00]],
            elem  = ['O','O'],
            pos   = [[ 0.,    0.,     0.575],
                     [ 0.,    0.,    -0.575]],
            )
        ref_in.CuO_prim = obj(
            units = 'A',
            axes  = [[ 2.34150000e+00, -1.71100000e+00, 0.00000000e+00],
                     [ 2.34150000e+00,  1.71100000e+00, 0.00000000e+00],
                     [-8.49894838e-01,  0.00000000e+00, 5.05708046e+00]],
            elem  = ['Cu','O','Cu','O'],
            pos   = [[ 1.17075   ,  0.8555  ,  0.        ],
                     [-0.21247371,  1.430396,  1.26427011],
                     [ 0.74580258,  2.5665  ,  2.52854023],
                     [ 1.70407887,  0.280604,  3.79281034]],
            )
        ref_in.Ca2CuO3_prim = obj(
            units = 'A',
            axes  = [[ 1.885,  1.625, -6.115],
                     [-1.885,  1.625,  6.115],
                     [ 1.885, -1.625,  6.115]],
            elem  = ['Cu','O','O','O','Ca','Ca'],
            pos   = [[0.000, 0.000,  0.000],
                     [1.885, 0.000,  0.000],
                     [0.000, 0.000,  1.960],
                     [0.000, 0.000, 10.270],
                     [0.000, 0.000,  4.290],
                     [0.000, 0.000,  7.940]],
            )
        ref_in.La2CuO4_prim = obj(
            units = 'A',
            axes  = [[ 1.9045,  1.9045, -6.5845],
                     [-1.9045,  1.9045,  6.5845],
                     [ 1.9045, -1.9045,  6.5845]],
            elem  = ['Cu','O','O','O','O','La','La'],
            pos   = [[ 0.      ,  0.      ,  0.      ],
                     [ 0.95225 ,  0.95225 , -3.29225 ],
                     [-0.95225 ,  0.95225 ,  3.29225 ],
                     [ 0.346619, -0.346619,  1.198379],
                     [-0.346619,  0.346619, -1.198379],
                     [ 0.689429, -0.689429,  2.383589],
                     [-0.689429,  0.689429, -2.383589]],
            )
        ref_in.graphene_prim = obj(
            units = 'A',
            axes  = [[ 2.462, 0.00000000,  0.00],
                     [-1.231, 2.13215454,  0.00],
                     [ 0.000, 0.00000000, 15.00]],
            elem  = ['C','C'],
            pos   = [[0.   , 0.        , 0. ],
                     [1.231, 0.71071818, 0. ]],
            )
        for k,v in ref_in.items():
            reference_inputs[k] = v
        #end for
    #end if
    return obj(reference_inputs)
#end def get_reference_inputs


def get_reference_structures():
    from generic import obj
    from structure import Structure
    if len(reference_structures)==0:
        ref_in = get_reference_inputs()
        ref = reference_structures
        for name,inputs in ref_in.items():
            ref[name] = Structure(**inputs)
        #end for
    #end if
    return obj(reference_structures)
#end def get_reference_structures


def get_generated_structures():
    from generic import obj
    from structure import generate_structure
    if len(generated_structures)==0:
        ref_in = get_reference_inputs()
        gen = generated_structures
        for name,inputs in ref_in.iteritems():
            gen[name] = generate_structure(**inputs)
        #end for
    #end if
    return obj(generated_structures)
#end def get_generated_structures


def get_crystal_structures():
    from generic import obj
    from structure import Crystal,generate_structure
    if len(crystal_structures)==0:
        crys = crystal_structures
        for (latt,cell),inputs in Crystal.known_crystals.iteritems():
            s = generate_structure(structure=latt,cell=cell)
            crys[latt+'_'+cell] = s
        #end for
    #end if
    return obj(crystal_structures)
#end def get_crystal_structures



def test_import():
    from structure import Structure,Crystal
    from structure import generate_structure
    from structure import read_structure
#end def test_import



def test_empty_init():
    from structure import Structure
    from structure import generate_structure
    s1 = Structure()
    s2 = generate_structure('empty')
    assert(object_eq(s1,s2))
#end def test_empty_init



def test_reference_inputs():
    from generic import obj
    ref_in = get_reference_inputs()
    assert(len(ref_in)>0)
#end def test_reference_inputs



def test_direct_init():
    ref = get_reference_structures()
    assert(len(ref)>0)
#end def test_direct_init



def test_generate_init():
    ref = get_reference_structures()
    gen = get_generated_structures()
    assert(len(gen)>0)

    for name,s in ref.items():
        assert(structure_same(s,gen[name]))
    #end for
#end def test_generate_init



def test_crystal_init():
    ref = get_reference_structures()
    crys = get_crystal_structures()
    assert(len(crys)>0)

    for name,s in ref.items():
        assert(structure_same(s,crys[name]))
    #end for
#end def test_crystal_init



def test_diagonal_tiling():
    ref = get_reference_structures()
    diag_tilings = [
        (1, 1, 1),
        (1, 2, 3),
        (1, 2, 4),
        (1, 3, 1),
        (1, 3, 2),
        (1, 3, 3),
        (1, 5, 3),
        (1, 5, 5),
        (1, 5, 6),
        (2, 1, 2),
        (2, 1, 3),
        (2, 2, 2),
        (2, 3, 1),
        (2, 3, 2),
        (2, 3, 3),
        (2, 4, 4),
        (2, 5, 1),
        (2, 5, 2),
        (2, 5, 3),
        (3, 1, 1),
        (3, 1, 2),
        (3, 1, 3),
        (3, 2, 1),
        (3, 2, 3),
        (3, 3, 6),
        (4, 6, 4),
        (5, 5, 1),
        (5, 5, 5),
        (6, 2, 6),
        (6, 3, 1),
        (6, 3, 6),
        (6, 4, 6),
        (6, 6, 4),
        ]
    for name,s in ref.items():
        for tvec in diag_tilings:
            st = s.tile(tvec)
            st.check_tiling()
        #end for
    #end for
#end test diagonal_tiling



def test_matrix_tiling():
    ref = get_reference_structures()
    matrix_tilings = [
        (-3,-2, 0,-2,-2,-3, 3, 2,-1),
        (-3,-1, 3, 0, 1, 1,-3, 3,-3),
        (-3, 0,-2,-3, 0, 3, 2,-3, 0),
        (-3, 1, 0,-3, 3,-3, 0,-2, 1),
        (-2,-1, 2, 3,-3,-3,-2, 2,-2),
        (-2, 0, 3,-2, 1,-1,-3, 0, 1),
        (-2, 1,-3,-3, 0, 3, 2, 0, 3),
        (-2, 2,-3,-1, 0, 1,-1,-2, 3),
        (-1, 3, 2,-3, 2,-2, 3,-3,-1),
        (-1, 3, 3,-3, 3, 0,-3,-1, 1),
        ( 0,-1,-2, 2, 3, 1, 3,-2, 0),
        ( 0, 1, 3,-2, 3, 1, 2, 1, 2),
        ( 0, 3,-1,-1,-1, 2, 2, 3, 2),
        ( 1,-3,-3,-3, 3, 0,-2, 2,-2),
        ( 1,-2, 2, 1, 1, 1, 3, 2,-1),
        ( 1, 0, 1,-3, 2, 1, 1, 2,-2),
        ( 1, 2, 1,-2, 1,-1,-2, 3, 2),
        ( 1, 2, 1, 2, 1,-1, 3, 2, 2),
        ( 1, 2, 2,-2, 1, 0, 1, 0,-1),
        ( 1, 2, 2, 1, 0, 1,-3, 0,-2),
        ( 1, 3,-1, 0, 1,-2, 1, 0, 2),
        ( 1, 3, 0,-1, 0,-3, 2, 0, 2),
        ( 1, 3, 0, 0,-2,-2,-1,-1, 1),
        ( 1, 3, 1,-3, 2, 0, 0, 0, 2),
        ( 2,-3, 0, 2, 2,-3,-1,-2, 0),
        ( 2,-2, 1, 3, 3,-2,-2, 1,-3),
        ( 2,-1, 2,-1, 0,-3, 2, 0, 1),
        ( 2,-1, 2, 0,-2, 0,-2, 0, 0),
        ( 2,-1, 3, 2, 1, 2, 1,-1,-3),
        ( 2, 1, 2, 0,-3,-2,-2,-3,-1),
        ( 2, 2, 1,-3,-2,-3, 2,-1,-2),
        ( 2, 3,-3, 1, 3, 1, 2,-2,-1),
        ( 3,-3, 0, 1,-2, 3, 2, 1, 1),
        ( 3,-2, 3, 3,-3,-3, 0, 0,-1),
        ( 3,-1, 2,-3, 2,-2, 3, 2,-1),
        ( 3, 0,-3, 0, 2,-1, 3,-2,-1),
        ( 3, 0, 2, 3,-2,-3,-3,-2, 1),
        ( 3, 1, 1,-2,-1,-3, 3, 1, 0),
        ( 3, 2,-3,-2, 2, 1, 1, 0,-1),
        ( 3, 2,-2,-1,-3,-1, 3, 3, 2),
        ( 3, 2, 3,-1, 2, 2, 0,-1, 0),
        ( 3, 3,-1, 0,-3, 2,-3,-1, 0),
        ( 3, 3, 1,-3,-2,-3, 2,-3, 3),
        ( 3, 3, 1, 1,-2, 2, 2,-2, 0),
        ]
    #for name in sorted(ref.keys()):
    for name in ['diamond_prim']:
        s = ref[name]
        npass = 0
        for tmat in matrix_tilings:
            tmat = np.array(tmat,dtype=int)
            tmat.shape = 3,3
            st = s.tile(tmat)
            st.check_tiling()
        #end for
    #end for
#end def test_matrix_tiling



def test_count_kshells():
    from test_physical_system import example_structure_h4
    s1 = example_structure_h4()
    kf = 1.465
    kcut = 5*kf
    nksh = s1.count_kshells(kcut)
    assert(nksh==13)
#end def test_count_kshells


def test_rinscribe():
    from structure import generate_structure
    s = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        tiling    = (4,4,1),
        )
    ri = s.rinscribe()
    assert(value_eq(float(ri),4.2643090882345795))
#end def test_rinscribe
