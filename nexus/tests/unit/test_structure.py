#!/env/bin/python

import numpy as np
import versions
import testing
from testing import value_eq as value_eq_orig
from testing import object_eq as object_eq_orig
from testing import object_diff as object_diff_orig


struct_atol = 1e-10

def value_eq(*args,**kwargs):
    if 'atol' not in kwargs:
        kwargs['atol'] = struct_atol
    #end if
    return value_eq_orig(*args,**kwargs)
#end def value_eq

def object_eq(*args,**kwargs):
    if 'atol' not in kwargs:
        kwargs['atol'] = struct_atol
    #end if
    return object_eq_orig(*args,**kwargs)
#end def object_eq

def object_diff(*args,**kwargs):
    if 'atol' not in kwargs:
        kwargs['atol'] = struct_atol
    #end if
    return object_diff_orig(*args,**kwargs)
#end def object_diff


associated_files     = dict()
reference_inputs     = dict()
reference_structures = dict()
generated_structures = dict()
crystal_structures   = dict()


def get_files():
    return testing.collect_unit_test_file_paths('structure',associated_files)
#end def get_files


def structure_diff(s1,s2):
    keys = ('units','elem','pos','axes','kpoints','kweights','kaxes')
    o1 = s1.obj(keys)
    o2 = s2.obj(keys)
    return object_diff(o1,o2,full=True)
#end def structure_diff


def structure_same(s1,s2):
    import numpy as np
    keys = ('units','elem','axes','kpoints','kweights','kaxes','frozen','mag')
    o1 = s1.obj(keys)
    o2 = s2.obj(keys)
    osame = object_eq(o1,o2)
    psame = value_eq(s1.pos,s2.pos)
    if osame and not psame and len(s1.pos)==len(s2.pos):
        if s1.all_periodic() and s2.all_periodic():
            # wrap points around box edges
            u1 = s1.pos_unit()
            u1[np.abs(u1-1.0)<1e-10] = 0.0
            u2 = s2.pos_unit()
            u2[np.abs(u2-1.0)<1e-10] = 0.0
            psame = value_eq(u1,u2)
            if not psame:
                # ensure nearest neighbors have zero min image distance
                nt,dt = s1.neighbor_table(s1.pos,s2.pos,distances=True)
                dsame = np.abs(dt[:,0]).max()<1e-10
                # ensure species of nearest neighbors match
                esame = True
                for i,j in enumerate(nt[:,0]):
                    esame &= s1.elem[i]==s2.elem[j]
                #end for
                psame = esame and dsame
            #end if
        #end if
    #end if
    return osame and psame
#end def structure_same


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
        ref_in.diamond_64 = obj(
            structure = 'diamond',
            cell      = 'prim',
            tiling    = [[ 2, -2,  2],
                         [ 2,  2, -2],
                         [-2,  2,  2]],
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
            if 'cell' not in inputs:
                ref[name] = Structure(**inputs)
            #end if
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
        for name,inputs in ref_in.items():
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
        for (latt,cell),inputs in Crystal.known_crystals.items():
            s = generate_structure(structure=latt,cell=cell)
            crys[latt+'_'+cell] = s
        #end for
    #end if
    return obj(crystal_structures)
#end def get_crystal_structures


def example_structure_h4():
    # hydrogen at rs=1.31
    from structure import Structure
    natom = 4
    alat = 3.3521298178767225
    axes = alat*np.eye(3)
    elem = ['H']*natom
    pos = np.array([
      [0, 0, 0], [alat/2., 0, 0], [0, alat/2, 0], [0, 0, alat/2]
    ])
    s1 = Structure(axes=axes, elem=elem, pos=pos, units='B')
    return s1
#end def example_structure_h4



def test_files():
    filenames = [
        'La2CuO4_ICSD69312.cif',
        'coronene.xyz',
        ]
    files = get_files()
    assert(set(files.keys())==set(filenames))
#end def test_files



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



def test_change_units():
    import numpy as np
    ref = get_reference_structures()
    s = ref.diamond_conv.copy()
    assert(value_eq(s.pos[-1],np.array([2.6775,2.6775,0.8925])))
    s.change_units('B')
    assert(value_eq(s.pos[-1],np.array([5.05974172,5.05974172,1.68658057])))
#end def test_change_units



def test_rotate():
    import numpy as np
    ref = get_reference_structures()
    s0 = ref.CuO_prim.copy()
    s1 = ref.CuO_prim.copy()

    # Test the various parameter choices in the case that rp is given
    # Perform rotation taking x-axis to x-axis (original positions should be found)
    s1.rotate('x','x')
    assert(value_eq(s1.pos,s0.pos))
    assert(value_eq(s1.axes,s0.axes))
    assert(value_eq(s1.kaxes,s0.kaxes))
    # Perform active rotation taking x-coords to y-coords
    s1.rotate([1,0,0],[0,1,0])
    assert(value_eq(s1.pos,np.array([-s0.pos[:,1],s0.pos[:,0],s0.pos[:,2]]).T))
    assert(value_eq(s1.axes,np.array([-s0.axes[:,1],s0.axes[:,0],s0.axes[:,2]]).T))
    assert(value_eq(s1.kaxes,np.array([-s0.kaxes[:,1],s0.kaxes[:,0],s0.kaxes[:,2]]).T))
    # Perform passive rotation taking x-axis to y-axis (original positions should be found)
    s1.rotate('x','y',passive=True)
    assert(value_eq(s1.pos,s0.pos))
    assert(value_eq(s1.axes,s0.axes))
    assert(value_eq(s1.kaxes,s0.kaxes))
    # Perform active rotation about z-axis by an angle pi/2
    s1.rotate('z',np.pi/2.0)
    assert(value_eq(s1.pos,np.array([-s0.pos[:,1],s0.pos[:,0],s0.pos[:,2]]).T))
    assert(value_eq(s1.axes,np.array([-s0.axes[:,1],s0.axes[:,0],s0.axes[:,2]]).T))
    assert(value_eq(s1.kaxes,np.array([-s0.kaxes[:,1],s0.kaxes[:,0],s0.kaxes[:,2]]).T))
    # Perform active rotation taking y-coords to x-coords (original positions should be found)
    s1.rotate('y',[1,0,0])
    assert(value_eq(s1.pos[-1],s0.pos[-1]))
    assert(value_eq(s1.axes[-1],s0.axes[-1]))
    assert(value_eq(s1.kaxes[-1],s0.kaxes[-1]))
    # Perform active rotation taking a0-coords to a2-coords
    s1.rotate(s0.axes[0],s0.axes[2])
    assert(value_eq(s1.pos[-1],np.array([-2.15536928,3.46035669,0.86507139])))
    assert(value_eq(s1.axes[-1],np.array([-3.91292278,3.02549423,-1.35344154])))
    assert(value_eq(s1.kaxes[-1],np.array([-0.90768302,0.83458438,-0.15254555])))
    # Perform active rotation taking a2-coords to a0-coords (original positions should be found)
    s1.rotate('a2','a0')
    assert(value_eq(s1.pos,s0.pos))
    assert(value_eq(s1.axes,s0.axes))
    assert(value_eq(s1.kaxes,s0.kaxes))

    # Test the case where rp is not given
    # Perform active rotation taking a2-coords to a0-coords
    R = [[0.2570157723433977, 0.6326366344635742,-0.7305571719594085], 
         [0.4370696746690278, 0.5981289557203555, 0.6717230469572912], 
         [0.8619240060767753,-0.4919478031900122,-0.12277771249328594]]
    s1.rotate(R)
    assert(value_eq(s1.pos[-1],np.array([-2.15536928,3.46035669,0.86507139])))
    assert(value_eq(s1.axes[-1],np.array([-3.91292278,3.02549423,-1.35344154])))
    assert(value_eq(s1.kaxes[-1],np.array([-0.90768302,0.83458438,-0.15254555])))

    # A final test which places the structure back into its original form
    # Perform passive rotation taking a2-coords to a0-coords (original positions should be found)
    s1.rotate([-0.5871158698555267, -0.8034668669004766, -0.09867091342903483],1.7050154439645373,passive=True)
    assert(value_eq(s1.pos,s0.pos))
    assert(value_eq(s1.axes,s0.axes))
    assert(value_eq(s1.kaxes,s0.kaxes))

#end def test_change_units



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



def test_gen_molecule():
    """
    Create an H2O molecule.
    """
    from structure import generate_structure

    h2o = generate_structure(
        elem  = ['O','H','H'], 
        pos   = [[0.000000, 0.000000, 0.000000],
                 [0.000000,-0.757160, 0.586260],
                 [0.000000, 0.757160, 0.586260]],
        units = 'A', # Angstroms
        )

    # print important internal attributes
    assert(tuple(h2o.elem)==tuple('OHH'))
    assert(value_eq(h2o.pos[2,2],0.586260))
    assert(h2o.units=='A')
#end def demo_h2o_gen



def test_gen_diamond_direct():
    """
    Create a conventional cell of diamond directly.
    """
    import numpy as np
    from structure import generate_structure

    diamond = generate_structure(
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
        )

    assert(diamond.units=='A')
    assert(tuple(diamond.bconds)==tuple('ppp'))
    assert(value_eq(diamond.axes,3.57*np.eye(3)))
    assert(value_eq(list(diamond.center),3*[1.785]))
    assert(value_eq(diamond.kaxes,1.75999588*np.eye(3)))
    assert(list(diamond.elem)==8*['C'])
    assert(value_eq(tuple(diamond.pos[-1]),(2.6775,2.6775,0.8925)))
    assert(len(diamond.kpoints)==0)
    assert(len(diamond.kweights)==0)
    assert(diamond.frozen is None)
    assert(diamond.folded_structure is None)
#end def test_gen_diamond_direct



def test_gen_diamond_lattice():
    """
    Create a conventional cell of diamond from lattice information.
    """
    from structure import generate_structure

    diamond = generate_structure(
        lattice   = 'cubic',        # cubic tetragonal orthorhombic rhombohedral
                                    # hexagonal triclinic monoclinic
        cell      = 'conventional', # primitive, conventional
        centering = 'F',            # P A B C I R F
        constants = 3.57,           # a,b,c,alpha,beta,gamma
        units     = 'A',            # A or B
        atoms     = 'C',            # species in primitive cell
        basis     = [[0,0,0],       # basis vectors (optional)
                     [.25,.25,.25]],
        )

    ref = generate_structure(
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
        )

    assert(structure_same(diamond,ref))

#end def test_gen_diamond_lattice



def test_gen_graphene():
    """
    Create a graphene cell using stored information.
    """
    from structure import generate_structure

    graphene = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        )

    ref = generate_structure(
        units = 'A',
        axes  = [[ 2.46200000e+00,  0.00000000e+00,  0.00000000e+00],
                 [-1.23100000e+00,  2.13215454e+00,  0.00000000e+00],
                 [ 9.18485099e-16,  1.59086286e-15,  1.50000000e+01]],
        elem  = ['C','C'],
        pos   = [[0.   ,      0.        , 0.        ],
                 [1.231,      0.71071818, 0.        ]],
        )
        
    assert(structure_same(graphene,ref))
#end def test_gen_graphene



def test_read_write():
    """
    Write/read conventional diamond cell to/from XYZ, XSF, and POSCAR formats.
    """
    import os
    from structure import generate_structure, read_structure

    tpath = testing.setup_unit_test_output_directory('structure','test_read_write')

    d8 = generate_structure(
        structure = 'diamond',
        cell      = 'conv',
        )

    # Write an XYZ file
    xyz_file = os.path.join(tpath,'diamond8.xyz')
    d8.write(xyz_file)
    
    # Write an XSF file
    xsf_file = os.path.join(tpath,'diamond8.xsf')
    d8.write(xsf_file)

    # Write a POSCAR file
    poscar_file = os.path.join(tpath,'diamond8.POSCAR')
    d8.write(poscar_file)

    # Read an XYZ file
    d8_xyz = read_structure(xyz_file) # no cell info
    assert(value_eq(d8_xyz.elem,d8.elem))
    assert(value_eq(d8_xyz.pos,d8.pos))

    # Read an XSF file
    d8_xsf = read_structure(xsf_file)
    assert(structure_same(d8_xsf,d8))

    # Read a POSCAR file
    d8_poscar = read_structure(poscar_file)

    assert(structure_same(d8_poscar,d8))
#end def test_read_write



if versions.pycifrw_available and versions.cif2cell_available:
    def test_read_cif():
        """
        Read La2CuO4 structure from a CIF file.
        """
        from structure import read_structure,generate_structure

        files = get_files()

        # Read from CIF file
        s = read_structure(files['La2CuO4_ICSD69312.cif'])

        ref = generate_structure(
            units = 'A',
            axes  = [[ 2.665,   0.    , -6.5525],
                     [ 0.   ,   5.4126,  0.    ],
                     [ 2.665,   0.    ,  6.5525]],
            elem  = 'La La La La Cu Cu O O O O O O O O'.split(),
            pos   = [[ 2.665 ,      5.37038172, -1.8071795 ],
                     [ 2.665 ,      2.74851828,  4.7453205 ],
                     [ 2.665 ,      2.66408172, -4.7453205 ],
                     [ 2.665 ,      0.04221828,  1.8071795 ],
                     [ 0.    ,      0.        ,  0.        ],
                     [ 2.665 ,      2.7063    ,  0.        ],
                     [ 1.3325,      1.35315   , -0.128429  ],
                     [ 3.9975,      4.05945   ,  0.128429  ],
                     [ 3.9975,      1.35315   , -0.128429  ],
                     [ 1.3325,      4.05945   ,  0.128429  ],
                     [ 2.665 ,      0.23490684, -4.1398695 ],
                     [ 2.665 ,      2.47139316,  2.4126305 ],
                     [ 2.665 ,      2.94120684, -2.4126305 ],
                     [ 2.665 ,      5.17769316,  4.1398695 ]],
            )

        assert(structure_same(s,ref))

    #end def test_read_cif
#end if



def test_bounding_box():
    import numpy as np
    from structure import generate_structure,read_structure

    files = get_files()

    h2o = generate_structure(
        elem  = ['O','H','H'], 
        pos   = [[0.000000, 0.000000, 0.000000],
                 [0.000000,-0.757160, 0.586260],
                 [0.000000, 0.757160, 0.586260]],
        units = 'A', # Angstroms
        )

    # add a box by hand to the water molecule
    h2o_diy = h2o.copy()
    h2o_diy.set_axes(8.0*np.eye(3))
    assert(value_eq(h2o_diy.axes,8.*np.eye(3)))

    # automatically add a bounding box to the water molecule
    h2o_auto = h2o.copy()
    h2o_auto.bounding_box(box='cubic',minsize=8.0)
    assert(value_eq(h2o_auto.axes,8.*np.eye(3)))
    assert(value_eq(tuple(h2o_auto.pos[-1]),(4.,4.75716,4.29313)))


    s = read_structure(files['coronene.xyz'])

    # make a bounding box that is at least 5 A from the nearest atom
    s.bounding_box(mindist=5.0)

    ref_axes = np.array(
        [[16.4035,  0.    ,  0.    ],
         [ 0.    , 16.3786,  0.    ],
         [ 0.    ,  0.    , 10.    ]])

    assert(value_eq(s.axes,ref_axes))
    assert(value_eq(s.pos[:,2].min(),5.0))
    assert(value_eq(s.pos[:,2].max(),5.0))
    
#end def test_bounding_box



def test_opt_tiling():
    import numpy as np
    from structure import generate_structure

    dprim = generate_structure(
        structure = 'diamond',
        cell      = 'prim',
        )

    dopt = dprim.tile_opt(4)

    tmatrix_ref = np.array([[ 1, -1,  1],
                            [ 1,  1, -1],
                            [-1,  1,  1]])
    assert(value_eq(dopt.tmatrix,tmatrix_ref))
    assert(len(dopt.elem)==4*len(dprim.elem))
    assert(value_eq(dopt.axes,3.57*np.eye(3)))
#end def test_opt_tiling



if versions.seekpath_available:
    def test_primitive_search():
        """
        Find the primitive cell given a supercell.
        """

        from structure import generate_structure

        d2 = generate_structure(
            structure = 'diamond',
            cell      = 'prim',
            )

        tmatrix = [[ 2, -2,  2],
                   [ 2,  2, -2],
                   [-2,  2,  2]]

        d64 = d2.tile(tmatrix)

        # Remove all traces of the 2 atom cell, supercell is all that remains
        d64.remove_folded()
        del d2

        # Find the primitive cell from the supercell
        dprim = d64.primitive()

        tmatrix = d64.tilematrix(dprim)

        axes_ref = np.array([[0.   , 1.785, 1.785],
                             [1.785, 0.   , 1.785],
                             [1.785, 1.785, 0.   ]])
        tmatrix_ref = np.array([[-2,  2,  2],
                                [ 2, -2,  2],
                                [ 2,  2, -2]])

        assert(value_eq(dprim.axes,axes_ref))
        assert(value_eq(tmatrix,tmatrix_ref))

    #end def test_primitive_search
#end if



def test_unit_coords():
    """
    Get atomic positions in unit coordinates.
    """
    import numpy as np
    from structure import generate_structure

    s = generate_structure(
        structure = 'diamond',
        cell      = 'prim',
        tiling    = (2,2,2),
        )

    upos_ref = np.array([
        [ 0.000, 0.000, 0.000 ],
        [ 0.125, 0.125, 0.125 ],
        [ 0.500, 0.000, 0.000 ],
        [ 0.625, 0.125, 0.125 ],
        [ 0.000, 0.500, 0.000 ],
        [ 0.125, 0.625, 0.125 ],
        [ 0.500, 0.500, 0.000 ],
        [ 0.625, 0.625, 0.125 ],
        [ 0.000, 0.000, 0.500 ],
        [ 0.125, 0.125, 0.625 ],
        [ 0.500, 0.000, 0.500 ],
        [ 0.625, 0.125, 0.625 ],
        [ 0.000, 0.500, 0.500 ],
        [ 0.125, 0.625, 0.625 ],
        [ 0.500, 0.500, 0.500 ],
        [ 0.625, 0.625, 0.625 ]])
        
    upos = s.pos_unit()

    upos[np.abs(upos-1.0)<1e-10] = 0.0

    assert(value_eq(upos,upos_ref))

#end def test_unit_coords



def test_monkhorst_pack_kpoints():
    """
    Add k-points according to Monkhorst-Pack grid.
    """
    import numpy as np
    from generic import obj
    from structure import generate_structure

    g11 = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        )

    g44g_gen = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        tiling    = (4,4,1),
        kgrid     = (2,2,1),
        kshift    = (0,0,0),
        )

    g44s_gen = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        tiling    = (4,4,1),
        kgrid     = (2,2,1),
        kshift    = (0.5,0.5,0),
        )

    g44g_ukp_ref = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.0, 0.5, 0.0],
        [0.5, 0.5, 0.0]])

    g44s_ukp_ref = np.array([
        [0.25, 0.25, 0.00],
        [0.75, 0.25, 0.00],
        [0.25, 0.75, 0.00],
        [0.75, 0.75, 0.00]])

    g11g_ukp_ref = np.array([
        [ 0.000, 0.000, 0.000 ],
        [ 0.125, 0.000, 0.000 ],
        [ 0.000, 0.125, 0.000 ],
        [ 0.125, 0.125, 0.000 ],
        [ 0.250, 0.000, 0.000 ],
        [ 0.375, 0.000, 0.000 ],
        [ 0.250, 0.125, 0.000 ],
        [ 0.375, 0.125, 0.000 ],
        [ 0.500, 0.000, 0.000 ],
        [ 0.625, 0.000, 0.000 ],
        [ 0.500, 0.125, 0.000 ],
        [ 0.625, 0.125, 0.000 ],
        [ 0.750, 0.000, 0.000 ],
        [ 0.875, 0.000, 0.000 ],
        [ 0.750, 0.125, 0.000 ],
        [ 0.875, 0.125, 0.000 ],
        [ 0.000, 0.250, 0.000 ],
        [ 0.125, 0.250, 0.000 ],
        [ 0.000, 0.375, 0.000 ],
        [ 0.125, 0.375, 0.000 ],
        [ 0.250, 0.250, 0.000 ],
        [ 0.375, 0.250, 0.000 ],
        [ 0.250, 0.375, 0.000 ],
        [ 0.375, 0.375, 0.000 ],
        [ 0.500, 0.250, 0.000 ],
        [ 0.625, 0.250, 0.000 ],
        [ 0.500, 0.375, 0.000 ],
        [ 0.625, 0.375, 0.000 ],
        [ 0.750, 0.250, 0.000 ],
        [ 0.875, 0.250, 0.000 ],
        [ 0.750, 0.375, 0.000 ],
        [ 0.875, 0.375, 0.000 ],
        [ 0.000, 0.500, 0.000 ],
        [ 0.125, 0.500, 0.000 ],
        [ 0.000, 0.625, 0.000 ],
        [ 0.125, 0.625, 0.000 ],
        [ 0.250, 0.500, 0.000 ],
        [ 0.375, 0.500, 0.000 ],
        [ 0.250, 0.625, 0.000 ],
        [ 0.375, 0.625, 0.000 ],
        [ 0.500, 0.500, 0.000 ],
        [ 0.625, 0.500, 0.000 ],
        [ 0.500, 0.625, 0.000 ],
        [ 0.625, 0.625, 0.000 ],
        [ 0.750, 0.500, 0.000 ],
        [ 0.875, 0.500, 0.000 ],
        [ 0.750, 0.625, 0.000 ],
        [ 0.875, 0.625, 0.000 ],
        [ 0.000, 0.750, 0.000 ],
        [ 0.125, 0.750, 0.000 ],
        [ 0.000, 0.875, 0.000 ],
        [ 0.125, 0.875, 0.000 ],
        [ 0.250, 0.750, 0.000 ],
        [ 0.375, 0.750, 0.000 ],
        [ 0.250, 0.875, 0.000 ],
        [ 0.375, 0.875, 0.000 ],
        [ 0.500, 0.750, 0.000 ],
        [ 0.625, 0.750, 0.000 ],
        [ 0.500, 0.875, 0.000 ],
        [ 0.625, 0.875, 0.000 ],
        [ 0.750, 0.750, 0.000 ],
        [ 0.875, 0.750, 0.000 ],
        [ 0.750, 0.875, 0.000 ],
        [ 0.875, 0.875, 0.000 ]])

    g11s_ukp_ref = np.array([
        [ 0.0625, 0.0625, 0.0000 ],
        [ 0.1875, 0.0625, 0.0000 ],
        [ 0.0625, 0.1875, 0.0000 ],
        [ 0.1875, 0.1875, 0.0000 ],
        [ 0.3125, 0.0625, 0.0000 ],
        [ 0.4375, 0.0625, 0.0000 ],
        [ 0.3125, 0.1875, 0.0000 ],
        [ 0.4375, 0.1875, 0.0000 ],
        [ 0.5625, 0.0625, 0.0000 ],
        [ 0.6875, 0.0625, 0.0000 ],
        [ 0.5625, 0.1875, 0.0000 ],
        [ 0.6875, 0.1875, 0.0000 ],
        [ 0.8125, 0.0625, 0.0000 ],
        [ 0.9375, 0.0625, 0.0000 ],
        [ 0.8125, 0.1875, 0.0000 ],
        [ 0.9375, 0.1875, 0.0000 ],
        [ 0.0625, 0.3125, 0.0000 ],
        [ 0.1875, 0.3125, 0.0000 ],
        [ 0.0625, 0.4375, 0.0000 ],
        [ 0.1875, 0.4375, 0.0000 ],
        [ 0.3125, 0.3125, 0.0000 ],
        [ 0.4375, 0.3125, 0.0000 ],
        [ 0.3125, 0.4375, 0.0000 ],
        [ 0.4375, 0.4375, 0.0000 ],
        [ 0.5625, 0.3125, 0.0000 ],
        [ 0.6875, 0.3125, 0.0000 ],
        [ 0.5625, 0.4375, 0.0000 ],
        [ 0.6875, 0.4375, 0.0000 ],
        [ 0.8125, 0.3125, 0.0000 ],
        [ 0.9375, 0.3125, 0.0000 ],
        [ 0.8125, 0.4375, 0.0000 ],
        [ 0.9375, 0.4375, 0.0000 ],
        [ 0.0625, 0.5625, 0.0000 ],
        [ 0.1875, 0.5625, 0.0000 ],
        [ 0.0625, 0.6875, 0.0000 ],
        [ 0.1875, 0.6875, 0.0000 ],
        [ 0.3125, 0.5625, 0.0000 ],
        [ 0.4375, 0.5625, 0.0000 ],
        [ 0.3125, 0.6875, 0.0000 ],
        [ 0.4375, 0.6875, 0.0000 ],
        [ 0.5625, 0.5625, 0.0000 ],
        [ 0.6875, 0.5625, 0.0000 ],
        [ 0.5625, 0.6875, 0.0000 ],
        [ 0.6875, 0.6875, 0.0000 ],
        [ 0.8125, 0.5625, 0.0000 ],
        [ 0.9375, 0.5625, 0.0000 ],
        [ 0.8125, 0.6875, 0.0000 ],
        [ 0.9375, 0.6875, 0.0000 ],
        [ 0.0625, 0.8125, 0.0000 ],
        [ 0.1875, 0.8125, 0.0000 ],
        [ 0.0625, 0.9375, 0.0000 ],
        [ 0.1875, 0.9375, 0.0000 ],
        [ 0.3125, 0.8125, 0.0000 ],
        [ 0.4375, 0.8125, 0.0000 ],
        [ 0.3125, 0.9375, 0.0000 ],
        [ 0.4375, 0.9375, 0.0000 ],
        [ 0.5625, 0.8125, 0.0000 ],
        [ 0.6875, 0.8125, 0.0000 ],
        [ 0.5625, 0.9375, 0.0000 ],
        [ 0.6875, 0.9375, 0.0000 ],
        [ 0.8125, 0.8125, 0.0000 ],
        [ 0.9375, 0.8125, 0.0000 ],
        [ 0.8125, 0.9375, 0.0000 ],
        [ 0.9375, 0.9375, 0.0000 ]])

    # Perform a 4x4 tiling of the primitive cell
    g44 = g11.tile(4,4,1)

    # Add a Gamma-centered 2x2 Monkhorst-Pack grid
    g44g = g44.copy()
    g11g = g44g.folded_structure

    g44g.add_kmesh(kgrid=(2,2,1),kshift=(0,0,0))
    assert(structure_same(g44g,g44g_gen))
    assert(structure_same(g11g,g44g_gen.folded_structure))
    assert(len(g44g.kpoints)==4)
    assert(len(g11g.kpoints)==64)
    assert(value_eq(g44g.kpoints_unit(),g44g_ukp_ref))
    assert(value_eq(g11g.kpoints_unit(),g11g_ukp_ref))

    # Add a shifted 2x2 Monkhorst-Pack grid
    g44s = g44.copy()
    g11s = g44s.folded_structure

    g44s.add_kmesh(kgrid=(2,2,1),kshift=(0.5,0.5,0))
    assert(structure_same(g44s,g44s_gen))
    assert(structure_same(g11s,g44s_gen.folded_structure))
    assert(len(g44s.kpoints)==4)
    assert(len(g11s.kpoints)==64)
    assert(value_eq(g44s.kpoints_unit(),g44s_ukp_ref))
    assert(value_eq(g11s.kpoints_unit(),g11s_ukp_ref))

    # Get the mapping between supercell and primitive cell k-points
    kmap_ref = obj({
        0 : set([0,32,4,48,8,60,12,44,16,40,20,56,24,52,28,36]),
        1 : set([1,61,5,49,9,13,45,17,37,25,53,41,57,33,29,21]),
        2 : set([2,50,54,6,26,10,62,34,14,18,46,22,58,38,42,30]),
        3 : set([35,3,51,7,63,23,47,15,27,43,19,11,55,39,59,31]),
        })

    kmap = g44s.kmap()
    assert(object_eq(kmap,kmap_ref))

#end def test_monkhorst_pack_kpoints



if versions.spglib_available:
    def test_symm_kpoints():
        """
        Add symmetrized Monkhorst-Pack kpoints.
        """
        from structure import generate_structure

        # Note: this demo requires spglib

        g44 = generate_structure(
            structure  = 'graphene',
            cell       = 'prim',
            tiling     = (4,4,1),
            kgrid      = (4,4,1),
            kshift     = (0,0,0),
            symm_kgrid = True,
            )

        g11 = g44.folded_structure

        g44_kw_ref = np.array([1,6,3,6],dtype=float)

        g44_ukp_ref = np.array([
            [ 0.00, 0.00, 0.00 ],
            [ 0.25, 0.00, 0.00 ],
            [ 0.50, 0.00, 0.00 ],
            [ 0.25, 0.25, 0.00 ]])

        g11_ukp_ref = np.array([
            [ 0.0000, 0.0000, 0.0000 ],
            [ 0.0625, 0.0000, 0.0000 ],
            [ 0.1250, 0.0000, 0.0000 ],
            [ 0.0625, 0.0625, 0.0000 ],
            [ 0.2500, 0.0000, 0.0000 ],
            [ 0.3125, 0.0000, 0.0000 ],
            [ 0.3750, 0.0000, 0.0000 ],
            [ 0.3125, 0.0625, 0.0000 ],
            [ 0.5000, 0.0000, 0.0000 ],
            [ 0.5625, 0.0000, 0.0000 ],
            [ 0.6250, 0.0000, 0.0000 ],
            [ 0.5625, 0.0625, 0.0000 ],
            [ 0.7500, 0.0000, 0.0000 ],
            [ 0.8125, 0.0000, 0.0000 ],
            [ 0.8750, 0.0000, 0.0000 ],
            [ 0.8125, 0.0625, 0.0000 ],
            [ 0.0000, 0.2500, 0.0000 ],
            [ 0.0625, 0.2500, 0.0000 ],
            [ 0.1250, 0.2500, 0.0000 ],
            [ 0.0625, 0.3125, 0.0000 ],
            [ 0.2500, 0.2500, 0.0000 ],
            [ 0.3125, 0.2500, 0.0000 ],
            [ 0.3750, 0.2500, 0.0000 ],
            [ 0.3125, 0.3125, 0.0000 ],
            [ 0.5000, 0.2500, 0.0000 ],
            [ 0.5625, 0.2500, 0.0000 ],
            [ 0.6250, 0.2500, 0.0000 ],
            [ 0.5625, 0.3125, 0.0000 ],
            [ 0.7500, 0.2500, 0.0000 ],
            [ 0.8125, 0.2500, 0.0000 ],
            [ 0.8750, 0.2500, 0.0000 ],
            [ 0.8125, 0.3125, 0.0000 ],
            [ 0.0000, 0.5000, 0.0000 ],
            [ 0.0625, 0.5000, 0.0000 ],
            [ 0.1250, 0.5000, 0.0000 ],
            [ 0.0625, 0.5625, 0.0000 ],
            [ 0.2500, 0.5000, 0.0000 ],
            [ 0.3125, 0.5000, 0.0000 ],
            [ 0.3750, 0.5000, 0.0000 ],
            [ 0.3125, 0.5625, 0.0000 ],
            [ 0.5000, 0.5000, 0.0000 ],
            [ 0.5625, 0.5000, 0.0000 ],
            [ 0.6250, 0.5000, 0.0000 ],
            [ 0.5625, 0.5625, 0.0000 ],
            [ 0.7500, 0.5000, 0.0000 ],
            [ 0.8125, 0.5000, 0.0000 ],
            [ 0.8750, 0.5000, 0.0000 ],
            [ 0.8125, 0.5625, 0.0000 ],
            [ 0.0000, 0.7500, 0.0000 ],
            [ 0.0625, 0.7500, 0.0000 ],
            [ 0.1250, 0.7500, 0.0000 ],
            [ 0.0625, 0.8125, 0.0000 ],
            [ 0.2500, 0.7500, 0.0000 ],
            [ 0.3125, 0.7500, 0.0000 ],
            [ 0.3750, 0.7500, 0.0000 ],
            [ 0.3125, 0.8125, 0.0000 ],
            [ 0.5000, 0.7500, 0.0000 ],
            [ 0.5625, 0.7500, 0.0000 ],
            [ 0.6250, 0.7500, 0.0000 ],
            [ 0.5625, 0.8125, 0.0000 ],
            [ 0.7500, 0.7500, 0.0000 ],
            [ 0.8125, 0.7500, 0.0000 ],
            [ 0.8750, 0.7500, 0.0000 ],
            [ 0.8125, 0.8125, 0.0000 ]])

        assert(value_eq(g44.kweights,g44_kw_ref))
        assert(value_eq(g44.kpoints_unit(),g44_ukp_ref))
        assert(value_eq(g11.kpoints_unit(),g11_ukp_ref))

    #end def test_symm_kpoints
#end if


def test_count_kshells():
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



def test_rwigner():
    from structure import generate_structure
    s = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        tiling    = (4,4,1),
        )
    rw = s.rwigner()
    assert(value_eq(float(rw),4.924))
#end def test_rwigner



def test_volume():
    gen = get_generated_structures()
    d64 = gen.diamond_64
    assert(value_eq(d64.volume(),363.994344))
#end def test_volume



if versions.scipy_available:
    def test_madelung():
        gen = get_generated_structures()
        d64 = gen.diamond_64
        assert(value_eq(d64.madelung(),-0.210284756321))
    #end def test_madelung



    def test_makov_payne():
        gen = get_generated_structures()
        d64 = gen.diamond_64
        assert(value_eq(d64.makov_payne(q=1,eps=5.68),0.0185109820705))
        assert(value_eq(d64.makov_payne(q=2,eps=5.68),0.074043928282))
    #end def test_makov_payne
#end if



def test_min_image_distances():
    """
    Compute minimum image distances between nearest neighbors.
    """
    import numpy as np
    from structure import generate_structure

    g = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        tiling    = (4,4,1),
        )

    # Get the neighbor (index) table, along with sorted distance 
    # and displacement tables.
    nt,dt,vt = g.neighbor_table(distances=True,vectors=True)

    # Restrict to 3 nearest neighbors (not including self at index 0)
    nt = nt[:,1:4]
    dt = dt[:,1:4]
    vt = vt[:,1:4]
    
    nt_ref = np.array([
        [ 1, 11,  9],
        [30,  0, 24],
        [ 3, 11, 13],
        [ 2, 26, 24],
        [ 5, 15, 13],
        [ 4, 26, 28],
        [ 7,  9, 15],
        [ 6, 30, 28],
        [ 9, 19, 17],
        [ 8,  6,  0],
        [11, 19, 21],
        [10,  0,  2],
        [13, 23, 21],
        [12,  2,  4],
        [15, 17, 23],
        [14,  4,  6],
        [17, 27, 25],
        [16,  8, 14],
        [19, 27, 29],
        [18, 10,  8],
        [21, 29, 31],
        [20, 12, 10],
        [23, 25, 31],
        [22, 12, 14],
        [25,  1,  3],
        [24, 22, 16],
        [27,  3,  5],
        [26, 18, 16],
        [29,  5,  7],
        [28, 18, 20],
        [31,  1,  7],
        [30, 22, 20]])

    for nti,nti_ref in zip(nt,nt_ref):
        assert(set(nti)==set(nti_ref))
    #end for

    dist = dt.ravel()
    assert(value_eq(dist.min(),1.42143636))
    assert(value_eq(dist.max(),1.42143636))

    vec = vt.ravel()
    vec.shape = np.prod(dt.shape),3
    vdist = np.linalg.norm(vec,axis=1)
    assert(value_eq(vdist,dist))

#end def test_min_image_distances



def test_freeze():
    """
    Freeze sets of atoms to prevent relaxation.
    """
    from structure import generate_structure

    g = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        tiling    = (6,6,1),
        )
    g.recenter((0,0,0))

    # Select atoms to be frozen
    outer_atoms = np.sqrt((g.pos**2).sum(1)) > 3.2

    # Freeze the atoms
    g.freeze(outer_atoms)

    # Get a mask array identifying the frozen atoms
    frozen = g.is_frozen()

    assert((frozen==outer_atoms).all())

#end def test_freeze



if versions.scipy_available:
    def test_embed():
        """
        Embed a "relaxed" structure in a larger pristine cell.
        """
        import numpy as np
        from structure import generate_structure

        center = (0,0,0)

        g = generate_structure(
            structure = 'graphene',
            cell      = 'prim',
            tiling    = (4,4,1),
            )
        g.recenter(center)

        # Represent the "relaxed" cell
        gr = g.copy()
        npos = len(gr.pos)
        dr = gr.min_image_vectors(center)
        dr.shape = npos,3
        r = np.linalg.norm(dr,axis=1)
        dilation = 2*r*np.exp(-r)
        for i in range(npos):
            if r[i]>0:
                gr.pos[i] += dilation[i]/r[i]*dr[i]
            #end if
        #end for

        # Represent the unrelaxed large cell
        gl = generate_structure(
            structure = 'graphene',
            cell      = 'rect',
            tiling    = (8,4,1),
            )
        gl.recenter(center)

        # Embed the relaxed cell in the large unrelaxed cell
        ge = gl.copy()
        ge.embed(gr)

        assert(len(ge.elem)==len(gl.elem))
        assert(len(ge.pos)==len(gl.pos))

        # check that the large local distortion made in the small cell
        # is present in the large cell after embedding
        rnn_max_ref = 2.1076122431022664
        rnn_max = np.linalg.norm(gr.pos[1]-gr.pos[30])
        assert(value_eq(rnn_max,rnn_max_ref))
        rnn_max = np.linalg.norm(ge.pos[1]-ge.pos[28])
        assert(value_eq(rnn_max,rnn_max_ref))

    #end test_embed
#end if



def test_interpolate():
    """
    Interpolate between two "relaxed" structures for NEB initialization.
    """
    import numpy as np
    from structure import generate_structure
    
    g = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        tiling    = (4,4,1),
        )
    g11 = g.folded_structure

    # Make chromium atom positions (pretend like they are above the surface)
    npos = np.dot((1./3,2./3,0),g11.axes)
    npos1 = npos+3*g11.axes[0]
    npos2 = npos1+g11.axes[0]+g11.axes[1]

    # "Relaxed" structure with additional atom on one ring
    gr1 = g.copy()
    gr1.add_atoms(['Cr'],[npos1])

    # "Relaxed" structure with additional atom on neighboring ring
    gr2 = g.copy()
    gr2.add_atoms(['Cr'],[npos2])
    gr2.recenter()

    # Interpolate between the two structures within min. image convention
    spath = gr1.interpolate(gr2,10)

    Cr_positions_ref = np.array([
        [ 7.386     ,  1.42143636,  0.        ],
        [ 7.49790909,  1.61526859,  0.        ],
        [ 7.60981818,  1.80910083,  0.        ],
        [ 7.72172727,  2.00293306,  0.        ],
        [ 7.83363636,  2.19676529,  0.        ],
        [ 7.94554545,  2.39059752,  0.        ],
        [ 8.05745455,  2.58442975,  0.        ],
        [ 8.16936364,  2.77826198,  0.        ],
        [-1.56672727,  2.97209421,  0.        ],
        [-1.45481818,  3.16592644,  0.        ],
        [-1.34290909,  3.35975868,  0.        ],
        [-1.231     ,  3.55359091,  0.        ]])

    Cr_positions = []
    for gr in spath:
        Cr_positions.append(gr.pos[-1])
    #end for
    Cr_positions = np.array(Cr_positions)

    assert(value_eq(Cr_positions,Cr_positions_ref))

#end def test_interpolate



if versions.spglib_available:
    def test_point_group_operations():
        from structure import generate_structure,Crystal

        nrotations = dict(
            Ca2CuO3    =  8,
            CaO        = 48,
            Cl2Ca2CuO2 = 16,
            CuO        =  2,
            CuO2_plane = 16,
            La2CuO4    =  2,
            NaCl       = 48,
            ZnO        =  6,
            calcium    = 48,
            copper     = 48,
            diamond    = 24,
            graphene   = 12,
            oxygen     =  4,
            rocksalt   = 48,
            wurtzite   =  6,
            )

        for struct,cell in sorted(Crystal.known_crystals.keys()):
            if cell!='prim':
                continue
            #end if

            s = generate_structure(
                structure = struct,
                cell      = cell,
                )
                
            rotations = s.point_group_operations()
            assert(struct in nrotations)
            assert(len(rotations)==nrotations[struct])

            valid = s.check_point_group_operations(rotations,exit=False)
            assert(valid)
        #end for

    #end def test_point_group_operations
#end if



if versions.spglib_available and versions.seekpath_available:
    def test_rmg_transform():
        from numpy import array
        from generic import obj
        from structure import generate_structure

        ref = obj({
            ('Ca2CuO3', 'conv') : obj(
                R          = array([[ 8.62068966e-01,  0.00000000e+00,  0.00000000e+00],
                                    [-5.27865000e-17,  1.16000000e+00,  0.00000000e+00],
                                    [-5.27865000e-17, -7.10295144e-17,  1.00000000e+00]]),
                bv         = 'orthorhombic_P',
                tmatrix    = None,
                rmg_inputs = obj(
                    a_length             = 3.2500000000000004,
                    b_length             = 3.7700000000000005,
                    bravais_lattice_type = 'Orthorhombic Primitive',
                    c_length             = 12.23,
                    length_units         = 'Angstrom',
                    ),
                ),
            ('Ca2CuO3', 'prim') : obj(
                R          = array([[ 1.38777878e-17,  1.00000000e+00, -7.13693767e-18],
                                    [ 2.77555756e-17,  3.61249492e-17,  3.76307692e+00],
                                    [ 2.65739984e-01, -5.79633098e-18, -2.39520632e-17]]),
                bv         = 'orthorhombic_P',
                tmatrix    = array([[0, 1, 1],
                                    [1, 0, 1],
                                    [1, 1, 0]]),
                rmg_inputs = obj(
                    a_length             = 3.2500000000000004,
                    b_length             = 3.77,
                    bravais_lattice_type = 'Orthorhombic Primitive',
                    c_length             = 12.23,
                    length_units         = 'Angstrom',
                    ),
                ),
            ('CaO', 'conv') : obj(
                R          = array([[ 1.000000e+00,  0.000000e+00,  0.000000e+00],
                                    [-6.123234e-17,  1.000000e+00,  0.000000e+00],
                                    [-6.123234e-17, -6.123234e-17,  1.000000e+00]]),
                bv         = 'cubic_P',
                tmatrix    = None,
                rmg_inputs = obj(
                    a_length             = 4.81,
                    b_length             = 4.81,
                    bravais_lattice_type = 'Cubic Primitive',
                    c_length             = 4.81,
                    length_units         = 'Angstrom',
                    ),
                ),
            ('CaO', 'prim') : obj(
                R          = array([[-1.28679696e-17,  1.00000000e+00,  1.28679696e-17],
                                    [ 1.28679696e-17,  1.28679696e-17,  1.00000000e+00],
                                    [ 1.00000000e+00, -1.28679696e-17, -1.28679696e-17]]),
                bv         = 'cubic_F',
                tmatrix    = None,
                rmg_inputs = obj(
                    a_length             = 4.81,
                    b_length             = 4.81,
                    bravais_lattice_type = 'Cubic Face Centered',
                    c_length             = 4.81,
                    length_units         = 'Angstrom',
                    ),
                ),
            ('Cl2Ca2CuO2', 'afm') : obj(
                R          = array([[ 0.70710678,  0.70710678,  0.        ],
                                    [-0.70710678,  0.70710678,  0.        ],
                                    [ 0.        ,  0.        ,  1.        ]]),
                bv         = 'tetragonal_P',
                tmatrix    = None,
                rmg_inputs = obj(
                    a_length             = 5.471592272821505,
                    b_length             = 5.471592272821505,
                    bravais_lattice_type = 'Tetragonal Primitive',
                    c_length             = 15.049999999999999,
                    length_units         = 'Angstrom',
                    ),
                ),
            ('Cl2Ca2CuO2', 'prim') : obj(
                R          = array([[-1.99922127e-16,  1.00000000e+00,  1.13566145e-16],
                                    [ 1.84864975e-17, -1.84864975e-17,  3.88989403e+00],
                                    [ 2.57076412e-01,  1.87078112e-18, -1.25038407e-17]]),
                bv         = 'tetragonal_P',
                tmatrix    = array([[0, 1, 1],
                                    [1, 0, 1],
                                    [1, 1, 0]]),
                rmg_inputs = obj(
                    a_length             = 3.8689999999999998,
                    b_length             = 3.869,
                    bravais_lattice_type = 'Tetragonal Primitive',
                    c_length             = 15.05,
                    length_units         = 'Angstrom',
                    ),
                ),
            ('CuO', 'conv') : obj(
                R          = None,
                bv         = 'monoclinic_P',
                tmatrix    = None,
                rmg_inputs = obj(
                    ),
                ),
            ('ZnO', 'conv') : obj(
                R          = array([[ 1.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                    [-2.51021563e-16,  1.00000000e+00,  0.00000000e+00],
                                    [-6.12323400e-17, -1.06057524e-16,  1.00000000e+00]]),
                bv         = 'hexagonal_P',
                tmatrix    = None,
                rmg_inputs = obj(
                    a_length             = 3.349999999999999,
                    b_length             = 3.349999999999999,
                    bravais_lattice_type = 'Hexagonal Primitive',
                    c_length             = 5.22,
                    length_units         = 'Angstrom',
                    ),
                ),
            })

        res = obj()
        for struct,cell in ref.keys():
             s = generate_structure(
                 structure = struct,
                 cell      = cell,
                 )
             st,rmg_inputs,R,tmatrix,bv = s.rmg_transform(
                 allow_tile    = True,
                 allow_general = True,
                 all_results   = True,
                 )
             res[struct,cell] = obj(
                 rmg_inputs = rmg_inputs,
                 R          = R,
                 tmatrix    = tmatrix,
                 bv         = bv,
                 )
        #end for

        assert(testing.check_object_eq(res,ref,atol=1e-12))
    #end def test_rmg_transform
#end if
