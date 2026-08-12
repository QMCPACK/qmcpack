import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.VASP_ANALYZER)

from ..generic import generic_settings,NexusError
generic_settings.raise_error = True

from nexus.nexus_base import nexus_core
from . import isolate_nexus_core, TEST_DIR
from .. import testing
from ..testing import object_eq,dict_serialize


serial_references = dict()

c_potcar_text = '''
This is not a real POTCAR file.

End of Dataset
'''

def generate_serial_references():
    import numpy as np
    from ..developer import obj

    ref = {
        'incar/encut' : 450.0,
        'incar/ibrion' : 2,
        'incar/icharg' : 2,
        'incar/isif' : 2,
        'incar/istart' : 0,
        'incar/nelect' : 64,
        'incar/nsw' : 5,
        'kpoints/centering' : 'monkhorst',
        'kpoints/kgrid' : (2, 2, 2),
        'kpoints/kshift' : (0, 0, 0),
        'kpoints/mode' : 'auto',
        'poscar/axes' : np.array([
            [3.57, 3.57, 0.0],
            [0.0, 3.57, 3.57],
            [3.57, 0.0, 3.57]]),
        'poscar/coord' : 'cartesian',
        'poscar/description' : None,
        'poscar/dynamic' : None,
        'poscar/elem' : ['C'],
        'poscar/elem_count' : [16],
        'poscar/pos' : np.array([
            [0.0, 0.0, 0.0],
            [0.8925, 0.8925, 0.8925],
            [1.785, 1.785, 0.0],
            [2.6775, 2.6775, 0.8925],
            [0.0, 1.785, 1.785],
            [0.8925, 2.6775, 2.6775],
            [1.785, 3.57, 1.785],
            [2.6775, 4.4625, 2.6775],
            [1.785, 0.0, 1.785],
            [2.6775, 0.8925, 2.6775],
            [3.57, 1.785, 1.785],
            [4.4625, 2.6775, 2.6775],
            [1.785, 1.785, 3.57],
            [2.6775, 2.6775, 4.4625],
            [3.57, 3.57, 3.57],
            [4.4625, 4.4625, 4.4625]]),
        'poscar/scale' : 1.0,
        'poscar/vel' : None,
        'poscar/vel_coord' : None,
        'potcar/files' : ['C.POTCAR'],
        'potcar/pseudos' : obj(),
        }
    
    ref_read = ref.copy()

    ref_read['incar/nelect']       = 64.0
    ref_read['kpoints/centering']  = 'monkhorst-pack'
    ref_read['kpoints/kgrid']      = np.array((2, 2, 2),dtype=int)
    ref_read['kpoints/kshift']     = np.array((0, 0, 0),dtype=float)
    ref_read['poscar/description'] = 'System cell and coordinates'
    ref_read['poscar/elem']        = np.array(['C'],dtype=str)
    ref_read['poscar/elem_count']  = np.array([16],dtype=int)
    ref_read['potcar/files']       = None
    del ref_read['potcar/pseudos']
    ref_read['potcar/pseudos/0']   = c_potcar_text

    serial_references['read']     = ref_read

    serial_references['write']    = ref_read.copy()

    serial_references['generate'] = ref.copy()

#end def generate_serial_references


def get_serial_references():
    if len(serial_references)==0:
        generate_serial_references()
    #end if
    return serial_references
#end def get_serial_references


def check_vs_serial_reference(gi,name):
    from ..developer import obj
    sr = obj(get_serial_references()[name])
    sg = dict_serialize(gi,dict_type=obj)
    same = object_eq(sg,sr)
    if not same:
        print('\n'+name+' differs')
        testing.print_diff(sr,sg)
    #end if
    assert(same)
#end def check_vs_serial_reference


TEST_FILES = {
    "d16bulk.POSCAR":  TEST_DIR / "test_vasp_input_files/d16bulk.POSCAR",
    "diamond_INCAR":   TEST_DIR / "test_vasp_input_files/diamond_INCAR",
    "diamond_KPOINTS": TEST_DIR / "test_vasp_input_files/diamond_KPOINTS",
    "diamond_POSCAR":  TEST_DIR / "test_vasp_input_files/diamond_POSCAR",
    "diamond_POTCAR":  TEST_DIR / "test_vasp_input_files/diamond_POTCAR",
    }

for file in TEST_FILES.values():
    assert(file.exists()), f"Test file not found! {file}"


def test_keyword_consistency():
    from ..vasp_input import Incar,Stopcar

    for cls in Incar,Stopcar:
        scalar_keywords = set()
        array_keywords = set()
        typed_keywords = set()
        for keyword_type in cls.kw_scalars+cls.kw_arrays:
            keywords = getattr(cls,keyword_type)
            assert(typed_keywords.isdisjoint(keywords))
            typed_keywords |= keywords
            if keyword_type in cls.kw_scalars:
                scalar_keywords |= keywords
            else:
                array_keywords |= keywords
            #end if
        #end for
        assert(scalar_keywords==cls.scalar_keywords)
        assert(array_keywords==cls.array_keywords)
        assert(typed_keywords==cls.keywords)
        assert(set(cls.type.keys())==cls.keywords)
        assert(set(cls.read_value.keys())==cls.keywords)
        assert(set(cls.write_value.keys())==cls.keywords)
        assert(set(cls.assign_value.keys())==cls.keywords)
    #end for
#end def test_keyword_consistency


def test_current_keyword_types():
    from ..vasp_input import Incar

    expected = {
        'ints': {
            'efermi_nedos', 'elph_nbands', 'iopt', 'libmbd_n_omega_grid',
            'ml_output_mode',
            },
        'reals': {
            'apaco', 'ddr', 'elph_kspacing', 'nupdown', 'spring',
            },
        'bools': {
            'elph_run', 'lglobal', 'ml_lib',
            },
        'strings': {
            'efermi', 'elph_driver', 'ml_mode', 'pyamff_model', 'xc',
            },
        'int_arrays': {
            'elph_nbands_sum', 'libmbd_k_grid',
            },
        'real_arrays': {
            'elph_selfen_temps', 'pomass', 'quad_efg', 'smearings',
            },
        'bool_arrays': {
            'fmp_active',
            },
        }

    for keyword_type,keywords in expected.items():
        for keyword in keywords:
            assert(Incar.type[keyword]==keyword_type)
        #end for
    #end for
#end def test_current_keyword_types


def test_current_keywords_roundtrip():
    import numpy as np
    from ..vasp_input import Incar

    incar = Incar()
    incar.assign(
        elph_driver       = 'el',
        elph_nbands       = 8,
        elph_run          = True,
        elph_selfen_temps = [0.0,300.0],
        fmp_active        = [True,False],
        iopt              = 7,
        libmbd_k_grid     = [2,2,1],
        lglobal           = True,
        spring            = -5.0,
        )

    reread = Incar()
    reread.read_text(incar.write_text())

    assert(reread.elph_driver=='el')
    assert(reread.elph_nbands==8)
    assert(reread.elph_run)
    assert(np.array_equal(reread.elph_selfen_temps,[0.0,300.0]))
    assert(np.array_equal(reread.fmp_active,[True,False]))
    assert(reread.iopt==7)
    assert(np.array_equal(reread.libmbd_k_grid,[2,2,1]))
    assert(reread.lglobal)
    assert(reread.spring==-5.0)
#end def test_current_keywords_roundtrip


def test_integer_assignment_validation():
    import numpy as np
    from ..vasp_input import Incar,generate_vasp_input

    incar = Incar()
    incar.assign(
        ibrion      = np.int64(2),
        images      = 3.0+5e-9,
        random_seed = [1,2.0-5e-9,np.int32(3)],
        )
    assert(incar.ibrion==2)
    assert(incar.images==3)
    assert(np.array_equal(incar.random_seed,[1,2,3]))

    generated = generate_vasp_input(ibrion=2.0+5e-9)
    assert(generated.incar.ibrion==2)

    invalid_scalars = (
        True,
        '3',
        2.0+2e-8,
        2.1,
        np.nan,
        np.inf,
        float(2**54),
        )
    for value in invalid_scalars:
        with pytest.raises(NexusError,match='assign failed for keyword ibrion'):
            Incar().assign(ibrion=value)
        #end with
    #end for

    invalid_arrays = (
        [1,2.1,3],
        [1,True,3],
        [1,'2',3],
        [[1,2],[3,4]],
        )
    for value in invalid_arrays:
        with pytest.raises(
            NexusError,
            match='assign failed for keyword random_seed',
            ):
            Incar().assign(random_seed=value)
        #end with
    #end for

    with pytest.raises(NexusError,match='element at index 1'):
        Incar().assign(random_seed=[1,2.1,3])
    #end with
    with pytest.raises(NexusError,match='assign failed for keyword ibrion'):
        generate_vasp_input(ibrion=2.1)
    #end with
#end def test_integer_assignment_validation


def test_keyword_syntax_parsing():
    import numpy as np
    from ..vasp_input import Incar

    text = (
        'ENCUT = 4D2; ISMEAR = -1 ! trailing comment\n'
        'MAGMOM = 2*1.5 \\\n'
        '         2*-1.5 # another comment\n'
        'LCHARG = .false.\n'
        'IBRION = 1; IBRION = 2\n'
        )
    incar = Incar()
    incar.read_text(text)

    assert(incar.encut==400.0)
    assert(incar.ismear==-1)
    assert(np.array_equal(incar.magmom,[1.5,1.5,-1.5,-1.5]))
    assert(not incar.lcharg)
    assert(incar.ibrion==2)
#end def test_keyword_syntax_parsing


def test_quoted_string_preservation():
    from ..vasp_input import Incar

    text = (
        'SYSTEM = "a; b # c ! d = e"; '
        'WANNIER90_WIN = "Begin; x=y # z!"\n'
        'ENCUT = 400\n'
        )
    incar = Incar()
    incar.read_text(text)

    assert(incar.system=='a; b # c ! d = e')
    assert(incar.wannier90_win=='Begin; x=y # z!')
    assert(incar.encut==400.0)
#end def test_quoted_string_preservation


def test_array_compression_roundtrip():
    import numpy as np
    from ..vasp_input import Incar

    incar = Incar()
    incar.assign(
        lattice_constraints = [True]*4,
        magmom               = [1.5]*4,
        random_seed          = [2]*4,
        )
    text = incar.write_text()

    assert('4*T' in text)
    assert('4*1.5' in text)
    assert('4*2' in text)

    reread = Incar()
    reread.read_text(text)
    assert(np.array_equal(reread.lattice_constraints,[True]*4))
    assert(np.array_equal(reread.magmom,[1.5]*4))
    assert(np.array_equal(reread.random_seed,[2]*4))
#end def test_array_compression_roundtrip


def test_invalid_keyword_input():
    from ..vasp_input import Incar

    with pytest.raises(NexusError,match='NOT_A_TAG is not a keyword'):
        Incar().assign(not_a_tag=1)
    #end with
    with pytest.raises(NexusError,match='read failed for keyword lcharg'):
        Incar().read_text('LCHARG = maybe')
    #end with
    with pytest.raises(NexusError,match='assign failed for keyword magmom'):
        Incar().assign(magmom=1.0)
    #end with
    with pytest.raises(NexusError,match='quotation marks.*not paired'):
        Incar().read_text('SYSTEM = "unfinished')
    #end with
    with pytest.raises(NexusError,match='incomplete line continuation'):
        Incar().read_text('ENCUT = 400 \\')
    #end with
    with pytest.raises(NexusError,match='repeat count must be non-negative'):
        Incar().read_text('MAGMOM = -2*1.0')
    #end with
#end def test_invalid_keyword_input


def test_keyword_catalog_integrity():
    from ..vasp_input import Incar

    keyword_sets = (
        Incar.ints,
        Incar.reals,
        Incar.bools,
        Incar.strings,
        Incar.int_arrays,
        Incar.real_arrays,
        Incar.bool_arrays,
        )
    assert(len(Incar.keywords)==684)
    assert(sum(map(len,keyword_sets))==len(Incar.keywords))
    assert(Incar.deprecated<=Incar.keywords)
    assert(Incar.unsupported.isdisjoint(Incar.keywords))
#end def test_keyword_catalog_integrity


def test_deprecated_keyword_compatibility():
    from ..vasp_input import Incar

    incar = Incar()
    incar.assign(
        ichain     = 0,
        jacobian   = 1.0,
        lclimb     = True,
        ldneb      = False,
        timestep   = 0.1,
        )
    reread = Incar()
    reread.read_text(incar.write_text())

    assert(reread.ichain==0)
    assert(reread.jacobian==1.0)
    assert(reread.lclimb)
    assert(not reread.ldneb)
    assert(reread.timestep==0.1)
#end def test_deprecated_keyword_compatibility


def test_keyword_file_roundtrip():
    import numpy as np
    from ..vasp_input import Incar

    text = 'SYSTEM = "first line\nsecond line"\nENCUT = 400\nLATTICE_CONSTRAINTS = t f .TRUE.\n'
    incar = Incar()
    incar.read_text(text)

    assert(incar.system=='first line\nsecond line')
    assert(incar.encut==400.0)
    assert(np.array_equal(incar.lattice_constraints,(True,False,True)))

    reread = Incar()
    reread.read_text(incar.write_text())
    assert(reread.system==incar.system)
    assert(reread.encut==incar.encut)
#end def test_keyword_file_roundtrip


def test_kpoints_tetrahedra_roundtrip():
    import numpy as np
    from ..developer import obj
    from ..vasp_input import Kpoints

    kpoints = Kpoints()
    kpoints.mode = 'explicit'
    kpoints.coord = 'reciprocal'
    kpoints.kpoints = np.array([[0.0,0.0,0.0]])
    kpoints.kweights = np.array([1.0])
    kpoints.tetrahedra = obj()
    kpoints.tetrahedra[0] = obj(
        volume     = 0.5,
        degeneracy = 1,
        corners    = np.array([1,1,1,1]),
        )

    reread = Kpoints()
    reread.read_text(kpoints.write_text())

    assert(reread.mode=='explicit')
    assert(len(reread.tetrahedra)==1)
    assert(reread.tetrahedra[0].volume==0.5)
    assert(np.array_equal(reread.tetrahedra[0].corners,(1,1,1,1)))
#end def test_kpoints_tetrahedra_roundtrip


def test_poscar_coordinates_and_structure_conversion():
    import numpy as np
    from ..vasp_input import Poscar,VaspInput

    poscar = Poscar()
    poscar.scale = 2.0
    poscar.axes = np.eye(3)
    poscar.elem = np.array(['H'])
    poscar.elem_count = np.array([1])
    poscar.coord = 'direct'
    poscar.pos = np.array([[0.25,0.25,0.25]])

    vasp_input = VaspInput()
    vasp_input.poscar = poscar
    poscar.change_specifier('cartesian',vasp_input)
    poscar.change_specifier('direct',vasp_input)
    assert(np.allclose(poscar.pos,[[0.25,0.25,0.25]]))
    assert(poscar.check_complete()=='')

    structure = vasp_input.return_system(structure_only=True)
    assert(np.allclose(structure.axes,2.0*np.eye(3)))
    assert(np.allclose(structure.pos,[[0.5,0.5,0.5]]))

    selective_text = '''selective dynamics
1.0
1 0 0
0 1 0
0 0 1
H
1
Selective dynamics
Direct
0 0 0 t f t
'''
    selective = Poscar()
    selective.read_text(selective_text)
    assert(np.array_equal(selective.dynamic,[[True,False,True]]))
#end def test_poscar_coordinates_and_structure_conversion


def test_potcar_construction_and_metadata(tmp_path):
    from ..vasp_input import Potcar

    empty = Potcar()
    assert(empty.write_text()=='')

    potcar_text = 'PAW_PBE H 01Jan2001\n1.0\nVRHFIN =H:\nEnd of Dataset'
    potcar_file = tmp_path / 'H.POTCAR'
    potcar_file.write_text(potcar_text)

    potcar = Potcar(potcar_file)
    assert(len(potcar.pseudos)==1)
    assert(potcar.write_text()==potcar_text)

    info = potcar.pot_info()
    assert(len(info)==1)
    assert(info[0].element=='H')
    assert(info[0].Zval==1)

    configured = Potcar(tmp_path,['H.POTCAR'])
    configured.load()
    assert(configured.pseudos[0]==potcar_text)
#end def test_potcar_construction_and_metadata


def test_input_filename_filtering(tmp_path):
    from ..vasp_input import VaspInput

    (tmp_path / 'INCAR').write_text('ENCUT = 100\n')
    (tmp_path / 'job_INCAR.in').write_text('ENCUT = 200\n')

    vasp_input = VaspInput(tmp_path,prefix='job_',postfix='.in')
    assert(vasp_input.incar.encut==200.0)
#end def test_input_filename_filtering



def test_empty_init():
    from ..vasp_input import VaspInput
    vi = VaspInput()
#end test_empty_init



def test_read():
    from ..vasp_input import VaspInput

    test_files_dir = TEST_FILES['diamond_INCAR'].parent
    vi = VaspInput(test_files_dir, prefix='diamond_')

    del vi.potcar.filepath

    check_vs_serial_reference(vi,'read')
#end test_read



def test_write(tmp_path):
    from ..vasp_input import VaspInput

    test_files_dir = TEST_FILES['diamond_INCAR'].parent

    vi_read = VaspInput(test_files_dir, prefix='diamond_')
    
    vi_read.write(tmp_path, prefix='write_diamond_')

    vi_write = VaspInput(tmp_path, prefix='write_diamond_')

    del vi_write.potcar.filepath

    check_vs_serial_reference(vi_write,'write')
#end test_write


@isolate_nexus_core
def test_generate(tmp_path):
    from ..nexus_base import nexus_noncore
    from ..physical_system import generate_physical_system
    from ..vasp_input import generate_vasp_input,VaspInput

    pseudo_dir = tmp_path / 'pseudopotentials'
    pseudo_dir.mkdir()

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]
    nexus_noncore.pseudo_dir = pseudo_dir

    (pseudo_dir / 'C.POTCAR').write_text(c_potcar_text)


    dia16 = generate_physical_system(
        structure = TEST_FILES['d16bulk.POSCAR'],
        C         = 4                  
        )

    vi = generate_vasp_input(      
        system       = dia16,            
        pseudos      = ['C.POTCAR'], 
        input_type   = 'generic',
        istart       = 0, 
        icharg       = 2,
        encut        = 450,
        nsw          = 5,
        ibrion       = 2,
        isif         = 2,
        kcenter      = 'monkhorst',
        kgrid        = (2,2,2),                
        kshift       = (0,0,0),              
        )

    assert(isinstance(vi,VaspInput))

    del vi.potcar.filepath

    check_vs_serial_reference(vi,'generate')

    vi_minimal = generate_vasp_input(system=dia16)
    assert(vi_minimal.incar.nelect==64)

    vi_title = generate_vasp_input(title='diamond')
    assert(vi_title.incar.system=='diamond')
#end def test_generate
