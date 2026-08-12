import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.VASP_ANALYZER)

from ..generic import generic_settings
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
        cls.check_consistency()
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
    assert(reread.lglobal)
    assert(reread.spring==-5.0)
#end def test_current_keywords_roundtrip


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
