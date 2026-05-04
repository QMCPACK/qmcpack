import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.VASP_ANALYZER)

from ..generic import generic_settings
generic_settings.raise_error = True

from nexus.nexus_base import nexus_core
from . import isolate_nexus_core, TEST_DIR
from .. import testing
from ..testing import object_eq


def format_value(v):
    import numpy as np
    s = ''
    if isinstance(v,np.ndarray):
        pad = 12*' '
        s = 'np.array([\n'
        if len(v.shape)==1:
            s += pad
            for vv in v:
                s += format_value(vv)+','
            #end for
            s = s[:-1]
        else:
            for vv in v:
                s += pad + format_value(list(vv))+',\n'
            #end for
            s = s[:-2]
        #end if
        s += '])'
    elif isinstance(v,(str,np.bytes_)):
        s = "'"+str(v)+"'"
    else:
        s = str(v)
    #end if
    return s
#end def format_value


def make_serial_reference(gi):
    s = gi.serial()
    ref = '    ref = {\n'
    for k in sorted(s.keys()):
        v = s[k]
        ref +="        '{}' : {},\n".format(k,format_value(v))
    #end for
    ref += '        }\n'
    return ref
#end def make_serial_reference


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
    sg = gi.serial()
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
    nexus_core.remote_directory = tmp_path
    nexus_core.file_locations = nexus_core.file_locations + [tmp_path]
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
#end def test_generate
