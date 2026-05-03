import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.VASP_SIMULATION)

from ..generic import generic_settings
generic_settings.raise_error = True

from pathlib import Path
from . import isolate_nexus_core, create_pseudo_files
from nexus.nexus_base import nexus_core
from ..testing import clear_all_sims
from ..testing import failed,FailedTest
from ..testing import value_eq,object_eq,check_object_eq

from .test_vasp_input import c_potcar_text, TEST_FILES


def setup_vasp_sim(path,identifier='vasp',copy_files=False):
    import shutil
    from ..nexus_base import nexus_core
    from ..machines import job
    from ..physical_system import generate_physical_system
    from ..vasp import generate_vasp,Vasp

    nexus_core.runs = ''

    dia16 = generate_physical_system(
        structure = TEST_FILES['d16bulk.POSCAR'],
        C         = 4                  
        )

    sim = generate_vasp(
        identifier   = identifier,
        path         = path,
        job          = job(machine='ws1',cores=1),
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

    assert(isinstance(sim,Vasp))

    if copy_files:
        vfiles = [
            'diamond_INCAR',
            'diamond_KPOINTS',
            'diamond_POSCAR',
            'diamond_POTCAR',
            ]
        for vfile in vfiles:
            shutil.copy2(TEST_FILES[vfile],path)
        #end for
    #end if
        
    return sim
#end def setup_vasp_sim



def test_minimal_init():
    from ..machines import job
    from ..vasp import Vasp,generate_vasp

    sim = generate_vasp(
        job = job(machine='ws1',cores=1),
        )

    assert(isinstance(sim,Vasp))

    clear_all_sims()
#end def test_minimal_init


@isolate_nexus_core(needs_tmp_path=True)
def test_check_result(tmp_path):
    
    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = tmp_path
    nexus_core.file_locations = nexus_core.file_locations + [tmp_path]
    create_pseudo_files(
        tmp_dir=tmp_path,
        pseudos=["C.POTCAR"],
        pseudo_strs=[c_potcar_text],
    )

    sim = setup_vasp_sim(tmp_path)

    assert(not sim.check_result('unknown',None))

    assert(sim.check_result('structure',None))

    clear_all_sims()
#end def test_check_result


@isolate_nexus_core(needs_tmp_path=True)
def test_get_result(tmp_path):
    import shutil
    from numpy import array
    from ..developer import obj, NexusError

    create_pseudo_files(
        tmp_dir=tmp_path,
        pseudos=["C.POTCAR"],
        pseudo_strs=[c_potcar_text],
    )

    sim = setup_vasp_sim(tmp_path, identifier='diamond', copy_files=True)

    try:
        sim.get_result('unknown',None)
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    pcfile = tmp_path / 'diamond_POSCAR'
    ccfile = tmp_path / (sim.identifier+'.CONTCAR')

    assert(Path(sim.locdir).resolve()==tmp_path)

    shutil.copy2(pcfile,ccfile)

    assert(ccfile.exists())

    result = sim.get_result('structure',None)

    result_ref = obj(
          structure = obj(
              axes            = array(
                                [[3.57, 3.57, 0.  ],
                                 [0.  , 3.57, 3.57],
                                 [3.57, 0.  , 3.57]],dtype=float),
              background_charge = 0,
              bconds          = array([]),
              center          = array([3.57,3.57,3.57],dtype=float),
              dim             = 3,
              elem            = array(['C','C','C','C','C','C','C','C',
                                       'C','C','C','C','C','C','C','C'],
                                      dtype=object),
              folded_structure = None,
              frozen          = None,
              kaxes           = array(
                                [[ 0.87999794,  0.87999794, -0.87999794],
                                 [-0.87999794,  0.87999794,  0.87999794],
                                 [ 0.87999794, -0.87999794,  0.87999794]],
                                 dtype=float),
              kpoints         = array([],dtype=float),
              kweights        = array([],dtype=float),
              mag             = None,
              pos             = array(
                                [[0.    , 0.    , 0.    ],
                                 [0.8925, 0.8925, 0.8925],
                                 [1.785 , 1.785 , 0.    ],
                                 [2.6775, 2.6775, 0.8925],
                                 [0.    , 1.785 , 1.785 ],
                                 [0.8925, 2.6775, 2.6775],
                                 [1.785 , 3.57  , 1.785 ],
                                 [2.6775, 4.4625, 2.6775],
                                 [1.785 , 0.    , 1.785 ],
                                 [2.6775, 0.8925, 2.6775],
                                 [3.57  , 1.785 , 1.785 ],
                                 [4.4625, 2.6775, 2.6775],
                                 [1.785 , 1.785 , 3.57  ],
                                 [2.6775, 2.6775, 4.4625],
                                 [3.57  , 3.57  , 3.57  ],
                                 [4.4625, 4.4625, 4.4625]],
                                 dtype=float),
              scale           = 1.0,
              units           = 'A',
              ),
        )

    assert(check_object_eq(result,result_ref))

    clear_all_sims()
#end def test_get_result


@isolate_nexus_core(needs_tmp_path=True)
def test_incorporate_result(tmp_path):
    import shutil
    from numpy import array
    from ..developer import obj

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = tmp_path
    nexus_core.file_locations = nexus_core.file_locations + [tmp_path]
    create_pseudo_files(
        tmp_dir=tmp_path,
        pseudos=["C.POTCAR"],
        pseudo_strs=[c_potcar_text],
    )

    sim = setup_vasp_sim(tmp_path,identifier='diamond',copy_files=True)

    pcfile = tmp_path / 'diamond_POSCAR'
    ccfile = tmp_path / (sim.identifier+'.CONTCAR')

    assert(Path(sim.locdir).resolve()==tmp_path)

    shutil.copy2(pcfile,ccfile)

    assert(ccfile.exists())

    result = sim.get_result('structure',None)

    sim2 = setup_vasp_sim(tmp_path)

    pid = id(sim2.input.poscar)

    sim2.incorporate_result('structure',result,None)

    assert(id(sim2.input.poscar)!=pid)

    poscar_ref = obj(
        axes            = array(
                          [[3.57, 3.57, 0.  ],
                           [0.  , 3.57, 3.57],
                           [3.57, 0.  , 3.57]],dtype=float),
        coord           = 'cartesian',
        description     = None,
        dynamic         = None,
        elem            = ['C'],
        elem_count      = [16],
        pos             = array(
                          [[0.    , 0.    , 0.    ],
                           [0.8925, 0.8925, 0.8925],
                           [1.785 , 1.785 , 0.    ],
                           [2.6775, 2.6775, 0.8925],
                           [0.    , 1.785 , 1.785 ],
                           [0.8925, 2.6775, 2.6775],
                           [1.785 , 3.57  , 1.785 ],
                           [2.6775, 4.4625, 2.6775],
                           [1.785 , 0.    , 1.785 ],
                           [2.6775, 0.8925, 2.6775],
                           [3.57  , 1.785 , 1.785 ],
                           [4.4625, 2.6775, 2.6775],
                           [1.785 , 1.785 , 3.57  ],
                           [2.6775, 2.6775, 4.4625],
                           [3.57  , 3.57  , 3.57  ],
                           [4.4625, 4.4625, 4.4625]],dtype=float),
        scale           = 1.0,
        vel             = None,
        vel_coord       = None,
        )

    assert(object_eq(sim2.input.poscar.to_obj(),poscar_ref))

    clear_all_sims()
#end def test_incorporate_result


@isolate_nexus_core(needs_tmp_path=True)
def test_check_sim_status(tmp_path):

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = tmp_path
    nexus_core.file_locations = nexus_core.file_locations + [tmp_path]
    create_pseudo_files(
        tmp_dir=tmp_path,
        pseudos=["C.POTCAR"],
        pseudo_strs=[c_potcar_text],
    )

    sim = setup_vasp_sim(tmp_path)

    assert(Path(sim.locdir).resolve()==tmp_path)

    assert(not sim.finished)
    assert(not sim.input.performing_neb())

    sim.check_sim_status()

    assert(not sim.finished)

    outcar_file = tmp_path / 'OUTCAR'
    outcar_text = 'General timing and accounting'
    outcar_file.write_text(outcar_text)
    assert(outcar_file.exists())
    assert(outcar_text in outcar_file.read_text())

    sim.check_sim_status()

    assert(sim.finished)
    assert(not sim.failed)

    clear_all_sims()
#end def test_check_sim_status


@isolate_nexus_core(needs_tmp_path=True)
def test_get_output_files(tmp_path):

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = tmp_path
    nexus_core.file_locations = nexus_core.file_locations + [tmp_path]
    create_pseudo_files(
        tmp_dir=tmp_path,
        pseudos=["C.POTCAR"],
        pseudo_strs=[c_potcar_text],
    )

    sim = setup_vasp_sim(tmp_path)

    vfiles = 'INCAR KPOINTS POSCAR CONTCAR OUTCAR'.split()
    for vfile in vfiles:
        (tmp_path / vfile).touch()
    #end for

    files = sim.get_output_files()

    assert(value_eq(files,vfiles))

    for vfile in vfiles:
        assert((tmp_path / (sim.identifier+'.'+vfile)))
    #end for

    clear_all_sims()
#end def test_get_output_files
