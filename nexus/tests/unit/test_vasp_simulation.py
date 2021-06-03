
import testing
from testing import divert_nexus,restore_nexus,clear_all_sims
from testing import failed,FailedTest
from testing import value_eq,object_eq,text_eq,check_object_eq

from test_vasp_input import c_potcar_text,get_files


pseudo_inputs = dict(
    pseudo_dir = 'pseudopotentials',
    pseudo_files_create = [('C.POTCAR',c_potcar_text)],
    )

def setup_vasp_sim(path,identifier='vasp',files=False):
    import shutil
    from nexus_base import nexus_core
    from machines import job
    from physical_system import generate_physical_system
    from vasp import generate_vasp,Vasp
    
    copy_files = files
    del files

    nexus_core.runs = ''

    files = get_files()

    dia16 = generate_physical_system(
        structure = files['d16bulk.POSCAR'],
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
            shutil.copy2(files[vfile],path)
        #end for
    #end if
        
    return sim
#end def setup_vasp_sim



def test_import():
    from vasp import Vasp,generate_vasp
#end def test_import



def test_minimal_init():
    from machines import job
    from vasp import Vasp,generate_vasp

    sim = generate_vasp(
        job = job(machine='ws1',cores=1),
        )

    assert(isinstance(sim,Vasp))

    clear_all_sims()
#end def test_minimal_init



def test_check_result():
    tpath = testing.setup_unit_test_output_directory('vasp_simulation','test_check_result',**pseudo_inputs)

    sim = setup_vasp_sim(tpath)

    assert(not sim.check_result('unknown',None))

    assert(sim.check_result('structure',None))

    clear_all_sims()
    restore_nexus()
#end def test_check_result



def test_get_result():
    import os
    import shutil
    from numpy import array
    from generic import obj,NexusError

    tpath = testing.setup_unit_test_output_directory('vasp_simulation','test_get_result',**pseudo_inputs)

    sim = setup_vasp_sim(tpath,identifier='diamond',files=True)

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

    pcfile = os.path.join(tpath,'diamond_POSCAR')
    ccfile = os.path.join(tpath,sim.identifier+'.CONTCAR')

    assert(sim.locdir.strip('/')==tpath.strip('/'))

    shutil.copy2(pcfile,ccfile)

    assert(os.path.exists(ccfile))

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
    restore_nexus()
#end def test_get_result



def test_incorporate_result():
    import os
    import shutil
    from numpy import array
    from generic import obj

    tpath = testing.setup_unit_test_output_directory('vasp_simulation','test_incorporate_result',**pseudo_inputs)

    sim = setup_vasp_sim(tpath,identifier='diamond',files=True)

    pcfile = os.path.join(tpath,'diamond_POSCAR')
    ccfile = os.path.join(tpath,sim.identifier+'.CONTCAR')

    assert(sim.locdir.strip('/')==tpath.strip('/'))

    shutil.copy2(pcfile,ccfile)

    assert(os.path.exists(ccfile))

    result = sim.get_result('structure',None)

    sim2 = setup_vasp_sim(tpath)

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
    restore_nexus()
#end def test_incorporate_result



def test_check_sim_status():
    import os

    tpath = testing.setup_unit_test_output_directory('vasp_simulation','test_check_sim_status',**pseudo_inputs)

    sim = setup_vasp_sim(tpath)

    assert(sim.locdir.strip('/')==tpath.strip('/'))

    assert(not sim.finished)
    assert(not sim.input.performing_neb())

    sim.check_sim_status()

    assert(not sim.finished)

    outcar_file = os.path.join(tpath,'OUTCAR')
    outcar_text = 'General timing and accounting'
    outcar = open(outcar_file,'w')
    outcar.write(outcar_text)
    outcar.close()
    assert(os.path.exists(outcar_file))
    assert(outcar_text in open(outcar_file,'r').read())

    sim.check_sim_status()

    assert(sim.finished)
    assert(not sim.failed)

    clear_all_sims()
    restore_nexus()
#end def test_check_sim_status



def test_get_output_files():
    import os

    tpath = testing.setup_unit_test_output_directory('vasp_simulation','test_get_output_files',**pseudo_inputs)

    sim = setup_vasp_sim(tpath)

    vfiles = 'INCAR KPOINTS POSCAR CONTCAR OUTCAR'.split()
    for vfile in vfiles:
        f = open(os.path.join(tpath,vfile),'w')
        f.write('')
        f.close()
    #end for

    files = sim.get_output_files()

    assert(value_eq(files,vfiles))

    for vfile in vfiles:
        assert(os.path.exists(os.path.join(tpath,sim.identifier+'.'+vfile)))
    #end for

    clear_all_sims()
    restore_nexus()
#end def test_get_output_files
