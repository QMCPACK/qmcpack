
import testing
from testing import divert_nexus,restore_nexus,clear_all_sims
from testing import failed,FailedTest
from testing import value_eq,object_eq,text_eq



pseudo_inputs = dict(
    pseudo_dir = 'pseudopotentials',
    pseudo_files_create = ['C.BFD.upf'],
    )


def get_system():
    from physical_system import generate_physical_system

    system = generate_physical_system(
        units  = 'A',
        axes   = [[ 1.785,  1.785,  0.   ],
                  [ 0.   ,  1.785,  1.785],
                  [ 1.785,  0.   ,  1.785]],
        elem   = ['C','C'],
        pos    = [[ 0.    ,  0.    ,  0.    ],
                  [ 0.8925,  0.8925,  0.8925]],
        tiling = (1,1,1),
        kgrid  = (1,1,1),
        kshift = (0,0,0),
        C      = 4
        )

    return system
#end def get_system


def get_pwscf_sim(type='scf'):
    from nexus_base import nexus_core
    from machines import job
    from pwscf import Pwscf,generate_pwscf

    nexus_core.runs = ''

    sim = None

    if type=='scf':
        sim = generate_pwscf(
            identifier   = 'scf',
            path         = 'scf',
            job          = job(machine='ws1',cores=1),
            input_type   = 'generic',
            calculation  = 'scf',
            input_dft    = 'lda', 
            ecutwfc      = 200,   
            conv_thr     = 1e-8, 
            nosym        = True,
            wf_collect   = True,
            system       = get_system(),
            pseudos      = ['C.BFD.upf'], 
            nogamma      = True,
            )
    else:
        failed()
    #end if

    assert(sim is not None)
    assert(isinstance(sim,Pwscf))

    return sim
#end def get_pwscf_sim



def test_import():
    from pwscf import Pwscf,generate_pwscf
#end def test_import



def test_minimal_init():
    from machines import job
    from pwscf import Pwscf,generate_pwscf

    sim = generate_pwscf(
        job    = job(machine='ws1',cores=1),
        system = get_system(),
        )

    assert(isinstance(sim,Pwscf))

    clear_all_sims()
#end def test_minimal_init



def test_check_result():
    tpath = testing.setup_unit_test_output_directory('pwscf_simulation','test_check_result',**pseudo_inputs)

    sim = get_pwscf_sim('scf')
    
    assert(not sim.check_result('unknown',None))
    assert(sim.check_result('charge_density',None))
    assert(sim.check_result('restart',None))
    assert(sim.check_result('orbitals',None))
    assert(not sim.check_result('structure',None))

    clear_all_sims()
    restore_nexus()
#end def test_check_result



def test_get_result():
    from generic import NexusError

    tpath = testing.setup_unit_test_output_directory('pwscf_simulation','test_get_result',**pseudo_inputs)

    sim = get_pwscf_sim('scf')

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

    result  = sim.get_result('charge_density',None)
    result2 = sim.get_result('restart',None)

    assert(object_eq(result,result2))

    result_ref = dict(
        locdir        = 'scf',
        outdir        = 'scf/pwscf_output',
        location      = 'scf/pwscf_output/pwscf.save/charge-density.dat',
        spin_location = 'scf/pwscf_output/pwscf.save/spin-polarization.dat',
        )

    assert(set(result.keys())==set(result_ref.keys()))
    for k,path in result.items():
        path = path.replace(tpath,'').lstrip('/')
        assert(path==result_ref[k])
    #end for

    result = sim.get_result('orbitals',None)

    result_ref = dict(
        location = 'scf/pwscf_output/pwscf.wfc1',
        )

    assert(set(result.keys())==set(result_ref.keys()))
    for k,path in result.items():
        path = path.replace(tpath,'').lstrip('/')
        assert(path==result_ref[k])
    #end for

    clear_all_sims()
    restore_nexus()
#end def test_get_result



def test_incorporate_result():
    from generic import obj
    tpath = testing.setup_unit_test_output_directory('pwscf_simulation','test_incorporate_result',**pseudo_inputs)

    sim = get_pwscf_sim('scf')

    sim_start = sim.to_obj().copy()

    assert(object_eq(sim.to_obj(),sim_start))

    # charge density
    result = obj(
        locdir = sim.locdir,
        )

    sim.incorporate_result('charge_density',result,None)

    assert(object_eq(sim.to_obj(),sim_start))

    # restart
    sim.incorporate_result('restart',result,None)

    c = sim.input.control
    assert(c.restart_mode=='restart')
    del c.restart_mode
    assert(object_eq(sim.to_obj(),sim_start))
    
    # structure
    altered_structure = sim.system.structure.copy()
    altered_structure.pos += 0.1

    result = obj(
        structure = altered_structure,
        )

    sim.incorporate_result('structure',result,None)

    sim_ref = sim_start.copy()
    pos_ref = sim_ref.system.structure.delete('pos')+0.1
    sim_ref.input.atomic_positions.delete('positions')

    pos = sim.system.structure.delete('pos')
    apos = sim.input.atomic_positions.delete('positions')

    assert(value_eq(pos,pos_ref))
    assert(value_eq(apos,pos_ref))
    assert(object_eq(sim.to_obj(),sim_ref))

    clear_all_sims()
    restore_nexus()
#end def test_incorporate_result



def test_check_sim_status():
    import os
    tpath = testing.setup_unit_test_output_directory('pwscf_simulation','test_check_sim_status',**pseudo_inputs)

    sim = get_pwscf_sim('scf')

    assert(sim.locdir.rstrip('/')==os.path.join(tpath,'scf').rstrip('/'))
    if not os.path.exists(sim.locdir):
        os.makedirs(sim.locdir)
    #end if

    assert(not sim.finished)

    try:
        sim.check_sim_status()
    except IOError:
        None
    except Exception as e:
        failed(str(e))
    #end try

    assert(not sim.finished)

    out_path = os.path.join(tpath,'scf',sim.outfile)

    out_text = ''
    outfile = open(out_path,'w')
    outfile.write(out_text)
    outfile.close()
    assert(os.path.exists(out_path))

    sim.check_sim_status()

    assert(not sim.finished)

    out_text = 'JOB DONE'
    outfile = open(out_path,'w')
    outfile.write(out_text)
    outfile.close()
    assert(out_text in open(out_path,'r').read())

    sim.check_sim_status()

    assert(sim.finished)
    assert(not sim.failed)

    clear_all_sims()
    restore_nexus()
#end def test_check_sim_status

