import testing
from testing import divert_nexus,restore_nexus,clear_all_sims
from testing import failed,FailedTest
from testing import value_eq,object_eq,text_eq


def get_pyscf_sim(**kwargs):
    from machines import job
    from pyscf_sim import Pyscf,generate_pyscf

    sim = generate_pyscf(
        job = job(machine='ws1',cores=1),
        **kwargs
        )

    assert(isinstance(sim,Pyscf))

    return sim
#end def get_pyscf_sim



def test_import():
    from pyscf_sim import Pyscf,generate_pyscf
#end def test_import



def test_minimal_init():
    from machines import job
    from pyscf_sim import Pyscf,generate_pyscf

    sim = generate_pyscf(
        job = job(machine='ws1',cores=1),
        )

    assert(isinstance(sim,Pyscf))

    clear_all_sims()
#end def test_minimal_init



def test_check_result():
    sim = get_pyscf_sim()
    
    assert(not sim.check_result('unknown',None))
    assert(not sim.check_result('orbitals',None))
    assert(not sim.check_result('wavefunction',None))

    sim.input.prefix   = 'scf'
    sim.input.save_qmc = True

    assert(sim.check_result('orbitals',None))
    assert(not sim.check_result('wavefunction',None))

    sim.input.prefix   = None
    sim.input.save_qmc = None
    sim.input.checkpoint = True

    assert(not sim.check_result('orbitals',None))
    assert(sim.check_result('wavefunction',None))

    clear_all_sims()
#end def test_check_result



def test_get_result():
    import os
    from generic import NexusError,obj
    from nexus_base import nexus_core

    tpath = testing.setup_unit_test_output_directory('pyscf_simulation','test_get_result',divert=True)

    nexus_core.runs = ''

    template_file = 'scf_template.py'
    template_text = 'template $chkfile'
    template_filepath = os.path.join(tpath,template_file)
    f = open(template_filepath,'w')
    f.write(template_text)
    f.close()

    sim = get_pyscf_sim(
        prefix     = 'scf',
        checkpoint = 'scf.chk',
        template   = template_filepath,
        )
    
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

    result = sim.get_result('orbitals',None)

    assert(result.h5_file.replace(tpath,'').lstrip('/')=='scf.h5')

    result = sim.get_result('wavefunction',None)

    assert(result.chkfile.replace(tpath,'').lstrip('/')=='scf.chk')

    clear_all_sims()
    restore_nexus()
#end def test_get_result



def test_check_sim_status():
    sim = get_pyscf_sim()

    assert(not sim.failed)
    assert(not sim.finished)

    sim.check_sim_status()

    assert(not sim.failed)
    assert(not sim.finished)

    sim.job.finished = True

    sim.check_sim_status()

    assert(not sim.failed)
    assert(sim.finished)

    clear_all_sims()
#end def test_check_sim_status
