import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.PROJECT_MANAGER)

from . import isolate_nexus_core
from ..testing import value_eq
from ..testing import failed,FailedTest


def test_init():
    from ..developer import obj
    from ..project_manager import ProjectManager

    pm = ProjectManager()

    def check(v):
        assert isinstance(v,obj)
        assert len(v)==0
    check(pm.simulations)
    check(pm.cascades)
    check(pm.progressing_cascades)
#end def test_init



def test_add_simulations():
    from ..simulation import Simulation
    from ..project_manager import ProjectManager

    from .test_simulation_module import get_test_workflow
    sims = get_test_workflow(1)

    pm = ProjectManager()
    pm.add_simulations(list(sims.values()))

    assert(list(pm.cascades.keys())==[sims.s1.simid])
    assert(list(pm.progressing_cascades.keys())==[sims.s1.simid])
    assert(id(pm.cascades[sims.s1.simid])==id(sims.s1))
    assert(id(pm.progressing_cascades[sims.s1.simid])==id(sims.s1))

    assert(len(pm.simulations)==3)
    n = 0
    for s in sims.values():
        assert(s.simid in pm.simulations)
        assert(id(pm.simulations[s.simid])==id(s))
        n+=1
    #end for
    assert(n==3)

    pm = ProjectManager()
    pm.add_simulations()
    assert(len(pm.simulations)==len(Simulation.all_sims))
    n = 0
    for s in Simulation.all_sims:
        assert(s.simid in pm.simulations)
        assert(id(pm.simulations[s.simid])==id(s))
        n+=1
    #end for
    assert(n>=3)

    Simulation.clear_all_sims()
#end def test_add_simulations



def test_traverse_cascades():
    from ..simulation import Simulation
    from ..project_manager import ProjectManager

    from .test_simulation_module import get_test_workflow,n_test_workflows

    sims = []
    for n in range(n_test_workflows):
        wf = get_test_workflow(n)
        sims.extend(list(wf.values()))
    #end for

    pm = ProjectManager()
    pm.add_simulations(sims)

    pm.traverse_cascades()
    
    def count_visits(sim,visit_counts):
        i = sim.simid
        if i not in visit_counts:
            visit_counts[i] = 1
        else:
            visit_counts[i] += 1
        #end if
    #end def count_visits

    counts = dict()

    pm.traverse_cascades(count_visits,counts)

    assert(len(counts)==len(sims))
    for s in sims:
        assert(s.simid in counts)
        assert(counts[s.simid]==1)
    #end for

    Simulation.clear_all_sims()
#end def test_traverse_cascades



def test_screen_fake_sims():
    from ..developer import NexusError
    from ..simulation import Simulation
    from ..project_manager import ProjectManager

    from .test_simulation_module import get_test_workflow,n_test_workflows

    sims = []
    for n in range(n_test_workflows):
        wf = get_test_workflow(n)
        sims.extend(list(wf.values()))
    #end for

    pm = ProjectManager()
    pm.add_simulations(sims)

    pm.screen_fake_sims()

    s = sims[len(sims)//2]
    s.fake_sim = True

    try:
        pm.screen_fake_sims()
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    Simulation.clear_all_sims()
#end def test_screen_fake_sims


@isolate_nexus_core
def test_resolve_file_collisions():
    from ..developer import NexusError
    from ..simulation import Simulation
    from ..project_manager import ProjectManager

    from .test_simulation_module import get_test_workflow,n_test_workflows

    sims = []
    for n in range(n_test_workflows):
        wf = get_test_workflow(n)
        sims.extend(list(wf.values()))
    #end for

    pm = ProjectManager()
    pm.add_simulations(sims)

    pm.resolve_file_collisions()

    s1 = sims[len(sims)//2]
    s2 = sims[len(sims)//2+1]

    s2.locdir  = s1.locdir
    s2.outfile = s1.outfile

    with pytest.raises(
        FileExistsError,
        match="file collisions found in directory",
        ):
        pm.resolve_file_collisions()

    Simulation.clear_all_sims()
#end def test_resolve_file_collisions



def test_propagate_blockages():
    from ..simulation import Simulation
    from ..project_manager import ProjectManager

    from .test_simulation_module import get_test_workflow,n_test_workflows

    sims = []
    for n in range(n_test_workflows):
        wf = get_test_workflow(n)
        sims.extend(list(wf.values()))
    #end for

    for s in sims:
        assert(not s.block)
        assert(not s.block_subcascade)
    #end for

    pm = ProjectManager()
    pm.add_simulations(sims)

    pm.propagate_blockages()

    nb = 5
    dn = len(sims)//nb
    for b in range(nb):
        s = sims[nb*b]
        s.block = True
    #end for

    pm.propagate_blockages()

    def assert_blocked(sim):
        assert(sim.block)
        assert(sim.block_subcascade)
    #end def assert_blocked

    for b in range(nb):
        s = sims[nb*b]
        s.traverse_full_cascade(assert_blocked)
    #end for

    Simulation.clear_all_sims()
#end def test_propagate_blockages



def test_load_cascades():
    from ..simulation import Simulation
    from ..project_manager import ProjectManager

    from .test_simulation_module import get_test_workflow

    sims = get_test_workflow(1)

    pm = ProjectManager()
    pm.add_simulations(list(sims.values()))

    idc = id(pm.cascades)
    idp = id(pm.progressing_cascades)

    pm.load_cascades()

    assert(list(pm.cascades.keys())==[sims.s1.simid])
    assert(list(pm.progressing_cascades.keys())==[sims.s1.simid])
    assert(id(pm.cascades[sims.s1.simid])==id(sims.s1))
    assert(id(pm.progressing_cascades[sims.s1.simid])==id(sims.s1))

    assert(id(pm.cascades)!=idc)
    assert(id(pm.progressing_cascades)!=idp)

    Simulation.clear_all_sims()
#end def test_load_cascades



def test_check_dependencies():
    from ..simulation import Simulation
    from ..project_manager import ProjectManager

    from .test_simulation_module import get_test_workflow

    sims = get_test_workflow(1)

    pm = ProjectManager()
    pm.add_simulations(list(sims.values()))

    idc = id(pm.cascades)
    idp = id(pm.progressing_cascades)

    pm.check_dependencies()

    Simulation.clear_all_sims()
#end def test_check_dependencies


@isolate_nexus_core
def test_write_simulation_status():
    from ..generic import generic_settings
    from ..nexus_base import nexus_core
    from ..simulation import Simulation
    from ..project_manager import ProjectManager

    from .test_simulation_module import get_test_workflow

    log = generic_settings.devlog

    sims = get_test_workflow(3)
    for id,sim in sims.items():
        sim.identifier = 'test_sim_'+id
    #end for

    pm = ProjectManager()
    pm.add_simulations(list(sims.values()))

    def status_log():
        log.reset()
        pm.write_simulation_status()
        s = log.contents()
        return '\n'.join(line.rstrip() for line in s.splitlines())
    #end def status_log

    status_ref = '''
  cascade status
    setup, sent_files, submitted, finished, got_output, analyzed, failed
    000000           ------    test_sim_s11  ./runs/
    000000           ------    test_sim_s21  ./runs/
    000000           ------    test_sim_s12  ./runs/
    000000           ------    test_sim_s22  ./runs/
    000000           ------    test_sim_s3  ./runs/
    000000           ------    test_sim_s4  ./runs/
    000000           ------    test_sim_s5  ./runs/
    setup, sent_files, submitted, finished, got_output, analyzed, failed
    '''
    assert(status_log().strip()==status_ref.strip())

    sim = sims.s11
    sim.setup      = True
    sim.sent_files = True
    sim.submitted  = True

    log.reset()
    pm.status_line(sim)
    assert(log.contents().strip()=='111000           ------    test_sim_s11  ./runs/')

    sim.finished  = True
    sim.got_output = True
    sim.analyzed   = True

    log.reset()
    pm.status_line(sim)
    assert(log.contents().strip()=='111111  SUCCESS  ------    test_sim_s11  ./runs/')

    sim.failed = True

    log.reset()
    pm.status_line(sim)
    assert(log.contents().strip()=='111111  FAILURE  ------    test_sim_s11  ./runs/')

    Simulation.clear_all_sims()
#end def test_write_simulation_status


def test_color_status_result(monkeypatch):
    from ..project_manager import color_status_result

    class TestLog:
        def __init__(self,is_tty):
            self.is_tty = is_tty
        #end def __init__

        def isatty(self):
            return self.is_tty
        #end def isatty
    #end class TestLog

    tty_log = TestLog(is_tty=True)
    file_log = TestLog(is_tty=False)

    monkeypatch.delenv('NO_COLOR',raising=False)
    assert(color_status_result('SUCCESS',tty_log)=='\033[42mSUCCESS\033[0m')
    assert(color_status_result('FAILURE',tty_log)=='\033[41mFAILURE\033[0m')
    assert(color_status_result('',tty_log)=='')
    assert(color_status_result('SUCCESS',file_log)=='SUCCESS')
    assert(color_status_result('FAILURE',file_log)=='FAILURE')

    monkeypatch.setenv('NO_COLOR','1')
    assert(color_status_result('SUCCESS',tty_log)=='SUCCESS')
    assert(color_status_result('FAILURE',tty_log)=='FAILURE')

#end def test_color_status_result


@isolate_nexus_core
def test_run_project(tmp_path):
    from ..generic import generic_settings
    from ..nexus_base import nexus_core
    from ..simulation import Simulation,input_template
    from ..project_manager import ProjectManager

    from .test_simulation_module import get_test_workflow,n_test_workflows

    # divert_nexus()
    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    nexus_core.sleep = 0.1

    flags = ['setup','sent_files','submitted','finished','got_output','analyzed']

    def finished(s):
        f = True
        for flag in flags:
            f &= s[flag]
        #end for
        f &= isinstance(s.process_id,int)
        return f
    #end def finished

    def empty(s):
        e = True
        for flag in flags:
            e &= not s[flag]
        #end for
        e &= s.process_id is None
        e &= s.job.system_id is None
        return e
    #end def empty

    sims = []
    for n in range(n_test_workflows):
        wf = get_test_workflow(n)
        sims.extend(list(wf.values()))
    #end for

    template = '''
name = "$name"
a    = $a
'''
    for s in sims:
        si = input_template(text=template)
        si.assign(name='input_name',a=1)
        s.input = si
    #end for

    for s in sims:
        assert(empty(s))
    #end for

    pm = ProjectManager()
    pm.machine = sims[0].job.get_machine()
    pm.add_simulations(sims)

    #pm.write_simulation_status()

    pm.run_project()

    #pm.write_simulation_status()
    #print(log.contents())

    for s in sims:
        assert(finished(s))
    #end for

    # restore_nexus()

    Simulation.clear_all_sims()
#end def test_run_project
