
import testing
from testing import value_eq,object_eq
from testing import failed,FailedTest
from testing import divert_nexus_log,restore_nexus_log


def test_import():
    from project_manager import ProjectManager
#end def test_import



def test_init():
    from generic import obj
    from nexus_base import nexus_core
    from project_manager import ProjectManager

    pm = ProjectManager()

    modes = nexus_core.modes
    assert(pm.persistent_modes==set([modes.submit,modes.all]))
    assert(pm.simulations==obj())
    assert(pm.cascades==obj())
    assert(pm.progressing_cascades==obj())
#end def test_init



def test_add_simulations():
    from simulation import Simulation
    from project_manager import ProjectManager

    from test_simulation import get_test_workflow
    sims = get_test_workflow(1)

    pm = ProjectManager()
    pm.add_simulations(sims.list())

    assert(list(pm.cascades.keys())==[sims.s1.simid])
    assert(list(pm.progressing_cascades.keys())==[sims.s1.simid])
    assert(id(pm.cascades[sims.s1.simid])==id(sims.s1))
    assert(id(pm.progressing_cascades[sims.s1.simid])==id(sims.s1))

    assert(len(pm.simulations)==3)
    n = 0
    for s in sims:
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
    from simulation import Simulation
    from project_manager import ProjectManager

    from test_simulation import get_test_workflow,n_test_workflows

    sims = []
    for n in range(n_test_workflows):
        sims.extend(get_test_workflow(n).list())
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
    from generic import NexusError
    from simulation import Simulation
    from project_manager import ProjectManager

    from test_simulation import get_test_workflow,n_test_workflows

    sims = []
    for n in range(n_test_workflows):
        sims.extend(get_test_workflow(n).list())
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



def test_resolve_file_collisions():
    from generic import NexusError
    from simulation import Simulation
    from project_manager import ProjectManager

    from test_simulation import get_test_workflow,n_test_workflows

    divert_nexus_log()

    sims = []
    for n in range(n_test_workflows):
        sims.extend(get_test_workflow(n).list())
    #end for

    pm = ProjectManager()
    pm.add_simulations(sims)

    pm.resolve_file_collisions()

    s1 = sims[len(sims)//2]
    s2 = sims[len(sims)//2+1]

    s2.locdir  = s1.locdir
    s2.outfile = s1.outfile

    try:
        pm.resolve_file_collisions()
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    restore_nexus_log()

    Simulation.clear_all_sims()
#end def test_resolve_file_collisions



def test_propagate_blockages():
    from generic import NexusError
    from simulation import Simulation
    from project_manager import ProjectManager

    from test_simulation import get_test_workflow,n_test_workflows

    sims = []
    for n in range(n_test_workflows):
        sims.extend(get_test_workflow(n).list())
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
    from simulation import Simulation
    from project_manager import ProjectManager

    from test_simulation import get_test_workflow

    sims = get_test_workflow(1)

    pm = ProjectManager()
    pm.add_simulations(sims.list())

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
    from simulation import Simulation
    from project_manager import ProjectManager

    from test_simulation import get_test_workflow

    divert_nexus_log()

    sims = get_test_workflow(1)

    pm = ProjectManager()
    pm.add_simulations(sims.list())

    idc = id(pm.cascades)
    idp = id(pm.progressing_cascades)

    pm.check_dependencies()

    restore_nexus_log()

    Simulation.clear_all_sims()
#end def test_check_dependencies
