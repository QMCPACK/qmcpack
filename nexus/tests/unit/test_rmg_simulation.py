
import testing
from testing import divert_nexus,restore_nexus,clear_all_sims
from testing import failed,FailedTest
from testing import value_eq,object_eq,text_eq



def test_import():
    from rmg import Rmg,generate_rmg
#end def test_import



def test_minimal_init():
    from machines import job
    from rmg import Rmg,generate_rmg

    sim = generate_rmg(
        job    = job(machine='ws1',cores=1),
        )

    assert(isinstance(sim,Rmg))

    clear_all_sims()
#end def test_minimal_init

