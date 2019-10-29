
import testing
from testing import failed,FailedTest


def test_import():
    from bundle import bundle
    from bundle import SimulationBundle
#end def test_import



def test_bundle():
    from generic import NexusError
    from bundle import bundle
    from bundle import SimulationBundle

    from test_simulation import get_test_workflow

    # empty init
    try:
        bundle()
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    # init with wrong types
    try:
        bundle([1,2,3])
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    sims = get_test_workflow(7)

    level1 = [sims.s0, sims.s1, sims.s2]
    level2 = [sims.s3, sims.s5, sims.s7, sims.s10]
    level3 = [sims.s4, sims.s6, sims.s8, sims.s11]
    level4 = [sims.s9, sims.s12]

    m = sims.s0.job.get_machine()
    m.batch_capable = True

    # attempt to bundle sims that depend on each other
    try:
        bundle(level2+[sims.s8],serial=True)
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try
    try:
        bundle(level2+[sims.s9],serial=True)
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    print '\n end test'
#end def test_bundle
