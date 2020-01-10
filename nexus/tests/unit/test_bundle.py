
import testing
from testing import failed,FailedTest
from testing import value_eq,object_eq


def test_import():
    from bundle import bundle
    from bundle import SimulationBundle
#end def test_import



def test_bundle():
    from generic import NexusError
    from machines import job,get_machine
    from bundle import bundle
    from bundle import SimulationBundle

    from test_simulation_module import get_test_workflow

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

    def get_workflow():
        machine = get_machine('theta')
        machine.account = 'ABC123'
        test_job = job(machine='theta',nodes=1,app_command='test.x')
        sims = get_test_workflow(7,job=test_job)

        levels = []
        levels.append([sims.s0, sims.s1, sims.s2])
        levels.append([sims.s3, sims.s5, sims.s7, sims.s10])
        levels.append([sims.s4, sims.s6, sims.s8, sims.s11])
        levels.append([sims.s9, sims.s12])

        return sims,levels
    #end def get_workflow

    sims,levels = get_workflow()

    # attempt to bundle sims that depend on each other
    try:
        bundle(levels[1]+[sims.s8])
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try
    try:
        bundle(levels[1]+[sims.s9])
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    sims,levels = get_workflow()

    for n,level in enumerate(levels):
        for s in level:
            assert(not s.bundled)
            assert(s.bundler is None)
            assert(not s.skip_submit)
        #end for

        b = bundle(level,identifier='bundle_'+str(n))

        for s in level:
            assert(s.bundled)
            assert(id(s.bundler)==id(b))
            assert(s.skip_submit)
        #end for

        assert(isinstance(b,SimulationBundle))
        assert(isinstance(b.sims,list))
        assert(b.system is None)
        assert(b.infile is None)
        assert(not b.allow_create_directories)
        assert(not b.allow_get_dependencies  )
        assert(not b.allow_write_inputs      )
        assert(not b.allow_send_files        )
        assert(not b.allow_submit            )
        assert(not b.allow_get_output        )
        assert(not b.allow_analyze           )

        for sb,s in zip(b.sims,level):
            assert(id(sb)==id(s))
        #end for

        bj = b.job
        j  = level[0].job
        ns = len(level)
        assert(bj.nodes==ns*j.nodes)
        assert(bj.cores==ns*j.cores)
        assert(object_eq(bj.get_time(),j.get_time()))

    #end for

    testing.clear_all_sims()
#end def test_bundle
