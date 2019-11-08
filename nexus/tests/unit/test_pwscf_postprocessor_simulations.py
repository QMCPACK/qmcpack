
import testing
from testing import divert_nexus,restore_nexus,clear_all_sims
from testing import failed,FailedTest
from testing import value_eq,object_eq,text_eq


def get_class_generators():
    from pwscf_postprocessors import PP,generate_pp
    from pwscf_postprocessors import Dos,generate_dos
    from pwscf_postprocessors import Bands,generate_bands
    from pwscf_postprocessors import Projwfc,generate_projwfc
    from pwscf_postprocessors import Cppp,generate_cppp
    from pwscf_postprocessors import Pwexport,generate_pwexport

    class_generators = [
        (PP,generate_pp),
        (Dos,generate_dos),
        (Bands,generate_bands),
        (Projwfc,generate_projwfc),
        (Cppp,generate_cppp),
        (Pwexport,generate_pwexport),
        ]

    return class_generators
#end def get_class_generators



def test_import():
    from pwscf_postprocessors import PP,generate_pp
    from pwscf_postprocessors import Dos,generate_dos
    from pwscf_postprocessors import Bands,generate_bands
    from pwscf_postprocessors import Projwfc,generate_projwfc
    from pwscf_postprocessors import Cppp,generate_cppp
    from pwscf_postprocessors import Pwexport,generate_pwexport
#end def test_import



def test_minimal_init():
    from machines import job

    for cls,gen in get_class_generators():
        sim = gen(job=job(machine='ws1',cores=1))
        assert(isinstance(sim,cls))
    #end for

    clear_all_sims()
#end def test_minimal_init



def test_check_result():
    from machines import job

    for cls,gen in get_class_generators():
        sim = gen(job=job(machine='ws1',cores=1))
        assert(isinstance(sim,cls))
        assert(not sim.check_result('anything',None))
    #end for

    clear_all_sims()
#end def test_check_result



def test_get_result():
    from generic import NexusError
    from machines import job

    for cls,gen in get_class_generators():
        sim = gen(job=job(machine='ws1',cores=1))
        assert(isinstance(sim,cls))
        try:
            sim.get_result('anything',None)
            raise FailedTest
        except NexusError:
            None
        except FailedTest:
            failed()
        except Exception as e:
            failed(str(e))
        #end try
    #end for

    clear_all_sims()
#end def test_get_result



def test_incorporate_result():
    from generic import NexusError
    from machines import job

    for cls,gen in get_class_generators():
        sim = gen(job=job(machine='ws1',cores=1))
        assert(isinstance(sim,cls))
        try:
            sim.incorporate_result('anything',None,None)
            raise FailedTest
        except NexusError:
            None
        except FailedTest:
            failed()
        except Exception as e:
            failed(str(e))
        #end try
    #end for

    clear_all_sims()
#end def test_incorporate_result



def test_check_sim_status():
    from machines import job

    for cls,gen in get_class_generators():
        sim = gen(job=job(machine='ws1',cores=1))
        assert(isinstance(sim,cls))
        sim.check_sim_status()
        assert(sim.finished)
        assert(not sim.failed)
    #end for

    clear_all_sims()
#end def test_check_sim_status
