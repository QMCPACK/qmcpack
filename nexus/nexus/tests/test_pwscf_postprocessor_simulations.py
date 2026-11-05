import pytest

from ..machines import job
from ..pwscf_postprocessors import (
    PP,
    Bands,
    Cppp,
    Dos,
    Projwfc,
    Pwexport,
    generate_bands,
    generate_cppp,
    generate_dos,
    generate_pp,
    generate_projwfc,
    generate_pwexport,
)
from ..testing import FailedTest, clear_all_sims, failed
from . import NexusTestOrder

pytestmark = pytest.mark.order(NexusTestOrder.PWSCF_POSTPROCESSOR_SIMULATIONS)

def get_class_generators():
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



def test_minimal_init():
    for cls,gen in get_class_generators():
        sim = gen(job=job(machine='ws1',cores=1))
        assert(isinstance(sim,cls))
    #end for

    clear_all_sims()
#end def test_minimal_init



def test_check_result():
    for cls,gen in get_class_generators():
        sim = gen(job=job(machine='ws1',cores=1))
        assert(isinstance(sim,cls))
        assert(not sim.check_result('anything',None))
    #end for

    clear_all_sims()
#end def test_check_result



def test_get_result():
    for cls,gen in get_class_generators():
        sim = gen(job=job(machine='ws1',cores=1))
        assert(isinstance(sim,cls))
        try:
            sim.get_result('anything',None)
            raise FailedTest
        except NotImplementedError:
            pass
        except FailedTest:
            failed()
        except Exception as e:
            failed(str(e))
        #end try
    #end for

    clear_all_sims()
#end def test_get_result



def test_incorporate_result():
    for cls,gen in get_class_generators():
        sim = gen(job=job(machine='ws1',cores=1))
        assert(isinstance(sim,cls))
        try:
            sim.incorporate_result('anything',None,None)
            raise FailedTest
        except NotImplementedError:
            pass
        except FailedTest:
            failed()
        except Exception as e:
            failed(str(e))
        #end try
    #end for

    clear_all_sims()
#end def test_incorporate_result



def test_check_sim_status():
    for cls,gen in get_class_generators():
        sim = gen(job=job(machine='ws1',cores=1))
        assert(isinstance(sim,cls))
        sim.check_sim_status()
        assert(sim.finished)
        assert(not sim.failed)
    #end for

    clear_all_sims()
#end def test_check_sim_status
