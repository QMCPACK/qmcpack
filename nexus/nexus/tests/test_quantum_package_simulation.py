import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.QUANTUM_PACKAGE_SIMULATION)

from ..generic import generic_settings
generic_settings.raise_error = True

from pathlib import Path
from . import isolate_nexus_core

from ..testing import failed,FailedTest
from ..testing import object_eq


def clear_all_sims():
    from ..quantum_package import QuantumPackage
    from nexus.simulation import Simulation
    QuantumPackage.qprc = None

    Simulation.clear_all_sims()
#end def clear_all_sims


def get_quantum_package_sim(**kwargs):
    from ..physical_system import generate_physical_system
    from ..machines import job
    from ..quantum_package import QuantumPackage,generate_quantum_package

    QuantumPackage.qprc = ''

    system = generate_physical_system(
        elem_pos = '''
            O  0.000000  0.000000  0.000000 
            H  0.000000  0.757160  0.586260
            H  0.000000  0.757160 -0.586260
            ''',
        )

    sim = generate_quantum_package(
        job      = job(machine='ws1',cores=1),
        system   = system,
        prefix   = 'h2o',
        run_type = 'scf',
        **kwargs
        )

    assert(isinstance(sim,QuantumPackage))
    
    return sim
#end def get_quantum_package_sim



def test_minimal_init():
    sim = get_quantum_package_sim()

    clear_all_sims()
#end def test_minimal_init



def test_check_result():
    sim = get_quantum_package_sim()
    
    assert(not sim.check_result('unknown',None))
    assert(not sim.check_result('orbitals',None))

    sim.input.run_control.save_for_qmcpack = True

    assert(sim.check_result('orbitals',None))

    clear_all_sims()
#end def test_check_result



def test_get_result():
    from ..developer import NexusError, obj

    sim = get_quantum_package_sim()
    
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

    sim.input.run_control.save_for_qmcpack = True

    result = sim.get_result('orbitals',None)

    result_ref = obj(
        outfile = './runs/qp_savewf.out',
        )

    assert(object_eq(result,result_ref))

    clear_all_sims()
#end def test_get_result



def test_incorporate_result():
    from ..developer import NexusError
    from ..machines import job
    from ..simulation import Simulation
    from ..gamess import generate_gamess,Gamess

    Gamess.ericfmt = ''

    other = Simulation()

    gms = generate_gamess(
        job = job(machine='ws1',cores=1),
        )

    Gamess.ericfmt = None

    sim = get_quantum_package_sim()
    
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
    
    try:
        sim.get_result('orbitals',other)
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    clear_all_sims()
#end def test_incorporate_result


@isolate_nexus_core
def test_check_sim_status(tmp_path):
    import os
    from ..nexus_base import nexus_core

    nexus_core.runs = ''
    nexus_core.local_directory  = tmp_path
    nexus_core.remote_directory = tmp_path
    nexus_core.file_locations = nexus_core.file_locations + [tmp_path]

    sim = get_quantum_package_sim()

    assert(Path(sim.locdir).resolve()==tmp_path)

    assert(not sim.finished)
    assert(not sim.failed)

    try:
        sim.check_sim_status()
        raise FailedTest
    except IOError:
        None
    except Exception as e:
        failed(str(e))
    #end try

    sim.create_directories()
    outfile = Path(sim.locdir).resolve() / sim.outfile
    outfile_text = '* SCF energy'
    outfile.write_text(outfile_text)
    assert(outfile_text in outfile.read_text())
    sim.job.finished = True

    sim.check_sim_status()

    assert(sim.finished)
    assert(not sim.failed)

    clear_all_sims()
#end def test_check_sim_status
