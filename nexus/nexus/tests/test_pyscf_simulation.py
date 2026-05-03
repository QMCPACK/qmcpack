import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.PYSCF_SIMULATION)

from ..generic import generic_settings
generic_settings.raise_error = True

from . import isolate_nexus_core
from ..testing import clear_all_sims
from ..testing import failed,FailedTest


def get_pyscf_sim(**kwargs):
    from ..machines import job
    from ..pyscf_sim import Pyscf,generate_pyscf

    sim = generate_pyscf(
        job = job(machine='ws1',cores=1),
        **kwargs
        )

    assert(isinstance(sim,Pyscf))

    return sim
#end def get_pyscf_sim



def test_minimal_init():
    from ..machines import job
    from ..pyscf_sim import Pyscf,generate_pyscf

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


@isolate_nexus_core(needs_tmp_path=True)
def test_get_result(tmp_path):
    from ..developer import NexusError
    from ..nexus_base import nexus_core

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = tmp_path
    nexus_core.file_locations = nexus_core.file_locations + [tmp_path]

    nexus_core.runs = ''

    template_file = 'scf_template.py'
    template_text = 'template $chkfile'
    template_filepath = tmp_path / template_file
    template_filepath.write_text(template_text)

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

    assert(result.h5_file.replace(str(tmp_path),'').lstrip('/')=='scf.h5')

    result = sim.get_result('wavefunction',None)

    assert(result.chkfile.replace(str(tmp_path),'').lstrip('/')=='scf.chk')

    clear_all_sims()
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
