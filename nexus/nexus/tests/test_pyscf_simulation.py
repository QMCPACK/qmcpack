import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.PYSCF_SIMULATION)

from ..generic import generic_settings, obj
generic_settings.raise_error = True

from pathlib import Path
from . import isolate_nexus_core, TEST_DIR
from nexus.nexus_base import nexus_core
from nexus.physical_system import generate_physical_system
from nexus.structure import generate_trimer_structure
from ..testing import clear_all_sims
from ..testing import failed,FailedTest

TEST_FILES = {
    "scf_template.py": TEST_DIR / "test_pyscf_simulation_files/scf_template.py",
}


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


@isolate_nexus_core
def test_get_result(tmp_path):
    from ..developer import NexusError
    from ..nexus_base import nexus_core

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

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



@isolate_nexus_core
def test_check_sim_status(tmp_path):

    nexus_core.runs = ''
    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]

    # Water
    structure = generate_trimer_structure(
        trimer     = ["O", "H", "H"],
        units      = "A",
        separation = [1.0, 1.0],
        angle      = 104.5,
        )

    system = generate_physical_system(
        structure = structure,
        )

    sim = get_pyscf_sim(
        identifier = 'scf',
        path       = 'scf',
        template   = TEST_FILES["scf_template.py"],
        mole       = obj(
            verbose  = 5,
            basis    = 'ccpvtz',
            symmetry = True,
            ),
        save_qmc   = True,
        system     = system,
        )

    sim_dir = Path(sim.locdir).resolve()
    assert(sim_dir == (tmp_path / 'scf').resolve())
    sim_dir.mkdir()

    assert(not sim.failed)
    assert(not sim.finished)

    try:
        sim.check_sim_status()
    except IOError:
        None
    except Exception as e:
        raise e

    assert(not sim.failed)
    assert(not sim.finished)

    out_path = sim_dir / sim.outfile
    out_path.touch()
    assert(out_path.exists())

    err_path = sim_dir / sim.errfile
    err_path.touch()
    assert(err_path.exists())

    sim.check_sim_status()

    assert(not sim.failed)
    assert(not sim.finished)

    sim.job.finished = True

    sim.check_sim_status()

    assert(not sim.failed)
    assert(sim.finished)

    err_text = """
Traceback (most recent call last):
  File "/dummy/path", line 1, in <module>
    from pyscf import scf
ModuleNotFoundError: No module named 'pyscf'
"""
    err_path.write_text(err_text)
    assert(err_text in err_path.read_text())

    sim.check_sim_status()

    assert(sim.finished)
    assert(sim.failed)

    clear_all_sims()
#end def test_check_sim_status
