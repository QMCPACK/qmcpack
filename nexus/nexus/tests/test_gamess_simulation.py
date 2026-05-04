import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.GAMESS_SIMULATION)

from ..generic import generic_settings
generic_settings.raise_error = True

from pathlib import Path
from . import isolate_nexus_core
from .. import testing
from ..testing import failed,FailedTest
from ..testing import object_eq


def clear_all_sims():
    from ..gamess import Gamess

    Gamess.ericfmt = None

    testing.clear_all_sims()
#end def clear_all_sims


def get_gamess_sim(type='rhf'):
    from ..machines import job
    from ..gamess import Gamess,generate_gamess,GamessInput
    from .test_gamess_input import TEST_FILES

    sim = None

    Gamess.ericfmt = ''

    if type=='rhf':
        gi_input = GamessInput(TEST_FILES['rhf.inp'])

        sim = generate_gamess(
            identifier = 'rhf',
            path       = 'rhf',
            job        = job(machine='ws1',cores=1),
            input      = gi_input,
            )
    else:
        failed()
    #end if

    assert(sim is not None)
    assert(isinstance(sim,Gamess))

    return sim
#end def get_gamess_sim



def test_minimal_init():
    from ..machines import job
    from ..gamess import Gamess,generate_gamess

    Gamess.ericfmt = ''

    sim = generate_gamess(
        job = job(machine='ws1',cores=1),
        )

    assert(isinstance(sim,Gamess))

    clear_all_sims()
#end def test_minimal_init



def test_check_result():

    sim = get_gamess_sim('rhf')
    
    assert(not sim.check_result('unknown',None))
    assert(sim.check_result('orbitals',None))

    clear_all_sims()
#end def test_check_result


@isolate_nexus_core
def test_get_result(tmp_path):
    from ..developer import obj, NexusError
    from ..nexus_base import nexus_core

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = tmp_path
    nexus_core.file_locations = nexus_core.file_locations + [tmp_path]
    nexus_core.runs = ''

    sim = get_gamess_sim('rhf')

    assert(Path(sim.locdir).resolve()==tmp_path / 'rhf')

    sim.create_directories()
    sim.write_inputs()
    Path(sim.imresdir).resolve().mkdir(exist_ok=True)
    analyzer = sim.analyzer_type(sim)
    analyzer.save(Path(sim.imresdir).resolve() / sim.analyzer_image)

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

    result_ref = obj(
        location        = 'rhf/rhf.out',
        mos             = 0,
        norbitals       = 0,
        outfile         = 'rhf/rhf.out',
        scftyp          = 'rohf',
        vec             = None,
        )

    result.location = result.location.replace(str(tmp_path),'').lstrip('/')
    result.outfile  = result.outfile.replace(str(tmp_path),'').lstrip('/')

    assert(object_eq(result,result_ref))

    clear_all_sims()
#end def test_get_result



def test_incorporate_result():

    from ..developer import NexusError, obj

    sim = get_gamess_sim('rhf')

    try:
        sim.incorporate_result('unknown',None,None)
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    result = obj(
        vec       = 'vec text',
        norbitals = 10,
        )

    input = sim.input

    assert('vec' not in input)
    assert('norb' not in input.guess)
    assert('prtmo' not in input.guess)

    sim.incorporate_result('orbitals',result,None)

    assert(input.vec.text=='vec text')
    assert(input.guess.guess=='moread')
    assert(input.guess.norb==10)
    assert(input.guess.prtmo==True)

    clear_all_sims()
#end def test_incorporate_result


@isolate_nexus_core
def test_check_sim_status(tmp_path):
    from ..nexus_base import nexus_core

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = tmp_path
    nexus_core.file_locations = nexus_core.file_locations + [tmp_path]
    nexus_core.runs = ''

    sim = get_gamess_sim('rhf')

    assert(Path(sim.locdir).resolve()==tmp_path / 'rhf')

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
    outfile_text = 'EXECUTION OF GAMESS TERMINATED NORMALLY'
    outfile.write_text(outfile_text)
    assert(outfile_text in outfile.read_text())

    sim.check_sim_status()

    assert(sim.finished)
    assert(not sim.failed)

    clear_all_sims()
#end def test_check_sim_status
