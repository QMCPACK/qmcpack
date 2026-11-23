import shutil
from pathlib import Path

import pytest

from ..developer import obj
from ..machines import job
from ..nexus_base import nexus_config
from ..rmg import Rmg, generate_rmg
from ..testing import clear_all_sims
from . import NexusTestOrder, isolate_nexus_core

pytestmark = pytest.mark.order(NexusTestOrder.RMG_SIMULATION)

def test_minimal_init():
    sim = generate_rmg(
        job    = job(machine='ws1',cores=1),
        )

    assert(isinstance(sim,Rmg))

    clear_all_sims()
#end def test_minimal_init


def get_rmg_sim(path,calculation_mode,**kwargs):
    sim = generate_rmg(
        identifier        = calculation_mode.split()[0].lower(),
        path              = path,
        job               = job(machine='ws1',cores=1),
        calculation_mode  = calculation_mode,
        wavefunction_grid = (16,16,16),
        **kwargs
        )
    assert(isinstance(sim,Rmg))
    return sim
#end def get_rmg_sim


@isolate_nexus_core
def test_check_result(tmp_path):
    nexus_config.runs             = ''
    nexus_config.local_directory  = str(tmp_path)
    nexus_config.remote_directory = str(tmp_path)

    scf  = get_rmg_sim('scf','Quench Electrons')
    nscf = get_rmg_sim('nscf','NSCF',write_qmcpack_restart=True)

    assert(not scf.check_result('unknown',None))
    assert(scf.check_result('charge_density',None))
    assert(not scf.check_result('orbitals',None))

    assert(not nscf.check_result('charge_density',None))
    assert(nscf.check_result('orbitals',None))

    assert(not scf.sync_from_scf)
    assert(nscf.sync_from_scf)

    clear_all_sims()
#end def test_check_result


@isolate_nexus_core
def test_get_result(tmp_path):
    nexus_config.runs             = ''
    nexus_config.local_directory  = str(tmp_path)
    nexus_config.remote_directory = str(tmp_path)

    scf = get_rmg_sim('scf','Quench Electrons',write_qmcpack_restart=True)

    with pytest.raises(NotImplementedError,match='ability to get result unknown'):
        scf.get_result('unknown',None)
    #end with

    result = scf.get_result('charge_density',None)
    assert(Path(result.locdir)==tmp_path/'scf')
    assert(Path(result.location)==tmp_path/'scf/Waves/wave.out')
    assert(Path(result.outdir)==tmp_path/'scf/Waves')
    assert(tuple(result.wavefunction_grid)==(16,16,16))
    assert(result.processor_grid is None)

    result = scf.get_result('orbitals',None)
    assert(Path(result.h5file)==tmp_path/'scf/Waves/wave.out.h5')

    # a custom restart file name is followed
    scf2 = get_rmg_sim('scf2','Quench Electrons',output_wave_function_file='Out/psi')
    result = scf2.get_result('orbitals',None)
    assert(Path(result.h5file)==tmp_path/'scf2/Out/psi.h5')

    clear_all_sims()
#end def test_get_result


@isolate_nexus_core
def test_incorporate_result(tmp_path):
    nexus_config.runs             = ''
    nexus_config.local_directory  = str(tmp_path)
    nexus_config.remote_directory = str(tmp_path)

    scf = get_rmg_sim('scf','Quench Electrons',write_qmcpack_restart=True)
    result = scf.get_result('charge_density',None)

    waves = tmp_path/'scf/Waves'
    waves.mkdir(parents=True)
    (waves/'wave.out_spin0_kpt0_gridpe0').write_text('restart')
    (waves/'wave.out.h5').write_text('scf orbitals')

    # only an NSCF run can use the charge density
    scf_b = get_rmg_sim('scf_b','Quench Electrons')
    with pytest.raises(ValueError,match='calculation_mode="NSCF"'):
        scf_b.incorporate_result('charge_density',result,scf)
    #end with

    # the restart files are copied, the SCF orbital file is not
    nscf = get_rmg_sim('nscf','NSCF',write_qmcpack_restart=True)
    if shutil.which('rsync') is not None:
        nscf.incorporate_result('charge_density',result,scf)
        assert((tmp_path/'nscf/Waves/wave.out_spin0_kpt0_gridpe0').read_text()=='restart')
        assert(not (tmp_path/'nscf/Waves/wave.out.h5').exists())
        assert((tmp_path/'nscf/Waves/nexus_sync_record').exists())
        assert('input_wave_function_file' not in nscf.input)
    #end if

    # without sync the NSCF run reads the SCF files in place
    nscf2 = get_rmg_sim('nscf2','NSCF',sync_from_scf=False)
    nscf2.incorporate_result('charge_density',result,scf)
    assert(nscf2.input.input_wave_function_file=='../scf/Waves/wave.out')
    assert(not (tmp_path/'nscf2/Waves').exists())

    # the restart needs the same real-space grid
    nscf3 = generate_rmg(
        identifier        = 'nscf',
        path              = 'nscf3',
        job               = job(machine='ws1',cores=1),
        calculation_mode  = 'NSCF',
        wavefunction_grid = (20,20,20),
        )
    with pytest.raises(ValueError,match='wavefunction_grid of the NSCF run does not match'):
        nscf3.incorporate_result('charge_density',result,scf)
    #end with

    with pytest.raises(NotImplementedError,match='ability to incorporate result orbitals'):
        nscf.incorporate_result('orbitals',obj(),scf)
    #end with

    clear_all_sims()
#end def test_incorporate_result
