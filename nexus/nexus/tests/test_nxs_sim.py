import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.NXS_SIM)

from ..generic import generic_settings
generic_settings.raise_error = True

import sys
from pathlib import Path
from . import isolate_nexus_core
from ..testing import clear_all_sims
from ..testing import execute,text_eq


@isolate_nexus_core(needs_tmp_path=True)
def test_sim(tmp_path):
    from ..nexus_base import nexus_core
    from .test_simulation_module import get_sim

    tmp_dir = tmp_path / "test_nxs_sim_output"
    tmp_dir.mkdir()

    nexus_core.local_directory  = tmp_dir
    nexus_core.remote_directory = tmp_dir
    nexus_core.file_locations = nexus_core.file_locations + [tmp_dir]

    nexus_core.runs    = ''
    nexus_core.results = ''
    
    exe = Path(__file__).parent.parent / "bin/nxs-sim"

    sim = get_sim()

    sim.create_directories()

    sim.save_image()

    simp_path = (tmp_dir / sim.imlocdir / 'sim.p').resolve()
    assert(simp_path.is_file())


    # initial simulation state
    command = sys.executable+' {} show {}'.format(exe,simp_path)

    out,err,rc = execute(command)

    out_ref = '''
        {}
        setup        0
        sent_files   0
        submitted    0
        finished     0
        failed       0
        got_output   0
        analyzed     0
        '''.format(simp_path)

    assert(text_eq(out,out_ref))


    # final simulation state
    command = sys.executable+' {} complete {}'.format(exe,simp_path)
    out,err,rc = execute(command)

    command = sys.executable+' {} show {}'.format(exe,simp_path)
    out,err,rc = execute(command)

    out_ref = '''
        {}
        setup        1
        sent_files   1
        submitted    1
        finished     1
        failed       0
        got_output   1
        analyzed     1
        '''.format(simp_path)

    assert(text_eq(out,out_ref))


    # intermediate simulation state 1
    command = sys.executable+' {} reset {}'.format(exe,simp_path)
    out,err,rc = execute(command)

    command = sys.executable+' {} set setup sent_files submitted {}'.format(exe,simp_path)
    out,err,rc = execute(command)

    command = sys.executable+' {} show {}'.format(exe,simp_path)
    out,err,rc = execute(command)

    out_ref = '''
        {}
        setup        1
        sent_files   1
        submitted    1
        finished     0
        failed       0
        got_output   0
        analyzed     0
        '''.format(simp_path)


    # intermediate simulation state 2
    command = sys.executable+' {} complete {}'.format(exe,simp_path)
    out,err,rc = execute(command)

    command = sys.executable+' {} unset got_output analyzed {}'.format(exe,simp_path)
    out,err,rc = execute(command)

    command = sys.executable+' {} show {}'.format(exe,simp_path)
    out,err,rc = execute(command)

    out_ref = '''
        {}
        setup        1
        sent_files   1
        submitted    1
        finished     1
        failed       0
        got_output   0
        analyzed     0
        '''.format(simp_path)

    assert(text_eq(out,out_ref))

    clear_all_sims()
#end def test_sim
