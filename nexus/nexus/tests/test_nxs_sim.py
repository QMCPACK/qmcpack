import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.NXS_SIM)

import sys
from . import isolate_nexus_core, TEST_DIR
from ..testing import clear_all_sims
from ..testing import execute,text_eq


@isolate_nexus_core
def test_sim(tmp_path):
    from ..nexus_base import NEXUS_CONFIG
    from .test_simulation_module import get_sim

    NEXUS_CONFIG.local_directory  = str(tmp_path)
    NEXUS_CONFIG.remote_directory = str(tmp_path)
    NEXUS_CONFIG.file_locations = NEXUS_CONFIG.file_locations + [str(tmp_path)]

    NEXUS_CONFIG.runs    = ''
    NEXUS_CONFIG.results = ''
    
    exe = TEST_DIR.parent / "bin/nxs-sim"

    sim = get_sim()

    sim.create_directories()

    sim.save_image()

    simp_path = (tmp_path / sim.imlocdir / 'sim.p').resolve()
    assert(simp_path.is_file())


    # initial simulation state
    command = sys.executable+f' {exe} show {simp_path}'

    out,err,rc = execute(command)

    out_ref = f'''
        {simp_path}
        setup        0
        sent_files   0
        submitted    0
        finished     0
        failed       0
        got_output   0
        analyzed     0
        '''

    assert(text_eq(out,out_ref))


    # final simulation state
    command = sys.executable+f' {exe} complete {simp_path}'
    out,err,rc = execute(command)

    command = sys.executable+f' {exe} show {simp_path}'
    out,err,rc = execute(command)

    out_ref = f'''
        {simp_path}
        setup        1
        sent_files   1
        submitted    1
        finished     1
        failed       0
        got_output   1
        analyzed     1
        '''

    assert(text_eq(out,out_ref))


    # intermediate simulation state 1
    command = sys.executable+f' {exe} reset {simp_path}'
    out,err,rc = execute(command)

    command = sys.executable+f' {exe} set setup sent_files submitted {simp_path}'
    out,err,rc = execute(command)

    command = sys.executable+f' {exe} show {simp_path}'
    out,err,rc = execute(command)

    out_ref = f'''
        {simp_path}
        setup        1
        sent_files   1
        submitted    1
        finished     0
        failed       0
        got_output   0
        analyzed     0
        '''


    # intermediate simulation state 2
    command = sys.executable+f' {exe} complete {simp_path}'
    out,err,rc = execute(command)

    command = sys.executable+f' {exe} unset got_output analyzed {simp_path}'
    out,err,rc = execute(command)

    command = sys.executable+f' {exe} show {simp_path}'
    out,err,rc = execute(command)

    out_ref = f'''
        {simp_path}
        setup        1
        sent_files   1
        submitted    1
        finished     1
        failed       0
        got_output   0
        analyzed     0
        '''

    assert(text_eq(out,out_ref))

    clear_all_sims()
#end def test_sim
