
import testing
from testing import divert_nexus,restore_nexus,clear_all_sims
from testing import text_eq



def test_sim():
    import os
    from execute import execute
    from nexus_base import nexus_core
    from test_simulation import get_sim

    tpath = testing.setup_unit_test_output_directory('nxs_sim','test_sim',divert=True)

    nexus_core.runs    = ''
    nexus_core.results = ''
    
    exe = testing.executable_path('nxs-sim')

    sim = get_sim()

    sim.create_directories()

    sim.save_image()

    simp_path = os.path.join(tpath,sim.imlocdir,'sim.p')
    assert(os.path.isfile(simp_path))


    # initial simulation state
    command = '{} show {}'.format(exe,simp_path)

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
    command = '{} complete {}'.format(exe,simp_path)
    out,err,rc = execute(command)

    command = '{} show {}'.format(exe,simp_path)
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
    command = '{} reset {}'.format(exe,simp_path)
    out,err,rc = execute(command)

    command = '{} set setup sent_files submitted {}'.format(exe,simp_path)
    out,err,rc = execute(command)

    command = '{} show {}'.format(exe,simp_path)
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
    command = '{} complete {}'.format(exe,simp_path)
    out,err,rc = execute(command)

    command = '{} unset got_output analyzed {}'.format(exe,simp_path)
    out,err,rc = execute(command)

    command = '{} show {}'.format(exe,simp_path)
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
    restore_nexus()
#end def test_sim
