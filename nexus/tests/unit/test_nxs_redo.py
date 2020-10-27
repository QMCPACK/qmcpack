
import testing
from testing import create_file,create_path,execute


def test_redo():
    import os

    tpath = testing.setup_unit_test_output_directory('nxs_redo','test_redo')
    
    exe = testing.executable_path('nxs-redo')

    command = '{} {}'.format(exe,tpath)


    # empty directory
    assert(os.listdir(tpath)==[])

    out,err,rc = execute(command)

    assert('no simulation directories found' in out)


    # directory w/ files, but not nexus simulation directory
    create_file('qmc.in.xml',tpath)

    out,err,rc = execute(command)

    assert('no simulation directories found' in out)

    assert(set(os.listdir(tpath))==set(['qmc.in.xml']))


    # nexus simulation directory
    create_path('sim_qmc',tpath)

    assert(set(os.listdir(tpath))==set(['qmc.in.xml','sim_qmc']))

    out,err,rc = execute(command)

    assert(set(os.listdir(tpath))==set(['attempt1']))


    # nexus simulation directory w/ prior attempt
    create_file('qmc.in.xml',tpath)
    create_path('sim_qmc',tpath)

    assert(set(os.listdir(tpath))==set(['attempt1','qmc.in.xml','sim_qmc']))

    out,err,rc = execute(command)

    assert(set(os.listdir(tpath))==set(['attempt1','attempt2']))
#end def test_redo
