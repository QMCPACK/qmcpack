
import testing
from testing import value_eq,object_eq
from testing import divert_nexus_log,restore_nexus_log


def test_import():
    import nexus_base
    from nexus_base import nexus_core,nexus_noncore,nexus_core_noncore
    from nexus_base import NexusCore
#end def test_import



def test_namespaces():
    from nexus_base import nexus_core,nexus_core_defaults
    from nexus_base import nexus_noncore,nexus_noncore_defaults
    from nexus_base import nexus_core_noncore,nexus_core_noncore_defaults
    
    assert('runs' in nexus_core_defaults)
    assert('basis_dir' in nexus_noncore_defaults)
    assert('pseudo_dir' in nexus_core_noncore_defaults)

    assert(object_eq(nexus_core,nexus_core_defaults))
    assert(object_eq(nexus_noncore,nexus_noncore_defaults))
    assert(object_eq(nexus_core_noncore,nexus_core_noncore_defaults))
#end def test_namespaces



def test_empty_init():
    from nexus_base import NexusCore
    nc = NexusCore()
#end def test_empty_init



def test_write_splash():
    from nexus_base import NexusCore
    log = divert_nexus_log()
    nc = NexusCore()
    assert(not nc.wrote_splash)
    nc.write_splash()
    assert('Nexus' in log.contents())
    assert('Please cite:' in log.contents())
    assert(nc.wrote_splash)
    restore_nexus_log()
#end def test_write_splash
    


def test_enter_leave():
    import os
    from nexus_base import NexusCore

    tpath = testing.setup_unit_test_output_directory('nexus_base','test_enter_leave')

    cwd = os.getcwd()

    log = divert_nexus_log()

    nc = NexusCore()

    nc.enter(tpath)
    tcwd = os.getcwd()
    assert(tcwd==tpath)
    assert('Entering' in log.contents())
    assert(tpath in log.contents())

    nc.leave()
    assert(os.getcwd()==cwd)

    restore_nexus_log()
#end def test_enter_leave
