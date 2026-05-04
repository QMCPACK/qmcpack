import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.NEXUS_BASE)

from ..generic import generic_settings
generic_settings.raise_error = True

from . import isolate_nexus_core
from ..testing import object_eq



def test_namespaces():
    from ..nexus_base import nexus_core,nexus_core_defaults
    from ..nexus_base import nexus_noncore,nexus_noncore_defaults
    from ..nexus_base import nexus_core_noncore,nexus_core_noncore_defaults
    
    assert('runs' in nexus_core_defaults)
    assert('basis_dir' in nexus_noncore_defaults)
    assert('pseudo_dir' in nexus_core_noncore_defaults)

    assert(object_eq(nexus_core,nexus_core_defaults))
    assert(object_eq(nexus_noncore,nexus_noncore_defaults))
    assert(object_eq(nexus_core_noncore,nexus_core_noncore_defaults))
#end def test_namespaces



def test_empty_init():
    from ..nexus_base import NexusCore
    nc = NexusCore()
#end def test_empty_init


@isolate_nexus_core
def test_write_splash():
    from ..nexus_base import NexusCore

    log = generic_settings.devlog
    nc = NexusCore()
    assert(not nc.wrote_splash)
    nc.write_splash()
    assert('Nexus' in log.contents())
    assert('Please cite:' in log.contents())
    assert(nc.wrote_splash)
#end def test_write_splash
    

@isolate_nexus_core
def test_enter_leave(tmp_path):
    import os
    from ..nexus_base import NexusCore

    tmp_dir = tmp_path / "test_nexus_base_output"
    tmp_dir.mkdir(exist_ok=True)

    cwd = os.getcwd()

    log = generic_settings.devlog

    nc = NexusCore()

    nc.enter(tmp_dir)
    tcwd = os.getcwd()
    assert(tcwd==str(tmp_dir))
    assert('Entering' in log.contents())
    assert(str(tmp_dir) in log.contents())

    nc.leave()
    assert(os.getcwd()==cwd)
#end def test_enter_leave
