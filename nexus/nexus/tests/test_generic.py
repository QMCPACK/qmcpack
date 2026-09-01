import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.GENERIC_OPERATION)

from . import isolate_nexus_core, FakeLog
from ..generic import warn, NexusDevWarning, NexusUserWarning, nxs_deprecate

@isolate_nexus_core
def test_logging():
    from ..generic import log,error
    from ..generic import generic_settings,NexusError

    logfile = generic_settings.devlog

    # test log
    #   simple message
    s = 'simple message'
    logfile.reset()
    log(s)
    assert(logfile.s==s+'\n')

    #   list of items
    items = ['a','b','c',1,2,3]
    logfile.reset()
    log(*items)
    assert(logfile.s=='a b c 1 2 3 \n')

    #   message with indentation
    s = 'a message\nwith indentation'
    logfile.reset()
    log(s,indent='  ')
    assert(logfile.s=='  a message\n  with indentation\n')

    logfile.reset()
    log(s,indent='msg: ')
    assert(logfile.s=='msg: a message\nmsg: with indentation\n')
    
    #   writing to separate log files
    logfile2 = FakeLog()
    s1 = 'message to log 1'
    s2 = 'message to log 2'
    logfile.reset()
    logfile2.reset()
    log(s1)
    assert(logfile.s==s1+'\n')
    assert(logfile2.s=='')

    logfile.reset()
    logfile2.reset()
    log(s2,logfile=logfile2)
    assert(logfile.s=='')
    assert(logfile2.s==s2+'\n')

    # test error
    with pytest.raises(NexusError, match="testing environment"):
        error('testing environment')

    # test error with header
    with pytest.raises(NexusError, match="header error:\ntesting environment"):
        error('testing environment', "header")
#end def test_logging


@isolate_nexus_core
def test_warn():
    with pytest.warns(NexusUserWarning, match="This is a test warning"):
        warn("This is a test warning", warn_type="user")

    with pytest.warns(NexusDevWarning, match="This is a developer warning"):
        warn("This is a developer warning", warn_type="dev")
#end def test_warn


@isolate_nexus_core
def test_nxs_deprecate():

    @nxs_deprecate(since="2.3.9", replacement="Some other function")
    def deprecated_function():
        pass

    with pytest.warns(
        DeprecationWarning,
        match="deprecated_function is deprecated as of Nexus version 2.3.9, and will be removed in a future update."
        ):
        deprecated_function()
