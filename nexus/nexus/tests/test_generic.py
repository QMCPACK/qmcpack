import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.GENERIC_OPERATION)

from . import isolate_nexus_core
from ..generic import warn, NexusDevWarning, NexusUserWarning, nxs_deprecate
from ..generic import nxs_print,error
from ..generic import NexusError


def test_logging(tmp_path, capsys):

    # test log
    #   simple message
    s = 'simple message'
    nxs_print(s)
    captured = capsys.readouterr()
    assert(captured.out==s+'\n')

    #   list of items
    items = ['a','b','c',1,2,3]
    nxs_print(*items)
    captured = capsys.readouterr()
    assert(captured.out=='a b c 1 2 3 \n')

    #   message with indentation
    s = 'a message\nwith indentation'
    nxs_print(s,indent='  ')
    captured = capsys.readouterr()
    assert(captured.out=='  a message\n  with indentation\n')

    nxs_print(s,indent='msg: ')
    captured = capsys.readouterr()
    assert(captured.out=='msg: a message\nmsg: with indentation\n')
    
    #   writing to separate log files
    logfile = tmp_path / "fake.log"
    logfile.touch()
    s1 = 'message to log 1'
    s2 = 'message to log 2'
    nxs_print(s1)
    captured = capsys.readouterr()
    assert(captured.out==s1+'\n')
    assert(logfile.read_text()=='')

    with open(logfile, "w") as lf:
        nxs_print(s2, logfile=lf)
    captured = capsys.readouterr()
    assert(captured.out=='')
    assert(logfile.read_text()==s2+'\n')

    # test error
    with pytest.raises(NexusError, match="testing environment"):
        error('testing environment')

    # test error with header
    with pytest.raises(NexusError, match="header error:\ntesting environment"):
        error('testing environment', "header")
#end def test_logging


def test_warn():
    with pytest.warns(NexusUserWarning, match="This is a test warning"):
        warn("This is a test warning", warn_type="user")

    with pytest.warns(NexusDevWarning, match="This is a developer warning"):
        warn("This is a developer warning", warn_type="dev")
#end def test_warn


def test_nxs_deprecate():

    @nxs_deprecate(since="2.3.9", replacement="Some other function")
    def deprecated_function():
        pass

    with pytest.warns(
        DeprecationWarning,
        match="deprecated_function is deprecated as of Nexus version 2.3.9, and will be removed in a future update."
        ):
        deprecated_function()
