try:
    import pytest
    from . import NexusTestOrder
    pytestmark = pytest.mark.order(NexusTestOrder.EXECUTE)
except ImportError:
    pass


def test_import():
    from ..execute import execute
#end def test_import



def test_execute():
    from ..execute import execute

    command = 'echo test'

    out,err,rc = execute(command,skip=True)
    assert(out=='')
    assert(err=='')
    assert(rc==0)

    out,err,rc = execute(command)
    assert(out.strip()=='test')
    assert(err.strip()=='')
    assert(rc==0)
#end def test_execute
