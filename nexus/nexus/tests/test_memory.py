import pytest

from ..memory import memory, resident, stacksize
from . import NexusTestOrder

pytestmark = pytest.mark.order(NexusTestOrder.MEMORY)

def test_memory():
    mem = memory()
    assert(isinstance(mem,float))
    assert(mem>-1e-6)

    mem = memory(children=True)
    assert(isinstance(mem,float))
    assert(mem>-1e-6)
#end def test_memory



def test_resident():
    mem = resident()
    assert(isinstance(mem,float))
    assert(mem>-1e-6)

    mem = resident(children=True)
    assert(isinstance(mem,float))
    assert(mem>-1e-6)
#end def test_resident



def test_stacksize():
    mem = stacksize()
    assert(isinstance(mem,float))
    assert(mem>-1e-6)

    mem = stacksize(children=True)
    assert(isinstance(mem,float))
    assert(mem>-1e-6)
#end def test_stacksize
