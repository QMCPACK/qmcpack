import pytest

from . import NexusTestOrder

pytestmark = pytest.mark.order(NexusTestOrder.MEMORY)

from ..memory import memory, resident, stacksize


def test_memory():
    self_mem, child_mem = memory()
    assert(isinstance(self_mem,float))
    assert(isinstance(child_mem,float))
    assert(self_mem > -1e-6)
    assert(child_mem > -1e-6)

    self_mem, child_mem = memory(children=True)
    assert(isinstance(self_mem,float))
    assert(isinstance(child_mem,float))
    assert(self_mem > -1e-6)
    assert(child_mem > -1e-6)
#end def test_memory


def test_resident():
    self_mem, child_mem = resident()
    assert(isinstance(self_mem,float))
    assert(isinstance(child_mem,float))
    assert(self_mem > -1e-6)
    assert(child_mem > -1e-6)

    self_mem, child_mem = resident(children=True)
    assert(isinstance(self_mem,float))
    assert(isinstance(child_mem,float))
    assert(self_mem > -1e-6)
    assert(child_mem > -1e-6)
#end def test_resident


def test_stacksize():
    self_mem, child_mem = stacksize()
    assert(isinstance(self_mem,float))
    assert(isinstance(child_mem,float))
    assert(self_mem > -1e-6)
    assert(child_mem > -1e-6)

    self_mem, child_mem = stacksize(children=True)
    assert(isinstance(self_mem,float))
    assert(isinstance(child_mem,float))
    assert(self_mem > -1e-6)
    assert(child_mem > -1e-6)
#end def test_stacksize
