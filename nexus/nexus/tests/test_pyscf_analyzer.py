import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.PYSCF_ANALYZER)




def test_empty_init():
    from ..pyscf_analyzer import PyscfAnalyzer

    PyscfAnalyzer(None)
#end def test_empty_init
