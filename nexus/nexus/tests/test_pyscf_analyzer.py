import pytest

from ..pyscf_analyzer import PyscfAnalyzer
from . import NexusTestOrder

pytestmark = pytest.mark.order(NexusTestOrder.PYSCF_ANALYZER)

def test_empty_init():
    PyscfAnalyzer(None)
#end def test_empty_init
