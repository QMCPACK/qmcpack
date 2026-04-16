try:
    import pytest
    from . import NexusTestOrder
    pytestmark = pytest.mark.order(NexusTestOrder.PYSCF_ANALYZER)
except ImportError:
    pass


def test_import():
    from ..pyscf_analyzer import PyscfAnalyzer
#end def test_import



def test_empty_init():
    from ..pyscf_analyzer import PyscfAnalyzer

    PyscfAnalyzer(None)
#end def test_empty_init
