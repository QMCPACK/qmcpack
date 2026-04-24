try:
    import pytest
    from . import NexusTestOrder
    pytestmark = pytest.mark.order(NexusTestOrder.PYSCF_ANALYZER)
    from ..generic import generic_settings
    generic_settings.raise_error = True
except ImportError:
    pass


def test_import():
    from ..pyscf_analyzer import PyscfAnalyzer
#end def test_import



def test_empty_init():
    from ..pyscf_analyzer import PyscfAnalyzer

    PyscfAnalyzer(None)
#end def test_empty_init
