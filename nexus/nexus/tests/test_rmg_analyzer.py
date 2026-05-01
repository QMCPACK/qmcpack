import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.RMG_ANALYZER)

from ..generic import generic_settings
generic_settings.raise_error = True



def test_empty_init():
    from ..rmg_analyzer import RmgAnalyzer

    RmgAnalyzer()
#end def test_empty_init
