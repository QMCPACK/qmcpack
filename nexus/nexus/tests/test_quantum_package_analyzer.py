import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.QUANTUM_PACKAGE_ANALYZER)

from ..generic import generic_settings
generic_settings.raise_error = True



def test_empty_init():
    from ..quantum_package_analyzer import QuantumPackageAnalyzer

    QuantumPackageAnalyzer(None)
#end def test_empty_init
