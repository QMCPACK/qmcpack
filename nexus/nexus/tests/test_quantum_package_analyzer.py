import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.QUANTUM_PACKAGE_ANALYZER)




def test_empty_init():
    from ..quantum_package_analyzer import QuantumPackageAnalyzer

    QuantumPackageAnalyzer(None)
#end def test_empty_init
