import pytest

from ..quantum_package_analyzer import QuantumPackageAnalyzer
from . import NexusTestOrder

pytestmark = pytest.mark.order(NexusTestOrder.QUANTUM_PACKAGE_ANALYZER)

def test_empty_init():
    QuantumPackageAnalyzer(None)
#end def test_empty_init
