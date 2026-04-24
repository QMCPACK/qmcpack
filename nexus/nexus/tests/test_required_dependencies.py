try:
    import pytest
    from . import NexusTestOrder
    pytestmark = pytest.mark.order(NexusTestOrder.REQUIRED_DEPENDENCIES)
    from ..generic import generic_settings
    generic_settings.raise_error = True
except ImportError:
    pass


def test_numpy_available():
    from .. import versions
    assert(versions.numpy_available)
#end def test_numpy_available


# skip this since the rest of the test set actually tells you if it is supported
#def test_numpy_supported():
#    from .. import versions
#    if versions.numpy_available:
#        assert(versions.numpy_supported)
#    #end if
##end def test_numpy_supported
