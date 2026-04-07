try:
    import pytest
    from . import NexusTestOrder
    pytestmark = pytest.mark.order(NexusTestOrder.OPTIONAL_DEPENDENCIES)
    from ..generic import generic_settings
    generic_settings.raise_error = True
except ImportError:
    pass


def test_scipy_available():
    try:
        import pytest
        scipy = pytest.importorskip("scipy")
    except ImportError:
        pass
    from .. import versions
    assert(versions.scipy_available)
#end def test_scipy_available


def test_h5py_available():
    try:
        import pytest
        h5py = pytest.importorskip("h5py")
    except ImportError:
        pass
    from .. import versions
    assert(versions.h5py_available)
#end def test_h5py_available


def test_matplotlib_available():
    try:
        import pytest
        matplotlib = pytest.importorskip("matplotlib")
    except ImportError:
        pass
    from .. import versions
    assert(versions.matplotlib_available)
#end def test_matplotlib_available


def test_pydot_available():
    try:
        import pytest
        pydot = pytest.importorskip("pydot")
    except ImportError:
        pass
    from .. import versions
    assert(versions.pydot_available)
#end def test_pydot_available


def test_spglib_available():
    try:
        import pytest
        spglib = pytest.importorskip("spglib")
    except ImportError:
        pass
    from .. import versions
    assert(versions.spglib_available)
#end def test_spglib_available


def test_pycifrw_available():
    try:
        import pytest
        pycifrw = pytest.importorskip("pycifrw")
    except ImportError:
        pass
    from .. import versions
    assert(versions.pycifrw_available)
#end def test_pycifrw_available


def test_seekpath_available():
    try:
        import pytest
        seekpath = pytest.importorskip("seekpath")
    except ImportError:
        pass
    from .. import versions
    assert(versions.seekpath_available)
#end def test_seekpath_available
