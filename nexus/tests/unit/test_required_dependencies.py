
def test_numpy_available():
    import versions
    assert(versions.numpy_available)
#end def test_numpy_available


# skip this since the rest of the test set actually tells you if it is supported
#def test_numpy_supported():
#    import versions
#    if versions.numpy_available:
#        assert(versions.numpy_supported)
#    #end if
##end def test_numpy_supported
