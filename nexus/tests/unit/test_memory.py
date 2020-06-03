

def test_import():
    from memory import memory,resident,stacksize
#end def test_import



def test_memory():
    from memory import memory

    mem = memory()
    assert(isinstance(mem,float))
    assert(mem>-1e-6)

    mem = memory(children=True)
    assert(isinstance(mem,float))
    assert(mem>-1e-6)
#end def test_memory



def test_resident():
    from memory import resident

    mem = resident()
    assert(isinstance(mem,float))
    assert(mem>-1e-6)

    mem = resident(children=True)
    assert(isinstance(mem,float))
    assert(mem>-1e-6)
#end def test_resident



def test_stacksize():
    from memory import stacksize

    mem = stacksize()
    assert(isinstance(mem,float))
    assert(mem>-1e-6)

    mem = stacksize(children=True)
    assert(isinstance(mem,float))
    assert(mem>-1e-6)
#end def test_stacksize
