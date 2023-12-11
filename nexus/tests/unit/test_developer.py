
from testing import failed,FailedTest

def test_import():
    from developer import log,error,warn
    from developer import ci,interact
    from developer import DevBase
    from developer import Void
    from developer import unavailable,available
#end def test_import



def test_unavailable():
    from generic import generic_settings,NexusError
    from developer import Void
    from developer import unavailable,available

    gre = generic_settings.raise_error
    generic_settings.raise_error = True

    try:
        import keyword
    except:
        keyword = unavailable('keyword')
    #end try

    assert(not isinstance(keyword,Void))
    assert(available(keyword))


    try:
        import an_unavailable_module
    except:
        an_unavailable_module = unavailable('an_unavailable_module')
    #end try

    assert(isinstance(an_unavailable_module,Void))
    assert(not available(an_unavailable_module))


    try:
        from an_unavailable_module import a,b,c,d,e,f,g
    except:
        a,b,c,d,e,f,g = unavailable('an_unavailable_module','a','b','c','d','e','f','g')
    #end try

    void_imports = an_unavailable_module,b,c,d,e,f,g

    for v in void_imports:
        assert(isinstance(v,Void))
        assert(not available(v))
    #end for

    operations = [
        dir,
        len,
        repr,
        str,
        complex,
        int,
        float,
        lambda v: v==0,
        lambda v: v!=0,
        lambda v: v>0,
        lambda v: v<0,
        lambda v: v>=0,
        lambda v: v<=0,
        lambda v: v(),
        lambda v: v.a,
        lambda v: v['a'],
        lambda v: setattr(v,'a',0),
        lambda v: getattr(v,'a'),
        lambda v: delattr(v,'a'),
        lambda v: v+0,
        lambda v: v-0,
        lambda v: v*0,
        lambda v: v/0,
        lambda v: v%0,
        lambda v: v&0,
        lambda v: v|0,
        ]
    for op in operations:
        for v in void_imports:
            try:
                op(v)
                raise FailedTest
            except NexusError:
                None
            except FailedTest:
                failed()
            except Exception as e:
                failed(str(e))
            #end try
        #end for
    #end for

    generic_settings.raise_error = gre

#end def test_unavailable
