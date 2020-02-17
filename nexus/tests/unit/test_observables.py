
import testing
from testing import value_eq,object_eq,text_eq,check_object_eq
from testing import FailedTest,failed


def test_import():
    import observables
    from observables import AttributeProperties,DefinedAttributeBase
    from observables import Observable 
    from observables import MomentumDistribution
#end def test_import



def test_defined_attribute_base():
    from generic import obj,NexusError
    from observables import AttributeProperties,DefinedAttributeBase

    # empty init
    p = AttributeProperties()
    o = DefinedAttributeBase()

    pref = obj(
        assigned        = set(),
        deepcopy        = False,
        default         = None,
        dest            = None,
        name            = None,
        no_default      = False,
        required        = False,
        type            = None,
        )
    assert(check_object_eq(p,pref))
    assert(len(o)==0)

    # init
    p = AttributeProperties(
        default  = 2,
        dest     = 'nest',
        deepcopy = True,
        )
    
    pref = obj(
        assigned        = {'dest', 'default', 'deepcopy'},
        deepcopy        = True,
        default         = 2,
        dest            = 'nest',
        name            = None,
        no_default      = False,
        required        = False,
        type            = None,
        )
    assert(check_object_eq(p,pref))


    # define attributes
    class DA(DefinedAttributeBase):
        None
    #end class DA

    da_attributes = obj(
        a = obj(
            default    = 1,
            ),         
        b = obj(       
            default    = 2,
            type       = int,
            required   = True,
            ),         
        c = obj(       
            dest       = 'nest',
            type       = str,
            ),
        d = obj(
            type       = dict,
            deepcopy   = True,
            no_default = True,
            ),
        nest = obj(
            type       = obj,
            default    = obj,
            ),
        )

    DA.define_attributes(**da_attributes)

    def get_class_dict(cls):
        o = obj()
        o.transfer_from(cls.__dict__)
        for k in list(o.keys()):
            if k.startswith('_'):
                del o[k]
            #end if
        #end for
        return o
    #end def get_class_dict

    o = get_class_dict(DA)

    oref = obj(
        deepcopy_attributes = {'d'},
        required_attributes = {'b'},
        sublevel_attributes = {'c'},
        toplevel_attributes = {'d', 'a', 'b', 'nest'},
        typed_attributes    = {'d', 'c', 'b', 'nest'},
        attribute_definitions = obj(
            a = obj(
                assigned        = {'default'},
                deepcopy        = False,
                default         = 1,
                dest            = None,
                name            = 'a',
                no_default      = False,
                required        = False,
                type            = None,
                ),
            b = obj(
                assigned        = {'required', 'default', 'type'},
                deepcopy        = False,
                default         = 2,
                dest            = None,
                name            = 'b',
                no_default      = False,
                required        = True,
                type            = int,
                ),
            c = obj(
                assigned        = {'dest', 'type'},
                deepcopy        = False,
                default         = None,
                dest            = 'nest',
                name            = 'c',
                no_default      = False,
                required        = False,
                type            = str,
                ),
            d = obj(
                assigned        = {'deepcopy', 'no_default', 'type'},
                deepcopy        = True,
                default         = None,
                dest            = None,
                name            = 'd',
                no_default      = True,
                required        = False,
                type            = dict,
                ),
            nest = obj(
                assigned        = {'default', 'type'},
                deepcopy        = False,
                default         = obj,
                no_default      = False,
                dest            = None,
                name            = 'nest',
                required        = False,
                type            = obj,
                ),
            ),
        )
    
    assert(check_object_eq(o,oref))


    class DA2(DA):
        None
    #end class DA2

    DA2.define_attributes(
        DA,
        e = obj(
            required = True,
            )
        )

    o = get_class_dict(DA)
    assert(check_object_eq(o,oref))

    o2 = get_class_dict(DA2)
    oref.required_attributes.add('e')
    oref.toplevel_attributes.add('e')
    oref.attribute_definitions.e = obj(
        assigned        = {'required'},
        deepcopy        = False,
        default         = None,
        dest            = None,
        name            = 'e',
        no_default      = False,
        required        = True,
        type            = None,
        )
    assert(check_object_eq(o2,oref))


    # empty init
    da = DA()
    assert(len(da)==0)
    assert(not da.check_attributes())


    # set_default_attributes
    da.set_default_attributes()
    
    da_ref = obj(
        a    = 1,
        b    = 2,
        nest = obj(
            c = None,
            ),
        )

    assert(check_object_eq(da,da_ref))
    assert(not da.check_attributes())


    # set_attributes/init
    da = DA(
        a = 2,
        b = 3,
        c = 'hi',
        d = dict(a=2),
        )

    da_ref = obj(
        a = 2,
        b = 3,
        d = {'a': 2},
        nest = obj(
            c = 'hi',
            ),
        )

    assert(check_object_eq(da,da_ref))
    assert(da.check_attributes())


    # set_attribute
    da = DA()

    assert('b' not in da)
    da.set_attribute('b',5)
    assert('b' in da)
    assert(da.b==5)

    try:
        da.set_attribute('unknown',None)
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    da.set_attribute('b',3)
    assert(da.b==3)

    try:
        da.set_attribute('b',3.5)
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try


    # get_attribute
    da = DA()

    da.set_attribute('b',5)
    assert('b' in da)
    assert(da.b==5)

    assert(da.get_attribute('b')==5)
    assert(da.get_attribute('a',None) is None)
    assert(da.get_attribute('c',None) is None)

    try:
        da.get_attribute('unknown')
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    try:
        da.get_attribute('a')
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    try:
        da.get_attribute('c')
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try


    # default values
    class DA_def(DefinedAttributeBase):
        None
    #end class DA_def
    DA_def.set_unassigned_default(None)

    class DA_def2(DA_def):
        None
    #end class DA_def2

    DA_def2.define_attributes(**da_attributes)

    assert(not DefinedAttributeBase.class_has('unassigned_default'))
    assert(DA_def.class_has('unassigned_default'))
    assert(DA_def2.class_has('unassigned_default'))
    assert(DA_def.unassigned_default is None)
    assert(DA_def2.unassigned_default is None)

    da = DA_def2()
    assert('a' not in da)
    assert('c' not in da)
    assert('nest' not in da)
    assert(da.get_attribute('a',None) is None)
    assert(da.get_attribute('c',None) is None)

    try:
        da.get_attribute('a')
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    try:
        da.get_attribute('c')
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    da.set_default_attributes()

    da_ref = obj(
        a = 1,
        b = 2,
        nest = obj(
            c = None,
            )
        )
    
    assert(check_object_eq(da,da_ref))

    da.b = None

    assert(da.get_attribute('a')==1)
    assert(da.get_attribute('b',assigned=False) is None)
    assert(da.get_attribute('c',assigned=False) is None)
    assert(da.get_attribute('a',2)==1)
    assert(da.get_attribute('a',assigned=False)==1)

    try:
        da.get_attribute('b')
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    try:
        da.get_attribute('c')
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

#end def test_defined_attribute_base
