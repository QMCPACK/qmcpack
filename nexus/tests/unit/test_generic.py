
import testing
from testing import failed,FailedTest
from testing import divert_nexus_log,restore_nexus_log,FakeLog
from testing import value_eq,object_eq,object_neq


def test_logging():
    from generic import log,message,warn,error
    from generic import generic_settings,NexusError

    # send messages to object rather than stdout
    divert_nexus_log()

    logfile = generic_settings.devlog

    # test log
    #   simple message
    s = 'simple message'
    logfile.reset()
    log(s)
    assert(logfile.s==s+'\n')

    #   list of items
    items = ['a','b','c',1,2,3]
    logfile.reset()
    log(*items)
    assert(logfile.s=='a b c 1 2 3 \n')

    #   message with indentation
    s = 'a message\nwith indentation'
    logfile.reset()
    log(s,indent='  ')
    assert(logfile.s=='  a message\n  with indentation\n')

    logfile.reset()
    log(s,indent='msg: ')
    assert(logfile.s=='msg: a message\nmsg: with indentation\n')
    
    #   writing to separate log files
    logfile2 = FakeLog()
    s1 = 'message to log 1'
    s2 = 'message to log 2'
    logfile.reset()
    logfile2.reset()
    log(s1)
    assert(logfile.s==s1+'\n')
    assert(logfile2.s=='')

    logfile.reset()
    logfile2.reset()
    log(s2,logfile=logfile2)
    assert(logfile.s=='')
    assert(logfile2.s==s2+'\n')
    

    # test warn
    logfile.reset()
    s = 'this is a warning'
    warn(s)
    so = '''
  warning:
    this is a warning
'''
    assert(logfile.s==so)
    logfile.reset()
    s = 'this is a warning'
    warn(s,header='Special')
    so = '''
  Special warning:
    this is a warning
'''
    assert(logfile.s==so)

    # test error
    #   in testing environment, should raise an error
    try:
        error('testing environment')
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try
    #   in standard/user environment, should print message
    generic_settings.raise_error = False
    logfile.reset()
    error('this is an error',header='User',exit=False)
    so = '''
  User error:
    this is an error
'''
    assert(logfile.s==so)
    generic_settings.raise_error = True
    #   in testing environment, should raise an error
    try:
        error('testing environment')
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    restore_nexus_log()

#end def test_logging



def test_intrinsics():
    # test object_interface functions
    import os
    from generic import obj,object_interface
    from generic import generic_settings,NexusError
    from numpy import array,bool_

    tpath = testing.setup_unit_test_output_directory('generic','test_intrinsics')

    # test object set/get
    # make a simple object
    o = obj()
    o.a = 1
    o.b = 'b'
    o['c'] = (1,1,1)
    o[3,4,5] = (5,6,7)

    # test member values
    assert(o.a==1)
    assert(o.b=='b')
    assert(o.c==(1,1,1))
    assert(o[3,4,5]==(5,6,7))

    # test member presence and length
    assert('a' in o)
    assert(2 not in o)
    assert(len(o)==4)
    del o.a
    assert('a' not in o)
    assert(len(o)==3)
    try:
        del o.d
        raise FailedTest
    except AttributeError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    # test add/access failures
    try: # add unhashable type
        o[{1:2,3:4}] = 5
        raise FailedTest
    except TypeError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try
    try: # access missing member
        v = o.d
        raise FailedTest
    except AttributeError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try
    try: # access missing member
        v = o['d']
        raise FailedTest
    except KeyError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    # test iterability
    # test list-like iterability
    l = list()
    for v in o:
        l.append(v)
    #end for
    l = set(l)
    l2 = {'b', (1, 1, 1), (5, 6, 7)}
    assert(l2==l)

    # test dict-like iterability
    d = dict()
    for k,v in o.items():
        d[k] = v
    #end for
    o2 = obj()
    o2.__dict__ = d
    assert(object_eq(o,o2))
    assert(set(o.keys())==set(d.keys()))
    assert(set(o.values())==set(d.values()))
    assert(set(o.items())==set(d.items()))

    # test repr
    ro = '''
  b                     str                 
  c                     tuple               
  (3, 4, 5)             tuple               
'''
    assert(repr(o)==ro[1:])
    o2 = obj(
        a = obj(a1=1,a2=2,a3=3),
        b = obj(b1='b1',b2='b2'),
        c = obj(c1=5,c2=('a',3,4)),
        d = array([3,4,5],dtype=int),
        )
    ro2 = '''
  a                     obj                 
  b                     obj                 
  c                     obj                 
  d                     ndarray             
'''
    assert(repr(o2)==ro2[1:])

    # test str
    so = '''
  b               = b
  c               = (1, 1, 1)
  (3, 4, 5)       = (5, 6, 7)
'''
    assert(str(o)==so[1:])
    so2 = '''
  d               = [3 4 5]
  a
    a1              = 1
    a2              = 2
    a3              = 3
  end a
  b
    b1              = b1
    b2              = b2
  end b
  c
    c1              = 5
    c2              = ('a', 3, 4)
  end c
'''
    assert(str(o2)==so2[1:])
    o3 = o2

    # test tree
    #   not committed to output, only check execution
    assert(isinstance(o2.tree(),str))
    assert(isinstance(o2.tree(depth=1),str))
    assert(isinstance(o2.tree(types=True),str))
    assert(isinstance(o2.tree(all=True),str))
    assert(isinstance(o2.tree(nindent=2),str))

    # test deepcopy
    o2 = o.copy()
    assert(id(o)!=id(o2))
    assert(object_eq(o,o2))
    o2.a=1
    assert(object_neq(o,o2))
    
    # test eq
    assert(o==o2)
    o4 = o3.copy()
    v = o3==o4
    assert(isinstance(v,bool_))
    assert(bool(v))
    assert(object_eq(o3,o4))

    # test clear
    o2.clear()
    assert(len(o2)==0)
    assert('a' not in o2)

    # test save/load
    save_file = os.path.join(tpath,'o.p')
    o.save(save_file)
    o2 = obj()
    o2.load(save_file)
    assert(object_eq(o,o2))

    # test class-level set/get methods
    class objint(object_interface):
        a = 1
        b = 'b'
        c = (1,1,1)
    #end class objint

    # test class_has
    assert(objint.class_has('c'))

    # test class_get
    assert(objint.class_get('c')==(1,1,1))
    try:
        val = objint.class_get('d')
        raise FailedTest
    except AttributeError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    # test class_keys
    ck = objint.class_keys()
    for v in 'abc':
        assert(v in ck)
    #end for

    # test class_set_single
    objint.class_set_single('d',1.34)
    assert(objint.class_has('d'))

    # test class_set
    objint.class_set(
        e = 45,
        f = 'a phrase',
        g = {4:6,'a':2}
        )
    for v in 'efg':
        assert(objint.class_has(v))
    #end for

    # test class_set_optional
    objint.class_set_optional(h=2)
    assert(objint.class_has('h'))
    assert(objint.class_get('h')==2)
    objint.class_set_optional(a=6)
    assert(objint.class_get('a')==1)
    

    # test logging functions

    # test open log, write, close log
    o = obj()
    o.open_log('log.out')
    s = 'log output'
    o.write(s)
    o.close_log()
    f = open('log.out','r')
    so = f.read()
    f.close()
    os.remove('log.out')
    assert(so==s)

    # send messages to object rather than stdout
    divert_nexus_log()
    logfile = object_interface._logfile

    #   simple message
    class DerivedObj(obj):
        None
    #end class DerivedObj
    o = DerivedObj()
    s = 'simple message'
    logfile.reset()
    o.log(s)
    assert(logfile.s==s+'\n')

    #   list of items
    items = ['a','b','c',1,2,3]
    logfile.reset()
    o.log(*items)
    assert(logfile.s=='a b c 1 2 3 \n')

    #   message with indentation
    s = 'a message\nwith indentation'
    logfile.reset()
    o.log(s,indent='  ')
    assert(logfile.s=='  a message\n  with indentation\n')

    logfile.reset()
    o.log(s,indent='msg: ')
    assert(logfile.s=='msg: a message\nmsg: with indentation\n')
    
    #   writing to separate log files
    logfile2 = FakeLog()
    s1 = 'message to log 1'
    s2 = 'message to log 2'
    logfile.reset()
    logfile2.reset()
    o.log(s1)
    assert(logfile.s==s1+'\n')
    assert(logfile2.s=='')

    logfile.reset()
    logfile2.reset()
    o.log(s2,logfile=logfile2)
    assert(logfile.s=='')
    assert(logfile2.s==s2+'\n')
    

    # test warn
    logfile.reset()
    s = 'this is a warning'
    o.warn(s)
    so = '''
  DerivedObj warning:
    this is a warning
'''
    assert(logfile.s==so)
    logfile.reset()
    s = 'this is a warning'
    o.warn(s,header='Special')
    so = '''
  Special warning:
    this is a warning
'''
    assert(logfile.s==so)

    # test error
    #   in testing environment, should raise an error
    try:
        o.error('testing environment')
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try
    #   in standard/user environment, should print message
    generic_settings.raise_error = False
    logfile.reset()
    o.error('this is an error',exit=False)
    so = '''
  DerivedObj error:
    this is an error
'''
    assert(logfile.s==so)
    logfile.reset()
    o.error('this is an error',header='User',exit=False)
    so = '''
  User error:
    this is an error
'''
    assert(logfile.s==so)
    generic_settings.raise_error = True
    #   in testing environment, should raise an error
    try:
        o.error('testing environment')
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    # restore logging function
    restore_nexus_log()

#end def test_intrinsics



def test_extensions():
    # test obj functions
    from generic import obj,NexusError

    # make a simple object
    o = obj(
        a = 1,
        b = 'b',
        c = (1,1,1),
        )
    o[3,4,5] = (5,6,7)

    # test member values
    assert(o.a==1)
    assert(o.b=='b')
    assert(o.c==(1,1,1))
    assert(o[3,4,5]==(5,6,7))

    # test member presence and length
    assert('a' in o)
    assert(2 not in o)
    assert(len(o)==4)

    # test list interface
    vals = [3,'t',6.4,(4,3,2)]
    l = list()
    lo = obj()
    for v in vals:
        l.append(v)
        lo.append(v)
    #end for
    assert(len(l)==len(lo))
    for i in range(len(l)):
        assert(l[i]==lo[i])
    #end for


    # test representations
    # test list representation
    l2 = lo.list()
    assert(isinstance(l2,list))
    assert(l==l2)
    l2 = lo.list_optional(1,3)
    assert(l2==['t',(4,3,2)])
    l2 = o.list_optional('b',(3,4,5))
    assert(l2==['b',(5,6,7)])

    # test tuple representation
    t = lo.tuple()
    assert(isinstance(t,tuple))
    assert(t==tuple(l))
    d = dict(
        a = 1,
        b = 'b',
        c = (1,1,1),
        )
    d[3,4,5] = (5,6,7)

    # test dict representation
    do = o.dict()
    assert(isinstance(do,dict))
    assert(do==d)
    d2 = d.copy()
    d2['d'] = d
    o2 = o.copy()
    o2.d = o
    d2o = o2.to_dict()
    assert(d2o==d2)

    # test obj representation
    o2 = o.obj()
    assert(isinstance(o2,obj))
    assert(id(o2)!=id(o))
    assert(object_eq(o2,o))
    o2 = o.copy().to_obj()
    assert(object_eq(o2,o))
    
    # test list extensions
    # test first
    assert(lo.first()==lo[0])

    # test last
    assert(lo.last()==lo[3])

    # test select_random
    v = lo.select_random()
    assert(v in lo.list())


    # test dict extensions
    # test random_key
    k = o.random_key()
    assert(k in o)
    o2 = obj()
    assert(o2.random_key() is None)

    # test set
    o2 = o.copy()
    o2.set(
        b = 'b2',
        d = ('a','b','c'),
        )
    assert(o2.b=='b2')
    assert(o2.d==tuple('abc'))
    o1 = obj(a=1,b=2)
    o2 = obj(c=3,d=4)
    o3 = obj(e=5,f=6)
    o4 = obj()
    o4.set(o1,o2,o3)
    for on in (o1,o2,o3):
        for k,v in on.items():
            assert(o4[k]==v)
        #end for
    #end for

    # test set optional
    o2 = o.copy()
    o2.set_optional(
        b = 'b2',
        d = ('a','b','c'),
        )
    assert(o2.b=='b')
    assert(o2.d==tuple('abc'))
    o1 = obj(a=1,b=2)
    o2 = obj(c=3,d=4)
    o3 = obj(e=5,f=6)
    o4 = obj()
    o4.set_optional(o1,o2,o3)
    for on in (o1,o2,o3):
        for k,v in on.items():
            assert(o4[k]==v)
        #end for
    #end for

    # test get
    assert(o.get('c')==(1,1,1))
    assert('d' not in o)
    assert(o.get('d') is None)

    # test get optional (identical to get)
    assert(o.get_optional('c')==(1,1,1))
    assert(o.get_optional('d') is None)

    # test get required
    assert(o.get_required('c')==(1,1,1))
    try:
        val = o.get_required('d')
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    # test delete
    o2 = o.copy()
    assert(o2.delete('c')==(1,1,1))
    assert('c' not in o2)
    keys = 'a','b','c',(3,4,5)
    vals = [1,'b',(1,1,1),(5,6,7)]
    o2 = o.copy()
    assert(o2.delete(*keys)==vals)
    assert(len(o2)==0)
    for k in keys:
        assert(k not in o2)
    #end for
    o2 = o.copy()
    assert(o2.delete(keys)==vals)
    assert(len(o2)==0)
    for k in keys:
        assert(k not in o2)
    #end for
    o2 = o.copy()
    try:
        o2.delete('a','d')
        raise FailedTest
    except KeyError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    # test delete optional
    o2 = o.copy()
    o2.delete_optional('c')
    assert('c' not in o2)
    assert('d' not in o2)
    o2.delete_optional('d')
    assert('d' not in o2)

    # test delete required
    o2 = o.copy()
    o2.delete_required('c')
    assert('c' not in o2)
    try:
        o2.delete_required('d')
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    # test add
    o2 = obj()
    o2.add('a',1)
    try: # add unhashable type
        o2.add([3,4,5],6)
        raise FailedTest
    except TypeError:
        None
    #end if
    
    # test add optional
    o2 = obj(a=1)
    o2.add_optional('a',2)
    o2.add_optional('b',3)
    assert(o2.a==1)
    assert(o2.b==3)


    # test transfer/copy/move functions
    dref = dict(
        a = 1,
        b = 'b',
        c = (1,1,1),
        d = dict(
            e = 5.4,
            f = [3.3,4.5],
            ),
        )
    dref[3,4,5] = (5,6,7)

    # test transfer from
    o = obj()
    o.transfer_from(dref)
    assert(o.to_dict()==dref)
    assert(id(o.d)==id(dref['d']))

    o = obj()
    o.transfer_from(dref,copy=True)
    assert(o.to_dict()==dref)
    assert(id(o.d)!=id(dref['d']))

    osmall = obj(b='b',c=(1,1,1))
    osmall[3,4,5] = (5,6,7)
    o = obj()
    oref = obj(dref)
    assert(oref.to_dict()==dref)
    o.transfer_from(oref,keys=['b','c',(3,4,5)])
    assert(object_eq(o,osmall))

    o = obj(a=6,b=7)
    o.transfer_from(oref,overwrite=False)
    assert(object_neq(o,oref))
    o.transfer_from(oref,overwrite=True)
    assert(object_eq(o,oref))
    assert(o.to_dict()==dref)

    o = obj()
    try:
        o.transfer_from(oref,keys=['a','x'])
        raise FailedTest
    except KeyError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    # test transfer to
    o = obj()
    oref.transfer_to(o)
    assert(object_eq(o,oref))
    assert(id(o.d)==id(oref.d))

    o = obj()
    oref.transfer_to(o,copy=True)
    assert(object_eq(o,oref))
    assert(id(o.d)!=id(oref.d))

    o = obj()
    oref.transfer_to(o,keys=['b','c',(3,4,5)])
    assert(object_eq(o,osmall))

    o = obj(a=6,b=7)
    oref.transfer_to(o,overwrite=False)
    assert(object_neq(o,oref))
    oref.transfer_to(o,overwrite=True)
    assert(object_eq(o,oref))

    o = obj()
    try:
        oref.transfer_to(o,keys=['a','x'])
        raise FailedTest
    except KeyError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    # test move from
    d2 = dref.copy()
    o = obj()
    o.move_from(d2)
    assert(len(d2)==0)
    assert(object_eq(o,oref))

    o2 = oref.copy()
    o = obj()
    o.move_from(o2)
    assert(len(o2)==0)
    assert(object_eq(o,oref))

    osmall2 = oref.copy()
    del osmall2.b
    del osmall2.c
    del osmall2[3,4,5]
    o2 = oref.copy()
    o = obj()
    o.move_from(o2,keys=['b','c',(3,4,5)])
    assert(object_eq(o,osmall))
    assert(object_eq(o2,osmall2))

    o2 = oref.copy()
    o = obj()
    o.move_from_optional(o2,keys=['b','c',(3,4,5),'alpha','beta'])
    assert(object_eq(o,osmall))
    assert(object_eq(o2,osmall2))

    o2 = oref.copy()
    o = obj()
    try:
        o.move_from(o2,keys=['a','x'])
        raise FailedTest
    except KeyError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    # test move to
    o2 = oref.copy()
    d = dict()
    o2.move_to(d)
    assert(len(o2)==0)
    assert(d==dref)

    o2 = oref.copy()
    o = obj()
    o2.move_to(o)
    assert(len(o2)==0)
    assert(object_eq(o,oref))

    o2 = oref.copy()
    o = obj()
    o2.move_to(o,keys=['b','c',(3,4,5)])
    assert(object_eq(o,osmall))
    assert(object_eq(o2,osmall2))

    o2 = oref.copy()
    o = obj()
    o2.move_to_optional(o,keys=['b','c',(3,4,5),'alpha','beta'])
    assert(object_eq(o,osmall))
    assert(object_eq(o2,osmall2))

    o2 = oref.copy()
    o = obj()
    try:
        o2.move_to(o,keys=['a','x'])
        raise FailedTest
    except KeyError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    # test copy from
    o = obj()
    o.copy_from(dref)
    assert(o.to_dict()==dref)
    assert(id(o.d)!=id(dref['d']))

    o = obj()
    o.copy_from(dref,deep=False)
    assert(o.to_dict()==dref)
    assert(id(o.d)==id(dref['d']))

    osmall = obj(b='b',c=(1,1,1))
    osmall[3,4,5] = (5,6,7)
    o = obj()
    oref = obj(dref)
    assert(oref.to_dict()==dref)
    o.copy_from(oref,keys=['b','c',(3,4,5)])
    assert(object_eq(o,osmall))

    o = obj()
    try:
        o.copy_from(oref,keys=['a','x'])
        raise FailedTest
    except KeyError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    # test copy to
    o = obj()
    oref.copy_to(o)
    assert(object_eq(o,oref))
    assert(id(o.d)!=id(oref.d))

    o = obj()
    oref.copy_to(o,deep=False)
    assert(object_eq(o,oref))
    assert(id(o.d)==id(oref.d))

    o = obj()
    oref.copy_to(o,keys=['b','c',(3,4,5)])
    assert(object_eq(o,osmall))

    o = obj()
    try:
        oref.copy_to(o,keys=['a','x'])
        raise FailedTest
    except KeyError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    # test extract
    o = oref.copy()
    o2 = o.extract()
    assert(len(o)==0)
    assert(object_eq(o2,oref))

    o = oref.copy()
    o2 = o.extract(['b','c',(3,4,5)])
    assert(object_eq(o2,osmall))
    assert(object_eq(o,osmall2))

    o = oref.copy()
    o2 = o.extract_optional(['b','c',(3,4,5),'alpha','beta'])
    assert(object_eq(o2,osmall))
    assert(object_eq(o,osmall2))


    # test check_required
    oref.check_required(['a','d',(3,4,5)])
    try:
        oref.check_required(['alpha','beta'])
        raise FailedTest
    except NexusError:
        None
    except FailedTest:
        failed()
    except Exception as e:
        failed(str(e))
    #end try

    # test check_types
    types = dict(
        a = int,
        b = str,
        c = tuple,
        d = dict,
        )
    types[3,4,5] = tuple
    oref.check_types(types)

    types['b'] = int
    try:
        oref.check_types(types)
        raise FailedTest
    except:
        None
    #end try

    types['b'] = str
    types['alpha'] = float
    types['beta'] = list
    oref.check_types_optional(types)

    # test shallow_copy
    class DerivedObj(obj):
        None
    #end class DerivedObj
    do = DerivedObj(
        a = 1,
        b = 'b',
        c = (1,1,1),
        )
    do[3,4,5] = (5,6,7)
    do2 = do.shallow_copy()
    assert(isinstance(do2,DerivedObj))
    assert(object_eq(do2,do))
    assert(id(do2.c)==id(do.c))

    # test inverse
    oi = do.inverse()
    assert(set(oi.keys())==set(do.values()))
    assert(set(oi.values())==set(do.keys()))

    assert(oi[1]=='a')
    assert(oi.b=='b')
    assert(oi[1,1,1]=='c')
    assert(oi[5,6,7]==(3,4,5))

    # test path operations
    # test path exists
    o2 = obj()
    o2.this = obj()
    o2.this.new = obj()
    o2.this.new.variable = 'here'
    path1 = ['this','new','variable']
    path2 = 'this/new/variable'
    assert(o2.path_exists(path1))
    assert(o2.path_exists(path2))

    # test set path
    o3 = obj()
    o3.set_path(path1,'here')
    assert(o3.path_exists(path1))
    assert(o3.path_exists(path2))
    assert(object_eq(o2,o3))
    o4 = obj()
    o4.set_path(path2,'here')
    assert(o4.path_exists(path1))
    assert(o4.path_exists(path2))
    assert(object_eq(o3,o4))
    
    # test get path
    assert(o2.get_path(path1)=='here')
    assert(o2.get_path(path2)=='here')

    # test serial
    o = obj(
        a = obj(
            a0 = 0,
            a1 = obj(
                a10 = 1,
                ),
            ),
        b = obj(
            b0 = 0,
            b1 = obj(
                b10 = 1,
                ),
            b2 = obj(
                b20 = obj(
                    b200 = 2,
                    ),
                b21 = obj(
                    b210 = obj(
                        b2100 = 3,
                        ),
                    ),
                ),
            ),
        c = obj(
            c0 = 0,
            c1 = obj(
                c10 = 1,
                ),
            c2 = obj(
                c20 = obj(
                    c200 = 2,
                    ),
                c21 = obj(
                    c210 = obj(
                        c2100 = 3,
                        ),
                    ),
                ),
            c3 = obj(
                c30 = obj(
                    c300 = obj(
                        c3000 = obj(
                            c30000 = 4,
                            ),
                        ),
                    c301 = obj(
                        c3010 = obj(
                            c30100 = obj(
                                c301000 = 5,
                                ),
                            ),
                        ),
                    ),
                c31 = obj(
                    c310 = obj(
                        c3100 = obj(
                            c31000 = obj(
                                c310000 = obj(
                                    c3100000 = 6,
                                    ),
                                ),
                            ),
                        ),
                    c311 = obj(
                        c3110 = obj(
                            c31100 = obj(
                                c311000 = obj(
                                    c3110000 = obj(
                                        c3110000 = 7,
                                        )
                                    ),
                                ),
                            ),
                        ),
                    ),
                ),
            ),
        )

    oref = obj({
        'a/a0':0,
        'a/a1/a10':1,
        'b/b0':0,
        'b/b1/b10':1,
        'b/b2/b20/b200':2,
        'b/b2/b21/b210/b2100':3,
        'c/c0':0,
        'c/c1/c10':1,
        'c/c2/c20/c200':2,
        'c/c2/c21/c210/c2100':3,
        'c/c3/c30/c300/c3000/c30000':4,
        'c/c3/c30/c301/c3010/c30100/c301000':5,
        'c/c3/c31/c310/c3100/c31000/c310000/c3100000':6,
        'c/c3/c31/c311/c3110/c31100/c311000/c3110000/c3110000':7,
        })

    for k,v in oref.items():
        assert(o.get_path(k)==v)
    #end for
    o2 = o.serial()
    assert(object_eq(o2,oref))

#end def test_extensions

