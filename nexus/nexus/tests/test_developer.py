import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.DEVELOPER)

from collections.abc import Mapping, MutableMapping



from ..testing import failed,FailedTest


def test_unavailable():
    from ..developer import Void, NexusError
    from ..developer import unavailable, available

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
            with pytest.raises(
                ImportError,
                match="this python module must be installed on your system to use this feature",
                ):
                op(v)
        #end for
    #end for
#end def test_unavailable


def test_valid_variable_name():
    from ..utilities import valid_variable_name

    assert(valid_variable_name('valid_variable_name'))
    assert(valid_variable_name('_valid_variable_name'))
    assert(valid_variable_name('__valid_variable_name'))
    assert(valid_variable_name('valid_variable_name__'))
    assert(valid_variable_name('validVariableName'))
    assert(not valid_variable_name('invalid variable name'))
    assert(not valid_variable_name('invalid_variable-name'))
    assert(not valid_variable_name('inval!d_variable_name'))
    assert(not valid_variable_name('inv@lid_variable_name'))
    assert(not valid_variable_name('in>alid_variable_name'))
    assert(not valid_variable_name('invalid_variable_name '))

#end def test_valid_variable_name


def check_dictlike(dict_type,*,check_repr_str=True,check_iter=True,check_copy=True):
    """Exercise the standard constructor and method interface of *dict_type*.

    The function returns ``None`` on success and raises ``AssertionError`` when
    an operation does not have normal ``dict`` semantics.  Each destructive
    group starts with a new instance so calling this function has no lasting
    side effects on an instance supplied by the caller.
    """
    def check(condition, operation):
        if not condition:
            raise AssertionError(
                '{} does not satisfy dict semantics for {}'.format(
                    getattr(dict_type, '__name__', repr(dict_type)), operation
                    )
                )

    # Constructors: empty, mapping, iterable, keywords, and mixed arguments.
    empty = dict_type()
    check(len(empty) == 0, 'dict_type()')
    check(not empty, 'empty truth value')
    check(dict_type({'a': 1, 'b': 2}) == {'a': 1, 'b': 2},
          'mapping constructor')
    check(dict_type([('a', 1), ('b', 2)]) == {'a': 1, 'b': 2},
          'iterable constructor')
    check(dict_type(a=1, b=2) == {'a': 1, 'b': 2},
          'keyword constructor')
    check(dict_type({'a': 1, 'b': 2}, a=3, c=4) ==
          {'a': 3, 'b': 2, 'c': 4}, 'mapping and keyword constructor')

    # Text representations should match dict for empty and populated values.
    if check_repr_str:
        check(repr(empty) == repr({}), '__repr__ of empty mapping')
        check(str(empty) == str({}), '__str__ of empty mapping')
        represented = dict_type([('a', 1), ('b', 'two')])
        reference = {'a': 1, 'b': 'two'}
        check(repr(represented) == repr(reference),
              '__repr__ of populated mapping')
        check(str(represented) == str(reference),
              '__str__ of populated mapping')

    # Core mapping protocol and dynamic views.
    mapping = dict_type({'a': 1, 'b': 2})
    check(len(mapping) == 2, '__len__')
    check(mapping['a'] == 1, '__getitem__')
    check('a' in mapping and 'missing' not in mapping, '__contains__')
    if check_iter:
        check(set(iter(mapping)) == {'a', 'b'}, '__iter__')
    check(set(mapping.keys()) == {'a', 'b'}, 'keys')
    check(set(mapping.items()) == {('a', 1), ('b', 2)}, 'items')
    check(sorted(mapping.values()) == [1, 2], 'values')
    keys = mapping.keys()
    items = mapping.items()
    mapping['c'] = 3
    check(mapping['c'] == 3 and 'c' in keys and ('c', 3) in items,
          '__setitem__ and dynamic views')
    del mapping['c']
    check('c' not in mapping, '__delitem__')

    # Lookup helpers.
    check(mapping.get('a') == 1, 'get of present key')
    check(mapping.get('missing') is None, 'get default of None')
    sentinel = object()
    check(mapping.get('missing', sentinel) is sentinel, 'get explicit default')
    check(mapping.setdefault('a', 9) == 1 and mapping['a'] == 1,
          'setdefault of present key')
    check(mapping.setdefault('c', 3) == 3 and mapping['c'] == 3,
          'setdefault of missing key')

    # update accepts the same input forms as the constructor.
    mapping = dict_type()
    check(mapping.update({'a': 1}) is None, 'update return value')
    mapping.update([('b', 2)])
    mapping.update(c=3)
    mapping.update({'a': 0, 'd': 4}, a=5)
    check(mapping == {'a': 5, 'b': 2, 'c': 3, 'd': 4}, 'update inputs')

    # Removal methods and their error/default behavior.
    check(mapping.pop('d') == 4 and 'd' not in mapping, 'pop present key')
    check(mapping.pop('missing', sentinel) is sentinel, 'pop default')
    try:
        mapping.pop('missing')
    except KeyError:
        pass
    else:
        check(False, 'pop missing key')
    before = len(mapping)
    popped = mapping.popitem()
    check(isinstance(popped, tuple) and len(popped) == 2,
          'popitem return value')
    check(len(mapping) == before - 1 and popped[0] not in mapping,
          'popitem removal')
    mapping.clear()
    check(len(mapping) == 0, 'clear')
    try:
        mapping.popitem()
    except KeyError:
        pass
    else:
        check(False, 'popitem on empty mapping')

    # copy must be shallow and independent at the mapping level.
    if check_copy:
        shared = []
        original = dict_type({'a': shared})
        copied = original.copy()
        check(copied == original and copied is not original, 'copy')
        check(copied['a'] is shared, 'copy is shallow')
        copied['b'] = 2
        check('b' not in original, 'copy is independent')

    # deepcopy must preserve the mapping type, recursively copy mutable values,
    # retain shared-reference relationships, and support recursive mappings.
    from copy import deepcopy
    shared = [1, 2]
    original = dict_type({'first': shared, 'second': shared})
    original['self'] = original
    copied = deepcopy(original)
    check(type(copied) is dict_type, 'deepcopy preserves mapping type')
    check(copied is not original, 'deepcopy creates a new mapping')
    check(copied['first'] == shared and copied['first'] is not shared,
          'deepcopy recursively copies values')
    check(copied['first'] is copied['second'],
          'deepcopy preserves shared references')
    check(copied['self'] is copied, 'deepcopy preserves recursive references')
    copied['first'].append(3)
    check(original['first'] == [1, 2], 'deepcopy is independent')

    # fromkeys is callable through the class and uses one shared default.
    default = []
    created = dict_type.fromkeys(('a', 'b'), default)
    check(created == {'a': default, 'b': default}, 'fromkeys with default')
    check(created['a'] is default and created['b'] is default,
          'fromkeys preserves default identity')
    check(dict_type.fromkeys(('a', 'b')) == {'a': None, 'b': None},
          'fromkeys default of None')

    # Missing and unhashable keys have the same error categories as dict.
    try:
        original['missing']
    except KeyError:
        pass
    else:
        check(False, '__getitem__ of missing key')
    try:
        original[[]] = 1
    except TypeError:
        pass
    else:
        check(False, 'unhashable key')
#end def check_dictlike


def check_dictlike_pair(dict_type1,dict_type2):
    """Check that two dict-like classes interoperate in both directions."""
    name1 = getattr(dict_type1, '__name__', repr(dict_type1))
    name2 = getattr(dict_type2, '__name__', repr(dict_type2))

    def check(condition, operation):
        if not condition:
            raise AssertionError(
                '{} and {} are incompatible for {}'.format(
                    name1, name2, operation
                    )
                )

    contents = [('a', 1), ('b', 'two')]
    left  = dict_type1(contents)
    right = dict_type2(reversed(contents))

    # Equality must be symmetric, independent of insertion order, and return
    # an actual bool.  Inequality must also agree in both directions.
    left_to_right = left == right
    right_to_left = right == left
    check(type(left_to_right) is bool and left_to_right,
          '{} == {}'.format(name1, name2))
    check(type(right_to_left) is bool and right_to_left,
          '{} == {}'.format(name2, name1))
    check(not (left != right) and not (right != left),
          'symmetric inequality of equal mappings')

    unequal_value = dict_type2([('a', 9), ('b', 'two')])
    unequal_keys  = dict_type2([('a', 1), ('c', 'two')])
    ue_pairs = [(unequal_value, 'different values'),
                (unequal_keys , 'different keys'  )]
    for unequal, description in ue_pairs:
        check(not (left == unequal) and not (unequal == left),
              'symmetric equality with {}'.format(description))
        check(left != unequal and unequal != left,
              'symmetric inequality with {}'.format(description))

    empty1 = dict_type1()
    empty2 = dict_type2()
    check(empty1 == empty2 and empty2 == empty1,
          'symmetric equality of empty mappings')

    # Each class must accept the other as a mapping constructor input.  This
    # conversion is shallow, just like construction from a dict.
    shared = []
    source1 = dict_type1({'shared': shared, 'value': 1})
    source2 = dict_type2({'shared': shared, 'value': 1})
    converted1 = dict_type1(source2)
    converted2 = dict_type2(source1)
    check(converted1 == source2 and converted1['shared'] is shared,
          '{} constructor from {}'.format(name1, name2))
    check(converted2 == source1 and converted2['shared'] is shared,
          '{} constructor from {}'.format(name2, name1))

    # update must consume the other mapping type and preserve normal overwrite
    # behavior without mutating the source.
    updated1 = dict_type1({'old': 0, 'a': 0})
    updated2 = dict_type2({'old': 0, 'a': 0})
    source1_before = dict(source1.items())
    source2_before = dict(source2.items())
    updated1.update(source2)
    updated2.update(source1)
    expected = {'old': 0, 'a': 0, 'shared': shared, 'value': 1}
    check(updated1 == expected, '{}.update({})'.format(name1, name2))
    check(updated2 == expected, '{}.update({})'.format(name2, name1))
    check(source1 == source1_before and source2 == source2_before,
          'update leaves source mappings unchanged')

    # Views from either implementation should describe the same mapping and
    # support the standard set-like keys/items interoperability.
    check(left.keys() == right.keys() and right.keys() == left.keys(),
          'keys-view equality')
    check(left.items() == right.items() and right.items() == left.items(),
          'items-view equality')
    check(sorted(left.values(), key=repr) ==
          sorted(right.values(), key=repr), 'values contents')

    # Generic mapping consumers use keys plus __getitem__.  Unpacking checks
    # that protocol directly rather than relying on iteration behavior.
    check({**left} == {**right} == dict(contents), 'mapping unpacking')

    # Nested mappings should retain symmetric, structural equality as well.
    nested1 = dict_type1({'outer': dict_type1(contents)})
    nested2 = dict_type2({'outer': dict_type2(reversed(contents))})
    check(nested1 == nested2 and nested2 == nested1,
          'symmetric nested equality')
#end def check_dictlike_pair




def test_dictlike_individual():
    from collections import UserDict
    from ..developer_tools import dotdict,obj
    check_dictlike(dict)
    check_dictlike(UserDict)
    check_dictlike(dotdict)
    check_dictlike(
        obj,
        check_repr_str = False, # always here
        )
#end def test_dictlike_individual



def test_dictlike_pairs():
    from collections import UserDict
    from ..developer_tools import dotdict,obj
    dictlike = [dict,UserDict,dotdict,obj]
    for n,dict_type1 in enumerate(dictlike):
        for dict_type2 in dictlike[n:]:
            check_dictlike_pair(dict_type1,dict_type2)
#end def test_dictlike_pairs



def test_dotdict_unique():
    """Check dotdict behavior beyond the ordinary dictionary interface."""
    from ..developer_tools import dotdict
    def check(condition, operation):
        if not condition:
            raise AssertionError(
                'dotdict does not satisfy its unique semantics for {}'.format(
                    operation
                    )
                )

    # String keys support dot access, assignment, and deletion in both
    # directions.
    mapping = dotdict(a=1)
    check(mapping.a == mapping['a'] == 1, 'key attribute access')
    mapping.b = 2
    check(mapping['b'] == 2, 'attribute assignment creates a key')
    mapping['c'] = 3
    check(mapping.c == 3, 'item assignment creates an attribute')
    del mapping.b
    check('b' not in mapping, 'attribute deletion removes a key')

    # Entries live in the dict base, not in the instance attribute dictionary.
    check(vars(mapping) == {}, 'storage outside __dict__')

    # Class attributes take precedence over keys with the same name.  The key
    # remains available through item access and removing it does not affect the
    # inherited method.
    mapping['update'] = 'stored update value'
    check(mapping['update'] == 'stored update value',
          'method-name key item access')
    check(callable(mapping.update), 'method wins over same-named key')
    mapping.update(extra=4)
    check(mapping.extra == 4, 'method remains callable after key collision')
    del mapping.update
    check('update' not in mapping and callable(mapping.update),
          'deleting colliding key preserves method')

    # The same precedence applies to special method names.  Python's implicit
    # special-method lookup and explicit attribute lookup both reach the class.
    mapping['__len__'] = 'stored length value'
    check(mapping['__len__'] == 'stored length value',
          'special-method-name key item access')
    check(callable(mapping.__len__) and len(mapping) == 4,
          'class special method wins over key')

    # __getattr__ delegates directly to item lookup, so missing dot access has
    # the same KeyError behavior as a missing item (unlike normal attributes).
    try:
        mapping.missing
    except KeyError as error:
        check(error.args == ('missing',), 'missing attribute KeyError contents')
    else:
        check(False, 'missing attribute raises KeyError')
    try:
        del mapping.missing
    except KeyError as error:
        check(error.args == ('missing',),
              'missing attribute deletion KeyError contents')
    else:
        check(False, 'missing attribute deletion raises KeyError')

    # Non-string keys remain valid mapping entries but cannot be dot-accessed.
    nonstring = (1, 2)
    mapping[nonstring] = 'tuple key'
    check(mapping[nonstring] == 'tuple key', 'non-string key support')

    # dotdict supplies its own memo-aware deepcopy implementation.  Verify
    # recursive structure and shared-reference preservation.
    from copy import deepcopy
    shared = []
    recursive = dotdict(first=shared, second=shared)
    recursive.self = recursive
    copied = deepcopy(recursive)
    check(type(copied) is dotdict and copied is not recursive,
          'deepcopy type and identity')
    check(copied.first is copied.second and copied.first is not shared,
          'deepcopy shared references')
    check(copied.self is copied, 'deepcopy recursive reference')
#end def test_dotdict_unique



def test_obj_unique():
    """Check obj behavior beyond the ordinary dictionary interface."""
    from ..developer_tools import obj
    def check(condition, operation):
        if not condition:
            raise AssertionError(
                'obj does not satisfy its unique semantics for {}'.format(
                    operation
                    )
                )

    # obj uses its instance attribute dictionary as mapping storage, making
    # item and attribute operations two views of exactly the same state.
    mapping = obj(a=1)
    check(mapping.a == mapping['a'] == 1, 'key attribute access')
    mapping.b = 2
    check(mapping['b'] == 2, 'attribute assignment creates a key')
    mapping['c'] = 3
    check(mapping.c == 3, 'item assignment creates an attribute')
    check(vars(mapping) is mapping.__dict__, '__dict__ storage identity')
    check(vars(mapping) == {'a': 1, 'b': 2, 'c': 3},
          '__dict__ contains mapping entries')
    del mapping.b
    check('b' not in mapping, 'attribute deletion removes a key')

    # Because entries are instance attributes, a key can shadow an ordinary
    # member function.  Item syntax still reaches the value; class dispatch can
    # reach the hidden method; deleting the key reveals the method again.
    mapping['update'] = 'stored update value'
    check(mapping.update == mapping['update'] == 'stored update value',
          'key shadows same-named method')
    type(mapping).update(mapping, {'extra': 4})
    check(mapping.extra == 4, 'class dispatch reaches shadowed method')
    del mapping.update
    check('update' not in mapping and callable(mapping.update),
          'deleting key restores method access')
    mapping.update(restored=5)
    check(mapping.restored == 5, 'restored method remains callable')

    # Explicit lookup of a same-named special method is shadowed too, but
    # implicit special-method lookup occurs on the class and remains functional.
    mapping['__len__'] = 'stored length value'
    check(mapping.__len__ == mapping['__len__'] == 'stored length value',
          'key shadows explicit special-method access')
    check(len(mapping) == 5, 'implicit special-method lookup uses class')
    del mapping.__len__
    check(callable(mapping.__len__) and len(mapping) == 4,
          'deleting key restores explicit special method')

    # Missing attributes follow normal object semantics rather than translating
    # the lookup into a missing mapping key.
    try:
        mapping.missing
    except AttributeError:
        pass
    else:
        check(False, 'missing attribute raises AttributeError')
    try:
        del mapping.missing
    except AttributeError:
        pass
    else:
        check(False, 'missing attribute deletion raises AttributeError')

    # Non-string keys can exist in the underlying __dict__ through item access,
    # although Python attribute syntax cannot express them.
    nonstring = (1, 2)
    mapping[nonstring] = 'tuple key'
    check(mapping[nonstring] == 'tuple key' and nonstring in vars(mapping),
          'non-string key support in __dict__')

    # obj deliberately provides its own structured representations rather than
    # using dict syntax.
    represented = obj(answer=42)
    check('answer' in repr(represented) and 'int' in repr(represented),
          'typed repr')
    check('answer' in str(represented) and '= 42' in str(represented),
          'pretty str')
#end def test_obj_unique



def test_obj_legacy():
    """Test legacy obj behavior that remains part of the new obj interface."""
    from copy import deepcopy
    from ..developer_tools import obj

    # construction from mappings, iterables, and keyword arguments
    assert(len(obj())==0)
    assert(obj({'a':1,'b':2})=={'a':1,'b':2})
    assert(obj([('a',1),('b',2)])=={'a':1,'b':2})
    assert(obj(a=1,b=2)=={'a':1,'b':2})
    assert(obj({'a':1},a=2,b=3)=={'a':2,'b':3})

    # item and attribute access share the same storage
    o = obj(a=1,b='b',c=(1,1,1))
    o[3,4,5] = (5,6,7)
    assert(o.a==1)
    assert(o['a']==1)
    assert(o.b=='b')
    assert(o.c==(1,1,1))
    assert(o[3,4,5]==(5,6,7))

    o.d = 4
    assert(o['d']==4)
    o['e'] = 5
    assert(o.e==5)
    assert('a' in o)
    assert('missing' not in o)
    assert(len(o)==6)

    del o.d
    del o['e']
    assert('d' not in o and 'e' not in o)
    assert(len(o)==4)

    with pytest.raises(AttributeError):
        _ = o.missing
    with pytest.raises(KeyError):
        _ = o['missing']
    with pytest.raises(AttributeError):
        del o.missing
    with pytest.raises(TypeError):
        o[[]] = 1

    # dictionary views and lookup helpers
    reference = {'a':1,'b':'b','c':(1,1,1),(3,4,5):(5,6,7)}
    assert(dict(o.items())==reference)
    assert(set(o.keys())==set(reference.keys()))
    assert(set(o.values())==set(reference.values()))
    assert(o.get('c')==(1,1,1))
    assert(o.get('missing') is None)
    assert(o.get('missing',7)==7)

    # representations remain available, though their format is new
    assert(isinstance(repr(o),str))
    assert(isinstance(str(o),str))

    # deepcopy produces an independent object and preserves derived types
    class DerivedObj(obj):
        pass

    original = DerivedObj(a=obj(value=[1,2]),b={'x':[3,4]})
    copied = deepcopy(original)
    assert(isinstance(copied,DerivedObj))
    assert(copied is not original)
    assert(copied.a is not original.a)
    assert(copied.a.value is not original.a.value)
    assert(copied.b is not original.b)
    assert(copied.b['x'] is not original.b['x'])
    assert(copied.a.value==original.a.value)
    assert(copied.b==original.b)

    copied.a.value.append(3)
    copied.b['x'].append(5)
    assert(original.a.value==[1,2])
    assert(original.b=={'x':[3,4]})

    o.clear()
    assert(len(o)==0)
    assert(not o)

#end def test_obj_legacy


def test_developer_tools_devbase(tmp_path):
    """Exercise the standalone DevBase implementation in developer_tools."""
    from copy import deepcopy

    from ..developer_tools import DevBase

    class DerivedDevBase(DevBase):
        def __init__(self,*args,**kwargs):
            self.update(dict(*args,**kwargs))
        #end def __init__
    #end class DerivedDevBase

    # Attribute and item access share the instance dictionary.
    value = DerivedDevBase(a=1,b='two')
    value.c = 3
    value['d'] = 4
    assert(value.a==value['a']==1)
    assert(value.c==value['c']==3)
    assert(value.d==value['d']==4)
    assert(len(value)==4)
    assert('a' in value and 'missing' not in value)
    assert(set(value.keys())=={'a','b','c','d'})
    assert(set(value.values())=={1,'two',3,4})
    assert(dict(value.items())=={'a':1,'b':'two','c':3,'d':4})

    value.update({'a':10},e=5)
    assert(value.a==10 and value.e==5)
    del value['d']
    assert('d' not in value)
    with pytest.raises(KeyError):
        _ = value['missing']
    with pytest.raises(RuntimeError,match='does not support bare iteration'):
        iter(value)

    # Pretty-printing remains available and handles nested DevBase values.
    nested = DerivedDevBase(inner=DerivedDevBase(x=1),plain='value')
    assert(isinstance(repr(nested),str))
    assert(isinstance(str(nested),str))
    assert('inner' in repr(nested) and 'plain' in str(nested))

    # Dictionary operations retain the derived type and normal mapping
    # semantics without relying on removed convenience functions.
    shared = []
    protected = DerivedDevBase(a=1,b=shared)
    copied = deepcopy(protected)
    assert(type(copied) is DerivedDevBase and copied is not protected)
    assert(dict(copied._items())==dict(protected._items()))
    assert(copied.b is not shared and copied.b==shared)
    assert(set(protected._keys())=={'a','b'})
    assert(list(protected._values())==[1,shared])
    assert(protected.a==1)
    assert('missing' not in protected)
    if 'c' not in protected:
        protected.c = 3
    #end if
    assert(protected.c==3)
    assert(protected._update(d=4) is None and protected.d==4)
    popped = protected.d
    del protected.d
    assert(popped==4 and 'd' not in protected)
    key = next(iter(protected.keys()))
    popped = protected[key]
    del protected[key]
    assert(key not in protected and popped is not None)

    fromkeys = DerivedDevBase({key:shared for key in ('x','y')})
    assert(type(fromkeys) is DerivedDevBase)
    assert(fromkeys.x is shared and fromkeys.y is shared)
    protected._clear()
    assert(len(protected)==0)

    # Save/load replaces the receiving object's complete state.
    filepath = tmp_path/'developer_tools_devbase.p'
    saved = DevBase()
    saved.a = [1,2]
    saved.b = {'x':3}
    saved.save(filepath)
    loaded = DevBase()
    loaded.stale = True
    loaded.load(filepath)
    assert('stale' not in loaded)
    assert(loaded.a==[1,2] and loaded.b=={'x':3})
    assert(loaded.a is not saved.a and loaded.b is not saved.b)

    value.clear()
    assert(len(value)==0)

#end def test_developer_tools_devbase


def test_obj_virtual_subclass():
    from ..developer_tools import obj
    from ..developer import obj_nexus
    from ..developer import obj_deprecated

    assert(issubclass(obj, MutableMapping))
    assert(issubclass(obj, Mapping))
    assert(isinstance(obj(), MutableMapping))
    assert(isinstance(obj(), Mapping))

    assert(issubclass(obj_nexus, MutableMapping))
    assert(issubclass(obj_nexus, Mapping))
    assert(isinstance(obj_nexus(), MutableMapping))
    assert(isinstance(obj_nexus(), Mapping))

    assert(not issubclass(obj_deprecated, MutableMapping))
    assert(not issubclass(obj_deprecated, Mapping))
    assert(not isinstance(obj_deprecated(), MutableMapping))
    assert(not isinstance(obj_deprecated(), Mapping))
