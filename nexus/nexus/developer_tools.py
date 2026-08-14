import copy
import pickle
from collections.abc import MutableMapping
from numbers import Number


# attempt to regain python 2 sorting
# code below is from https://stackoverflow.com/questions/26575183/how-can-i-get-2-x-like-sorting-behaviour-in-python-3-x
#===========================

# decorator for type to function mapping special cases
def per_type_cmp(type_):
    try:
        mapping = per_type_cmp.mapping
    except AttributeError:
        mapping = per_type_cmp.mapping = {}
    #end try
    def decorator(cmpfunc):
        mapping[type_] = cmpfunc
        return cmpfunc
    #end def decorator
    return decorator
#ned def per_type_cmp

class python2_sort_key(object):
    _unhandled_types = frozenset({complex})

    def __init__(self, ob):
       self._ob = ob
    #end def __init__

    def __lt__(self, other):
        _unhandled_types = self._unhandled_types
        self, other = self._ob, other._ob  # we don't care about the wrapper

        # default_3way_compare is used only if direct comparison failed
        try:
            return self < other
        except TypeError:
            pass
        #end try

        # hooks to implement special casing for types, dict in Py2 has
        # a dedicated __cmp__ method that is gone in Py3 for example.
        for type_, special_cmp in per_type_cmp.mapping.items():
            if isinstance(self, type_) and isinstance(other, type_):
                return special_cmp(self, other)
            #end if
        #end for

        # explicitly raise again for types that won't sort in Python 2 either
        if type(self) in _unhandled_types:
            raise TypeError('no ordering relation is defined for {}'.format(
                type(self).__name__))
        #end if
        if type(other) in _unhandled_types:
            raise TypeError('no ordering relation is defined for {}'.format(
                type(other).__name__))
        #end if

        # default_3way_compare from Python 2 as Python code
        # same type but no ordering defined, go by id
        if type(self) is type(other):
            return id(self) < id(other)
        #end if

        # None always comes first
        if self is None:
            return True
        #end if
        if other is None:
            return False
        #end if

        # Sort by typename, but numbers are sorted before other types
        self_tname = '' if isinstance(self, Number) else type(self).__name__
        other_tname = '' if isinstance(other, Number) else type(other).__name__

        if self_tname != other_tname:
            return self_tname < other_tname
        #end if

        # same typename, or both numbers, but different type objects, order
        # by the id of the type object
        return id(type(self)) < id(type(other))
    #end def __lt__
#end class python2_sort_key

@per_type_cmp(dict)
def dict_cmp(a, b, _s=object()):
    if len(a) != len(b):
        return len(a) < len(b)
    #end if
    adiff = min((k for k in a if a[k] != b.get(k, _s)), key=python2_sort_key, default=_s)
    if adiff is _s:
        # All keys in a have a matching value in b, so the dicts are equal
        return False
    #end if
    bdiff = min((k for k in b if b[k] != a.get(k, _s)), key=python2_sort_key)
    if adiff != bdiff:
        return python2_sort_key(adiff) < python2_sort_key(bdiff)
    #end if
    return python2_sort_key(a[adiff]) < python2_sort_key(b[bdiff])
#end def dict_cmp

def sorted_py2(iterable):
    return sorted(iterable,key=python2_sort_key)
#end def sorted_py2
#===========================



def save(o,filepath):
    with open(filepath,'wb') as f:
        binary = pickle.HIGHEST_PROTOCOL
        pickle.dump(o,f,binary)

def load(filepath):
    with open(filepath,'rb') as f:
        dl = pickle.load(f)
    return dl



def _pp_repr(self):
    s=''
    for k in sorted_py2(self.keys()):
        v = self.__dict__[k]
        if hasattr(v,'__class__'):
            s+='  {0:<20}  {1:<20}\n'.format(str(k),v.__class__.__name__)
        else:
            s+='  {0:<20}  {1:<20}\n'.format(str(k),type(v))
    return s


def _pp_str(self,nindent=1):
    pad = '  '
    npad = nindent*pad
    s=''
    normal = []
    qable  = []
    for k,v in self.items():
        if isinstance(v,(obj,DevBase)):
            qable.append(k)
        else:
            normal.append(k)
    normal = sorted_py2(normal)
    qable  = sorted_py2(qable)
    indent = npad+18*' '
    for k in normal:
        v = self[k]
        vstr = str(v).replace('\n','\n'+indent)
        s+=npad+'{0:<15} = '.format(str(k))+vstr+'\n'
    for k in qable:
        v = self[k]
        s+=npad+str(k)+'\n'
        s+=v.__str__(nindent+1)
        if isinstance(k,str):
            s+=npad+'end '+k+'\n'
    return s



class dotdict(dict):
    """Dictionary supporting attribute-style access to keys.

    ``dotdict`` is a subclass of :class:`dict` and supports the standard
    dictionary constructor and interface.  String keys that do not collide
    with class attributes can additionally be read, assigned, and deleted
    with the dot operator.

    Parameters
    ----------
    *args
        Positional arguments accepted by :class:`dict`.
    **kwargs
        Keyword arguments accepted by :class:`dict`.

    Notes
    -----
    Dictionary methods take precedence over keys during dot access.  For
    example, a key named ``"update"`` remains available with item access,
    while ``mapping.update`` continues to refer to :meth:`dict.update`.
    Missing keys accessed with the dot operator raise :class:`KeyError`, as
    they do with item access.

    Examples
    --------
    Keys and attributes provide interchangeable access when their names do
    not collide with dictionary attributes.

    >>> data = dotdict(answer=42)
    >>> data.answer
    42
    >>> data.units = "eV"
    >>> data["units"]
    'eV'

    Colliding keys are retained but require item access.

    >>> data["update"] = "stored value"
    >>> data["update"]
    'stored value'
    >>> callable(data.update)
    True
    >>> data.update(extra=1)
    >>> data.extra
    1
    """
    def __getattr__(self, item):
        return self[item]
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    def copy(self): return self.__class__(self)
    def __deepcopy__(self, memo):
        result = self.__class__.__new__(self.__class__)
        memo[id(self)] = result
        for key, value in self.items():
            result[copy.deepcopy(key, memo)] = copy.deepcopy(value, memo)
        return result



class obj:
    """Dictionary-like object sharing storage between keys and attributes.

    ``obj`` implements the usual mutable dictionary interface without
    inheriting from :class:`dict`.  Its entries are stored directly in the
    instance attribute dictionary, so item and dot access are two views of
    the same data.

    Parameters
    ----------
    *args
        Positional arguments accepted by :class:`dict`.
    **kwargs
        Keyword arguments accepted by :class:`dict`.

    Notes
    -----
    Unlike :class:`dotdict`, an entry whose key matches an ordinary ``obj``
    method shadows that method during dot access.  The hidden method can
    still be invoked through the class, but ``obj`` is most safely used as a
    drop-in dictionary replacement when keys do not overlap its methods.

    ``repr`` displays the type of each stored value, and ``str`` provides an
    indented, value-oriented representation.  These representations
    intentionally differ from dictionary syntax.  Missing attributes raise
    :class:`AttributeError`, while missing item access raises
    :class:`KeyError`.

    Examples
    --------
    Item and attribute operations act on the same entries.

    >>> data = obj(answer=42)
    >>> data.answer
    42
    >>> data["units"] = "eV"
    >>> data.units
    'eV'

    Keys can shadow methods because entries are also instance attributes.

    >>> data["update"] = "stored value"
    >>> data.update
    'stored value'
    >>> type(data).update(data, extra=1)
    >>> data.extra
    1

    The custom text representations emphasize values and their types.

    >>> o = obj(i=1,r=2.3,s='hello',o2=obj(b={1,2,3},c=[5,6,7]))
    >>> print(repr(o))
      i                     int
      o2                    obj
      r                     float
      s                     str
    >>> print(o)
      i               = 1
      r               = 2.3
      s               = hello
      o2
        b               = {1, 2, 3}
        c               = [5, 6, 7]
      end o2
    """
    # dict interface
    @classmethod
    def fromkeys(cls, keys, value=None):
        return cls(dict.fromkeys(keys, value))

    def __init__(self,*args,**kwargs):   self.__dict__.update(dict(*args,**kwargs))
    def items(self):              return self.__dict__.items()
    def clear(self):              return self.__dict__.clear()
    def copy(self):               return self.__class__(self.__dict__)
    def get(self,*a,**kw):        return self.__dict__.get(*a,**kw)
    def keys(self):               return self.__dict__.keys()
    def pop(self,*a,**kw):        return self.__dict__.pop(*a,**kw)
    def values(self):             return self.__dict__.values()
    def popitem(self,*a,**kw):    return self.__dict__.popitem(*a,**kw)
    def setdefault(self,*a,**kw): return self.__dict__.setdefault(*a,**kw)
    def update(self,*a,**kw):     return self.__dict__.update(*a,**kw)

    # basic functions, including dot access
    def __len__(self):               return len(self.__dict__)
    def __contains__(self,key):      return key in self.__dict__
    def __getitem__(self,key):       return self.__dict__[key]
    def __setitem__(self,key,value): self.__dict__[key]=value
    def __delitem__(self,key):       del self.__dict__[key]
    def __eq__(self,other):          return self.__dict__==other
    def __iter__(self):              return iter(self.__dict__)

    # pretty print
    __repr__ = _pp_repr
    __str__  = _pp_str
#end class obj



class DevBase:
    # similar to/same as dict
    def __len__(self):               return len(self.__dict__)
    def __contains__(self,key):      return key in self.__dict__
    def __getitem__(self,key):       return self.__dict__[key]
    def __setitem__(self,key,value): self.__dict__[key]=value
    def __delitem__(self,key):       del self.__dict__[key]
    def keys(self):                  return self.__dict__.keys()
    def values(self):                return self.__dict__.values()
    def items(self):                 return self.__dict__.items()
    def update(self,*a,**kw):        return self.__dict__.update(*a,**kw)
    def clear(self):                 return self.__dict__.clear()

    # protect against bare iteration
    def __iter__(self):
        raise RuntimeError('DevBase does not support bare iteration')

    # pretty print
    __repr__ = _pp_repr
    __str__  = _pp_str

    # protected dict interface
    def _keys(self):               return self.__dict__.keys()
    def _values(self):             return self.__dict__.values()
    def _items(self):              return self.__dict__.items()
    def _update(self,*a,**kw):     return self.__dict__.update(*a,**kw)
    def _clear(self):              return self.__dict__.clear()


    # save and load
    def save(self,filepath=None):
        if filepath is None:
            filepath='./'+self.__class__.__name__+'.p'
        save(self,filepath)

    def load(self,filepath=None):
        if filepath is None:
            filepath='./'+self.__class__.__name__+'.p'
        tmp = load(filepath)
        d = self.__dict__
        d.clear()
        for k,v in tmp.__dict__.items():
            d[k] = v

#end class DevBase



def to_obj(d):
    o = obj()
    for k,v in d.items():
        if hasattr(v,'__dict__'):
            o[k] = to_obj(v)
        else:
            o[k] = v
    return o
#end def to_obj

MutableMapping.register(obj)
