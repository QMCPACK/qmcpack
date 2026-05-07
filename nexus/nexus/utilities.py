
# attempt to regain python 2 sorting
# code below is from https://stackoverflow.com/questions/26575183/how-can-i-get-2-x-like-sorting-behaviour-in-python-3-x
#===========================
from numbers import Number


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
    _unhandled_types = {complex}

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




def to_str(s):
    '''Convert a value to a string'''
    if isinstance(s,bytes):
        return str(s,encoding='utf-8')
    else:
        return str(s)
    #end if
#end def to_str



import string
import os
from pathlib import Path



def is_valid_path(path):
    '''Screen out paths with invalid characters'''
    if isinstance(path,Path):
        path = str(path)
    if not isinstance(path,str):
        raise ValueError('path must be of type Path or str')
    if not hasattr(is_valid_path,'invalid_chars'):
        unprintable = [chr(c) for c in range(128) if chr(c) not in string.printable]
        special = '!@#$%^&*;|?\`",()[]{}<>' + "'"
        invalid = set(unprintable) | set(special) | set(string.whitespace)
        is_valid_path.invalid_chars = invalid
    invalid_chars = is_valid_path.invalid_chars
    is_valid = len( set(path) & invalid_chars )==0
    return is_valid
#end def is_valid_path



def is_relative_path(path):
    '''Determine if a path is relative to some current working directory'''
    if isinstance(path,Path):
        path = str(path)
    if not isinstance(path,str):
        raise ValueError('path must be of type Path or str')
    absolute = len(path)>0 and (path[0]=='/' or path[0]=='~')
    relative = not absolute
    return relative
#end def is_relative_path



def path_string(path,check=True):
    """Convert a path to a string.

    Parameters
    ----------
    path : str or Path
        A file path or directory path. 
    check: bool, default=True
      Check if a path contains only valid characters.
      ValueError is raised for invalid paths.

    Returns
    -------
    path_out : str
        The path as a string.
    """
    if isinstance(path,str):
        path_out = path
    elif isinstance(path,Path):
        path_out = str(path)
    else:
        raise ValueError('path must be of type "str" or "Path"')
    if check and not is_valid_path(path_out):
        raise ValueError('path contains invalid characters:\n'+path_out)
    return path_out
#end def path_string
