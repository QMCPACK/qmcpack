
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



import os
from pathlib import Path


def relative_path(path_start,path_end):
    '''Find the relative path between two paths'''
    rel_path = os.path.relpath(path_end,path_start)
    return rel_path
#end def relative_path


def is_relative_path(path, ret_bare=False):
    '''Determine if a path is relative to some current working directory'''
    if isinstance(path,Path):
        path = str(path)
    relative  = path.startswith('..') 
    relative |= path.startswith('./') 
    relative |= path=='.'
    relative |= path==''
    bare = False
    if not relative:
        bare = not path.startswith('/')
        relative |= bare
    elif path=='':
        bare = True
    if not ret_bare:
        return relative
    else:
        return relative,bare
#end def is_relative_path


def path_object(path            ,
                relative = None ,
                strict   = False,
                ret_rel  = False,
                ):
    """Convert a path to a ``pathlib.Path`` object.

    This function should be called in a pair with ``revert_path()``.

    Parameters
    ----------
    path : str or Path
        A file path or directory path. 
    relative: bool or None, default=None
        If bool, ensure the path is relative or not.
        If None, take no action.
    strict : bool, default=False
        If ``True``, require that ``path`` is a ``str``.
    ret_rel: bool, default=False
        If True, return whether the path is actually relative or not.
    is_rel : bool
        Indicates whether the path is actually relative or not.
        Returned only when ret_rel==True.

    Returns
    -------
    path : Path
        The path as a ``pathlib.Path`` object.
    """
    is_str  = isinstance(path,str)
    is_Path = isinstance(path,Path)
    if strict and not is_str:
        raise ValueError('path must be "str" type')
    elif not (is_str or is_Path):
        raise ValueError('path must be of type "str" or "Path"')
    is_rel = is_relative_path(path)
    if relative is not None:
        if relative and not is_rel:
            raise ValueError('path must be relative but it is non-relative')
        elif not relative and is_rel:
            raise ValueError('path must be non-relative but it is relative')
    if not is_Path:
        path = Path(path)
    if not ret_rel:
        return path
    else:
        return path,is_rel
#end def path_object


def path_string(path, 
                relative = None,
                ret_rel  = False,
                ):
    """Convert a path to a string.

    Parameters
    ----------
    path : str or Path
        A file path or directory path. 
    relative: bool or None, default=None
        If bool, ensure the path is relative or not.
        If None, take no action.
    ret_rel: bool, default=False
        If True, return whether the path is actually relative or not.

    Returns
    -------
    path : str
        The path as a ``str``.
    is_rel : bool
        Indicates whether the path is actually relative or not.
        Returned only when ret_rel==True.
    """
    is_str  = isinstance(path, str)
    is_Path = isinstance(path, Path)
    if not (is_str or is_Path):
        raise ValueError('path must be of type "str" or "Path"')
    is_rel,is_bare = is_relative_path(path, ret_bare=True)
    if relative is not None:
        if not isinstance(relative, bool):
            raise ValueError('relative must be "bool" type')
        if relative and not is_rel:
            raise ValueError('path must be relative but it is non-relative')
        elif not relative and is_rel:
            raise ValueError('path must be non-relative but it is relative')
    if is_Path:
        path = str(path)
    if path=='' or path=='.':
        path = '.'
        is_rel = True
    elif is_rel and is_bare:
        path = './'+path
    if not ret_rel:
        return path
    else:
        return path,is_rel
#end def path_string
