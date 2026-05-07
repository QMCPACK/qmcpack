
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



def relative_path(path_start,path_end):
    '''Find the relative path between two paths'''
    rel_path = os.path.relpath(path_end,path_start)
    return rel_path
#end def relative_path



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



def path_string(path, leading=None):
    """Convert a path to a string.

    Parameters
    ----------
    path : str or Path
        A file path or directory path. 
    leading: None, '', or './', default=None
      Only used if path is a Path object.
      Restores information lost when using Path(p) externally.

    Returns
    -------
    path_out : str
        The path as a string.
    """
    if leading not in (None,'','./'):
        raise ValueError("leading must be None, '', or './'")
    if isinstance(path, str):
        if leading is not None:
            raise ValueError('leading cannot be supplied for str paths.')
        path_out = path
    elif isinstance(path, Path):
        path_out = str(path)
        if leading is not None:
            if path_out=='.' and leading=='':
                path_out = ''
            elif leading=='./' and not path_out.startswith('/'):
                path_out = './'+path_out
    else:
        raise ValueError('path must be of type "str" or "Path"')
    return path_out
#end def path_string



def path_object(path                  ,
                strict         = False,
                guard          = True ,
                return_leading = False,
                ):
    """Convert a path to a pathlib.Path object.

    This function should be called in a pair with revert_path().

    Parameters
    ----------
    path : str or Path
        A file path or directory path. 
    strict : bool, default=False
        If True, require that path is a str.
    guard : bool, default=True
        Raise an exception if conversion to Path loses information.
    return_leading : bool, default=False
        Return information lost when converting to Path.

    Returns
    -------
    path_out : Path
        The path as a pathlib.Path object.
    leading  : None, '', './'
        Stores any information lost when converting to Path.
        Only returned if return_leading is True.
    """
    guard &= not return_leading
    is_str  = isinstance(path,str)
    is_Path = isinstance(path,Path)
    if strict and not is_str:
        raise ValueError('path must be "str" type')
    elif not (is_str or is_Path):
        raise ValueError('path must be of type "str" or "Path"')
    leading = None
    if is_Path:
        path_out = path
    else:
        path_out = Path(path)
        if path=='':
            leading = ''
        else:
            po = str(path_out)
            prepend  = path.startswith('./')
            prepend &= po!='.'
            prepend &= not po.startswith('..')
            if prepend:
                leading = './'
    if guard and leading is not None:
        raise RuntimeError('information has been lost in conversion to Path object')
    if not return_leading:
        return path_out
    else:
        return path_out, leading
#end def path_object
