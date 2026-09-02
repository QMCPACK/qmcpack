import string
from pathlib import Path



def to_str(s):
    '''Convert a value to a string'''
    if isinstance(s,bytes):
        return str(s,encoding='utf-8')
    else:
        return str(s)
    #end if
#end def to_str


def valid_variable_name(s):
    """Check if a variable name contains invalid characters."""
    if not any([i in ('!"#$%&\'()*+,-./:;<=>?@[\\]^`{|}-\n\t ') for i in s]):
        return True
    else:
        return False
    #end if
#end def valid_variable_name


def _path_to_str(path: str | bytes | Path) -> str:
    '''Simple conversion from bytes/Path types to str.'''
    if isinstance(path, str):
        pass
    elif isinstance(path, bytes):
        path = str(path, encoding='utf-8')
    elif isinstance(path, Path):
        path = str(path)
    else:
        raise TypeError(
            f'path must be of type "str", "bytes" or "Path". Type received: {path.__class__.__name__}'
            
            )
    return path
#end def _path_to_str


def is_valid_path(path: str | bytes | Path) -> bool:
    '''Screen out paths with invalid characters.'''
    path = _path_to_str(path)
    if not hasattr(is_valid_path,'invalid_chars'):
        unprintable = [chr(c) for c in range(128) if chr(c) not in string.printable]
        special = r'!@#$%^&*;|?\`",()[]{}<>' + r"'"
        whitespace = set(string.whitespace) - set([' '])
        invalid = set(unprintable) | set(special) | whitespace
        is_valid_path.invalid_chars = invalid
    invalid_chars = is_valid_path.invalid_chars
    is_valid = len(set(path) & invalid_chars) == 0
    return is_valid
#end def is_valid_path


def is_valid_filename(filename: str | bytes | Path) -> bool:
    '''Screen out filenames with invalid characters.'''
    filename = _path_to_str(filename)
    is_valid = True
    if len(filename) == 0:
        is_valid = False
    elif filename == '.':
        is_valid = False
    elif filename.startswith('..'):
        is_valid = False
    elif '/' in filename:
        is_valid = False
    elif not is_valid_path(filename):
        is_valid = False

    return is_valid
#end def is_valid_filename


def is_relative_path(path: str | bytes | Path):
    '''Determine if a path is relative to some current working directory.'''
    path     = _path_to_str(path)
    absolute = len(path) > 0 and (path[0] == '/' or path[0] == '~')
    relative = not absolute
    return relative
#end def is_relative_path


def path_string(
    path:     str | bytes | Path,
    *,
    strict:   bool = False,
    relative: bool = False,
    check:    bool = False,
    ) -> str:
    """Convert a path to a string.

    Parameters
    ----------
    path : str, bytes or Path
        A file path or directory path. 
    strict : bool, default=False
        Require inputted path to be str type.
        Raises ValueError otherwise.
    relative : bool, default=False
        Check if path is a relative path.
        Raises ValueError otherwise.
    check : bool, default=True
        Check if a path contains only valid characters.
        ValueError is raised for invalid paths.

    Returns
    -------
    path_out : str
        The path as a string.
    """
    if relative and not is_relative_path(path):
        raise ValueError('path must be relative')
    if strict:
        if not isinstance(path, str):
            raise ValueError('path must strictly be str type')
        path_out = path
    else:
        path_out = _path_to_str(path)
    if relative and not is_relative_path(path_out):
        raise ValueError('path_string converted relative path into non-relative path')
    if check and not is_valid_path(path_out):
        raise ValueError('path contains invalid characters:\n'+path_out)
    return path_out
#end def path_string
