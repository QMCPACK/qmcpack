"""Docstrings for :mod:`nexus.utilities`."""


to_str = """Convert a value to a string"""


_path_to_str = """Simple conversion from bytes/Path types to str."""


is_valid_path = """Screen out paths with invalid characters."""


is_valid_filename = """Screen out filenames with invalid characters."""


is_relative_path = """Determine if a path is relative to some current working directory."""


path_string = """Convert a path to a string.

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
