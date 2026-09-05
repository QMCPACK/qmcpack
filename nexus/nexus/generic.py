##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  generic.py                                                        #
#    Base class for all Nexus classes (obj).  Support for hidden     #
#    data UI (hidden).                                               #
#                                                                    #
#  Content summary:                                                  # 
#    obj                                                             #
#      Base class for all Nexus classes.                             #
#      Inherits from AllAbilities and wraps all functions for UI.    #
#      Also basic working object/class for generic use.              #
#      Can function like a standard dict, also mixes in parts of     #
#        the list interface.                                         #
#                                                                    #
#====================================================================#


import functools
import os
import sys
import traceback
import warnings
from typing import NoReturn, TextIO, ClassVar, TypeAlias

VersionStr: TypeAlias = str

class generic_settings:
    devlog         = sys.stdout

    # Warnings for trying to reference or set `raise_error`
    @property
    def raise_error(self):
        warn("Referencing `raise_error` is deprecated!", warn_type="dev")
    @raise_error.setter
    def raise_error(self, _):
        warn("Setting `raise_error` is deprecated!", warn_type="dev")
#end class generic_settings


class NexusError(Exception):
    """Exception for errors that are caused by a bug in Nexus."""
#end class NexusError


class FileFormatError(Exception):
    """Exception raised when a file is not formatted as expected."""
#end class FileFormatError


class NotAnElementError(ValueError):
    """Exception raised when a value can not be interpreted as an element."""
#end class FileFormatError


class NexusDevWarning(Warning):
    def __init__(self, msg: str, indent: str = "    ", cls: str | None = None):
        self.message = msg
        self.indent  = indent
        self.cls     = cls
#end class NexusDevWarning


class NexusUserWarning(NexusDevWarning):
    pass
#end class NexusUserWarning


# Hook for replacing `warnings.showwarning`
def __nexus_showwarning(message, category, filename, lineno, file=None, line=None):
    if file is None:
        file = generic_settings.devlog

    indent = ""
    cls    = ""
    if isinstance(message, DeprecationWarning):
        message = str(message)
        # To pass indentation through `DeprecationWarning` we prepend it
        # to the message along with a pipe character "|".
        # Here we strip the whitespace indentation, then make sure the
        # message starts with the pipe character before splitting it to
        # retrieve the indentation.
        # This reduces the chance of an accidental split.
        if message.lstrip().startswith("|"):
            indent, msg = message.split("|", maxsplit=1)
        else:
            msg = message
    elif isinstance(message, NexusDevWarning | NexusUserWarning):
        msg    = message.message
        indent = message.indent
        cls    = f"{message.cls}:" if message.cls is not None else ""
    else:
        msg = str(message)

    if os.path.exists(filename):
        filename = os.path.realpath(filename)

    msg = (indent*2)+msg.strip().replace("\n", "\n"+(indent*2))

    file.write(
        f"{indent}{filename}:{lineno}:{cls} {category.__name__}:\n"
        f"{msg}\n"
        )

warnings.showwarning = __nexus_showwarning


exit_call = sys.exit


def nocopy(value):
    return value
#end def nocopy


def nxs_print(
    *items,
    indent: str | None = None,
    logfile: TextIO | None = None,
    n: int = 0
    ) -> None:
    logfile = logfile if logfile is not None else generic_settings.devlog
    if n!=0:
        if indent is None:
            indent = n*'  '
        else:
            indent = n*indent

    if len(items)==1 and isinstance(items[0],str):
        s = items[0]
    else:
        s=''
        for item in items:
            s+=str(item)+' '
        #end for
    #end if
    if len(s)>0:
        if isinstance(indent,str):
            s=indent+s.replace('\n','\n'+indent)
        #end if
        s += '\n'
    #end if
    logfile.write(s)
#end def log


def message(msg,header=None,post_header=' message:',indent='    ',logfile=None):
    if logfile is None:
        logfile = generic_settings.devlog
    #end if
    if header is None:
        header = post_header.lstrip()
    else:
        header += post_header
    #end if
    nxs_print('\n  '+header,logfile=logfile)
    nxs_print(msg.rstrip(),indent=indent,logfile=logfile)
#end def message


def warn(
    msg      : str,
    indent   : str        = "    ",
    warn_type: str        = "user",
    cls      : str | None = None,
    ) -> None:
    """Report a warning.

    Parameters
    ----------
    msg : str
        The message for the warning.
    indent : str. default="    "
        The indentation for the warning message.
    warn_type : {"user", "dev", "deprecate", "class"}
        See the notes for a description of when you should use each
        warning type.
    cls : str, optional
        Meant to be passed through ``self.warn()`` in a class that
        derives from ``object_interface``.

    Notes
    -----
    If you are a developer who wants to use this function but are unsure
    of what kind of warning type to use, here are some basic guidelines
    for choosing:

    - If the function you're putting a warning in will be directly
      called by a user, you should use ``"user"``, as this will show
      them where in their own script they've called the function.
    - If the function you're putting a warning in will be indirectly
      called by a function that a user has called, **and** if the
      problem is due to something the user has done you should use the
      class's ``self.warn()`` functionality if it inherits from
      ``DevBase``, otherwise use the warning type ``"class"`` and pass
      in the name of the class and/or function you're calling from.
    - If the function is called indirectly by a user and it is not
      necessarily their fault that something has gone awry, you should
      use ``"dev"``. This should likely only be used in places where
      you are *mostly* sure there won't be any problems but the code
      has not been thoroughly tested.
    - If you are deprecating a function, use the ``@nxs_deprecate``
      decorator and supply the version and the replacement that a user
      should use instead of the soon-to-be-removed function.
    """
    match warn_type:
        case "dev":
            msg = NexusDevWarning(msg, indent, cls)
            warnings.warn(msg, stacklevel=2)
        case "user":
            msg = NexusUserWarning(msg, indent, cls)
            # len(traceback.format_stack()) brings you all the way back
            # to the user file no matter where the warning was raised.
            warnings.warn(msg, stacklevel=len(traceback.format_stack()))
        case "deprecate":
            msg = DeprecationWarning(f"{indent}|{msg}")
            warnings.warn(msg, stacklevel=3)
        case "class":
            msg = NexusUserWarning(msg, indent, cls)
            warnings.warn(msg, stacklevel=3)
#end def warn


def nxs_deprecate(
    since: VersionStr,
    replacement: str,
    indent: str = "    "
    ):
    """Decorator for signaling the deprecation of a Nexus function.

    This should itself be deprecated when the minimum Python version is
    3.13, when the official ``warnings.deprecated()`` decorator is added.

    Parameters
    ----------
    since : Version
        The version of Nexus when the function was deprecated.
    replacement : str
        The replacement function for the function being deprecated.
    indent : str, default = "    "
        Indentation to go before the warning.
    """
    def inner(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            warn_msg = f"{f.__qualname__} is deprecated as of Nexus version {since}, and will be removed in a future update."
            if replacement is not None:
                warn_msg += f"\nPlease use {replacement} instead."

            warn(msg=warn_msg, indent=indent, warn_type="deprecate")
            return f(*args, **kwargs)
        return wrapper
    return inner
#end def nxs_deprecate


def error(msg: str, header: str | None = None) -> NoReturn:
    """Raise a ``NexusError``

    Parameters
    ----------
    msg : str
        The message to write
    header : str, optional
        Header to print before the message. Will have ``" error:"`` appended to
        it.
    """
    header = f"{header} error:" if header is not None else "error:"
    msg = f"{header}\n{msg}"
    raise NexusError(msg)
#end def error
