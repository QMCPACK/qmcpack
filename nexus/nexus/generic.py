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


def log(
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
    log('\n  '+header,logfile=logfile)
    log(msg.rstrip(),indent=indent,logfile=logfile)
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


class Void(object):
    void_items: ClassVar[dict] = dict()

    @classmethod
    def _unavailable(cls,self):
        sid = id(self)
        if sid in Void.void_items:
            module,item = Void.void_items[id(self)]
        else:
            module,item = None,None
        #end if
        if module is None and item is None:
            msg = 'encountered a void item from an unavailable module'
        elif module is None:
            msg = 'item '+str(item)+' is from an unavailable module'
        elif item is None:
            msg = (
                'encountered a void item from unavailable module '+str(module)+'  \n'
                'this python module must be installed on your system to use this feature'
                )
        else:
            msg = (
                'item '+str(item)+' is from unavailable module '+str(module)+'  \n'
                'this python module must be installed on your system to use this feature'
                )
        #end if
        raise ImportError(msg)
    #end def _unavailable


    @classmethod
    def _class_unavailable(cls):
        msg = 'encountered a void item from an unavailable module'
        raise ImportError(msg)
    #end def _class_unavailable


    def __init__(self,module=None,item=None):
        Void.void_items[id(self)] = module,item
    #end def __init__


    #list of magic functions taken from the following sources
    #  http://web.archive.org/web/20110131211638/http://diveintopython3.org/special-method-names.html
    #  http://www.ironpythoninaction.com/magic-methods.html


    #class methods
    @classmethod
    def __instancecheck__(cls,*args,**kwargs):       Void._class_unavailable()
    @classmethod
    def __subclasscheck__(cls,*args,**kwargs):       Void._class_unavailable()
    @classmethod
    def __subclasshook__(cls,*args,**kwargs):        Void._class_unavailable()
    

    #member methods
    #def __new__(self,*args,**kwargs):
    #    Void._unavailable(self)
    def __eq__(self,*args,**kwargs):           Void._unavailable(self)
    def __ne__(self,*args,**kwargs):           Void._unavailable(self)
    def __lt__(self,*args,**kwargs):           Void._unavailable(self)
    def __le__(self,*args,**kwargs):           Void._unavailable(self)
    def __gt__(self,*args,**kwargs):           Void._unavailable(self)
    def __ge__(self,*args,**kwargs):           Void._unavailable(self)
    def __nonzero__(self,*args,**kwargs):      Void._unavailable(self)
    def __subclasses__(self,*args,**kwargs):   Void._unavailable(self)
    def __call__(self,*args,**kwargs):         Void._unavailable(self)
    def __hash__(self,*args,**kwargs):         Void._unavailable(self)
    #def __del__(self,*args,**kwargs):         Void._unavailable(self)
    def __dir__(self,*args,**kwargs):          Void._unavailable(self)
    def __getitem__(self,*args,**kwargs):      Void._unavailable(self)
    def __setitem__(self,*args,**kwargs):      Void._unavailable(self)
    def __delitem__(self,*args,**kwargs):      Void._unavailable(self)
    def __len__(self,*args,**kwargs):          Void._unavailable(self)
    def __contains__(self,*args,**kwargs):     Void._unavailable(self)
    def __iter__(self,*args,**kwargs):         Void._unavailable(self)
    def __reversed__(self,*args,**kwargs):     Void._unavailable(self)
    def __missing__(self,*args,**kwargs):      Void._unavailable(self)
    def __length_hint__(self,*args,**kwargs):  Void._unavailable(self)
    def __repr__(self,*args,**kwargs):         Void._unavailable(self)
    def __str__(self,*args,**kwargs):          Void._unavailable(self)
    def __unicode__(self,*args,**kwargs):      Void._unavailable(self)
    def __getattr__(self,*args,**kwargs):      Void._unavailable(self)
    def __setattr__(self,*args,**kwargs):      Void._unavailable(self)
    def __delattr__(self,*args,**kwargs):      Void._unavailable(self)
    def __getattribute__(self,*args,**kwargs): Void._unavailable(self)
    def __add__(self,*args,**kwargs):          Void._unavailable(self)
    def __sub__(self,*args,**kwargs):          Void._unavailable(self)
    def __mul__(self,*args,**kwargs):          Void._unavailable(self)
    def __floordiv__(self,*args,**kwargs):     Void._unavailable(self)
    def __div__(self,*args,**kwargs):          Void._unavailable(self)
    def __truediv__(self,*args,**kwargs):      Void._unavailable(self)
    def __mod__(self,*args,**kwargs):          Void._unavailable(self)
    def __divmod__(self,*args,**kwargs):       Void._unavailable(self)
    def __pow__(self,*args,**kwargs):          Void._unavailable(self)
    def __lshift__(self,*args,**kwargs):       Void._unavailable(self)
    def __rshift__(self,*args,**kwargs):       Void._unavailable(self)
    def __and__(self,*args,**kwargs):          Void._unavailable(self)
    def __xor__(self,*args,**kwargs):          Void._unavailable(self)
    def __or__(self,*args,**kwargs):           Void._unavailable(self)
    def __neg__(self,*args,**kwargs):          Void._unavailable(self)
    def __pos__(self,*args,**kwargs):          Void._unavailable(self)
    def __abs__(self,*args,**kwargs):          Void._unavailable(self)
    def __invert__(self,*args,**kwargs):       Void._unavailable(self)
    def __complex__(self,*args,**kwargs):      Void._unavailable(self)
    def __int__(self,*args,**kwargs):          Void._unavailable(self)
    def __float__(self,*args,**kwargs):        Void._unavailable(self)
    def __oct__(self,*args,**kwargs):          Void._unavailable(self)
    def __hex__(self,*args,**kwargs):          Void._unavailable(self)
    def __index__(self,*args,**kwargs):        Void._unavailable(self)
    def __enter__(self,*args,**kwargs):        Void._unavailable(self)
    def __exit__(self,*args,**kwargs):         Void._unavailable(self)
    def __get__(self,*args,**kwargs):          Void._unavailable(self)
    def __set__(self,*args,**kwargs):          Void._unavailable(self)
    def __delete__(self,*args,**kwargs):       Void._unavailable(self)
    def __doc__(self,*args,**kwargs):          Void._unavailable(self)
    def __dict__(self,*args,**kwargs):         Void._unavailable(self)
    #def __slots__(self,*args,**kwargs):       Void._unavailable(self)
    def __class__(self,*args,**kwargs):        Void._unavailable(self)
    def __bases__(self,*args,**kwargs):        Void._unavailable(self)
    def __name__(self,*args,**kwargs):         Void._unavailable(self)
    def __all__(self,*args,**kwargs):          Void._unavailable(self)
    def __file__(self,*args,**kwargs):         Void._unavailable(self)
    def __module__(self,*args,**kwargs):       Void._unavailable(self)
    #def __metaclass__(self,*args,**kwargs):   Void._unavailable(self)
    def __import__(self,*args,**kwargs):       Void._unavailable(self)
    def __radd__(self,*args,**kwargs):         Void._unavailable(self)
    def __rsub__(self,*args,**kwargs):         Void._unavailable(self)
    def __rmul__(self,*args,**kwargs):         Void._unavailable(self)
    def __rtruediv__(self,*args,**kwargs):     Void._unavailable(self)
    def __rfloordiv__(self,*args,**kwargs):    Void._unavailable(self)
    def __rmod__(self,*args,**kwargs):         Void._unavailable(self)
    def __rdivmod__(self,*args,**kwargs):      Void._unavailable(self)
    def __rpow__(self,*args,**kwargs):         Void._unavailable(self)
    def __rlshift__(self,*args,**kwargs):      Void._unavailable(self)
    def __rrshift__(self,*args,**kwargs):      Void._unavailable(self)
    def __rand__(self,*args,**kwargs):         Void._unavailable(self)
    def __rxor__(self,*args,**kwargs):         Void._unavailable(self)
    def __ror__(self,*args,**kwargs):          Void._unavailable(self)
    def __iadd__(self,*args,**kwargs):         Void._unavailable(self)
    def __isub__(self,*args,**kwargs):         Void._unavailable(self)
    def __imul__(self,*args,**kwargs):         Void._unavailable(self)
    def __itruediv__(self,*args,**kwargs):     Void._unavailable(self)
    def __ifloordiv__(self,*args,**kwargs):    Void._unavailable(self)
    def __imod__(self,*args,**kwargs):         Void._unavailable(self)
    def __ipow__(self,*args,**kwargs):         Void._unavailable(self)
    def __ilshift__(self,*args,**kwargs):      Void._unavailable(self)
    def __irshift__(self,*args,**kwargs):      Void._unavailable(self)
    def __iand__(self,*args,**kwargs):         Void._unavailable(self)
    def __ixor__(self,*args,**kwargs):         Void._unavailable(self)
    def __ior__(self,*args,**kwargs):          Void._unavailable(self)
    def __round__(self,*args,**kwargs):        Void._unavailable(self)
    def __ceil__(self,*args,**kwargs):         Void._unavailable(self)
    def __floor__(self,*args,**kwargs):        Void._unavailable(self)
    def __trunc__(self,*args,**kwargs):        Void._unavailable(self)
    def __bool__(self,*args,**kwargs):         Void._unavailable(self)
    def __copy__(self,*args,**kwargs):         Void._unavailable(self)
    def __deepcopy__(self,*args,**kwargs):     Void._unavailable(self)
    def __getstate__(self,*args,**kwargs):     Void._unavailable(self)
    def __reduce__(self,*args,**kwargs):       Void._unavailable(self)
    def __reduce_ex__(self,*args,**kwargs):    Void._unavailable(self)
    def __getnewargs__(self,*args,**kwargs):   Void._unavailable(self)
    def __setstate__(self,*args,**kwargs):     Void._unavailable(self)
    def __bytes__(self,*args,**kwargs):        Void._unavailable(self)
    def __format__(self,*args,**kwargs):       Void._unavailable(self)
    def __next__(self,*args,**kwargs):         Void._unavailable(self)
#end class Void


def unavailable(module,*items):
    voids = []
    if len(items)==0:
        voids.append(Void(module))
    #end if
    for item in items:
        voids.append(Void(module,item))
    #end for
    if len(voids)==1:
        return voids[0]
    else:
        return voids
    #end if
#end def unavailable


def available(*items):
    for item in items:
        if isinstance(item,Void):
            return False
        #end if
    #end for
    return True
#end def available
