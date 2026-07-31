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


import os
import sys
import traceback
from typing import TextIO, ClassVar
from copy import deepcopy
import pickle
from pickle import UnpicklingError
from random import randint
from pathlib import Path
import warnings
from typing import TypeAlias
import functools

from .utilities import path_string
from .developer_tools import sorted_py2 as sorted_generic

VersionStr: TypeAlias = str

class generic_settings:
    devlog         = sys.stdout
    raise_error    = False
#end class generic_settings


class NexusError(Exception):
    None
#end class NexusError


class NexusDevWarning(Warning):
    def __init__(self, msg: str, indent: str = "    ", cls: str | None = None):
        self.message = msg
        self.indent  = indent
        self.cls     = cls
#end class NexusDevWarning


class NexusUserWarning(NexusDevWarning):
    None
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



nexus_modules = [mod.stem for mod in Path(__file__).parent.iterdir() if mod.suffix == ".py"]

class NexusUnpickler(pickle.Unpickler):
    """This class is designed for backwards compatibility with pickles generated
    before Nexus was packaged (PR #5700, December 20, 2025). 
    It shouldn't touch anything but old Nexus pickles.
    """
    def find_class(self, module, name):
        if module in nexus_modules and "nexus." not in module:
            module = "nexus." + module
        if module == "nexus.generic":
            name = {
                "obj"     : "obj_deprecated",
                "DevBase" : "DevBaseDeprecated",
                }.get(name,name)
        return super().find_class(module, name)


def log(*items,**kwargs):
    indent  = None
    logfile = generic_settings.devlog
    if len(kwargs)>0:
        indent  = kwargs.pop('indent' ,None   )
        logfile = kwargs.pop('logfile',logfile)
        n       = kwargs.pop('n',0)
        if n!=0:
            if indent is None:
                indent = n*'  '
            else:
                indent = n*indent
            #end if
        #end if
        if len(kwargs)>0:
            valid = 'indent logfile n'.split()
            error('Invalid keyword arguments provided.\nInvalid keywords: {0}\nValid options are: {1}'.format(sorted(kwargs.keys()),valid),'log')
        #end if
    #end if
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


def error(
    msg     : str,
    header  : str | None    = None,
    exit    : bool          = True,
    trace   : int | None    = -1,
    indent  : str           = '    ',
    logfile : TextIO | None = None,
    ):
    """Report an error.

    Parameters
    ----------
    msg : str
        The message to write
    header : str, optional
        Header to print before the message. Will have ``" error:"`` appended to
        it.
    exit : bool, default=True
        Whether or not to force an exit through ``sys.exit()``. See Notes.
    trace : int or None, default=-1
        How much of the traceback to print.

        ``None`` prints the full traceback, positive integers limit to a number
        of frames relative to this function, and negative integers limit to a
        number of frames up to this function. Setting this to ``0`` will disable
        the traceback completely.

        The default value is ``-1``, which will print only the traceback up to
        the location that this function was called, but will not show the
        location of this function.
    indent : str, default="    "
        A string that will be used to indent the error message.
    logfile : TextIO, default=generic_settings.devlog
        The already-opened file to write to. Typically this is ``sys.stdout``.

    Raises
    ------
    NexusError
        If ``generic_settings.raise_error = True``.

    Notes
    -----
    If you set ``exit=True``, then this will call ``sys.exit(1)``, which will
    force an exit that can not be caught or handled.

    See Also
    --------
    generic_settings
    """
    if generic_settings.raise_error:
        raise NexusError(msg)
    #end if
    if logfile is None:
        logfile = generic_settings.devlog
    #end if
    post_header=' error:'
    message(msg,header,post_header,indent,logfile)
    if exit:
        log('  exiting.\n')
        if trace is True: # Preserve old behavior
            trace = None
        traceback.print_stack(limit=trace)
        #end if
        exit_call()
    #end if
#end def error



class object_interface(object):
    _logfile = sys.stdout

    def __len__(self):
        return len(self.__dict__)
    #end def __len__

    def __contains__(self,name):
        return name in self.__dict__
    #end def

    def __getitem__(self,name):
        return self.__dict__[name]
    #end def __getitem__

    def __setitem__(self,name,value):
        self.__dict__[name]=value
    #end def __setitem__

    def __delitem__(self,name):
        del self.__dict__[name]
    #end def __delitem__

    def __iter__(self):
        for item in self.__dict__:
            yield self.__dict__[item]
        #end for
    #end def __iter__

    def __repr__(self):
        s=''
        for k in sorted_generic(self._keys()):
            if not isinstance(k,str) or k[0]!='_':
                v=self.__dict__[k]
                if hasattr(v,'__class__'):
                    s+='  {0:<20}  {1:<20}\n'.format(str(k),v.__class__.__name__)
                else:
                    s+='  {0:<20}  {1:<20}\n'.format(str(k),type(v))
                #end if
            #end if
        #end for
        return s
    #end def __repr__

    def __str__(self,nindent=1):
        pad = '  '
        npad = nindent*pad
        s=''
        normal = []
        qable  = []
        for k,v in self._items():
            if not isinstance(k,str) or k[0]!='_':
                if isinstance(v,object_interface):
                    qable.append(k)
                else:
                    normal.append(k)
                #end if
            #end if
        #end for
        normal = sorted_generic(normal)
        qable  = sorted_generic(qable)
        indent = npad+18*' '
        for k in normal:
            v = self[k]
            vstr = str(v).replace('\n','\n'+indent)
            s+=npad+'{0:<15} = '.format(str(k))+vstr+'\n'
        #end for
        for k in qable:
            v = self[k]
            s+=npad+str(k)+'\n'
            s+=v.__str__(nindent+1)
            if isinstance(k,str):
                s+=npad+'end '+k+'\n'
            #end if
        #end for
        return s
    #end def __str__

    def __eq__(self,other): 
        if not hasattr(other,'__dict__'):
            return False
        #end if
        eq = True
        for sname in self.__dict__:
            if sname not in other.__dict__:
                return False
            #end if
            svar  =  self.__dict__[sname]
            ovar  = other.__dict__[sname]
            stype = type(svar)
            otype = type(ovar)
            if stype!=otype:
                return False
            #end if
            eqval = svar == ovar
            if isinstance(eqval, bool):
                eq &= eqval
            else:
                try: # accommodate numpy arrays implicitly
                    eq &= eqval.all()
                except:
                    return False
                #end try
            #end if
        #end for
        return eq
    #end def __eq__

    def tree(self,depth=None,all=False,types=False,nindent=1):
        if depth==nindent-1:
            return ''
        #end if
        pad = '  '
        npad = nindent*pad
        s=''
        normal = []
        qable  = []
        for k,v in self._items():
            if not isinstance(k,str) or k[0]!='_':
                if isinstance(v,object_interface):
                    qable.append(k)
                else:
                    normal.append(k)
                #end if
            #end if
        #end for
        normal.sort()
        qable.sort()
        indent = npad+18*' '
        if all:
            for k in normal:
                v = self[k]
                if types:
                    s+=npad+'{0:<15} = '.format(k)
                    if hasattr(v,'__class__'):
                        s+='{0:<20}'.format(v.__class__.__name__)
                    else:
                        s+='{0:<20}'.format(type(v))
                    #end if
                else:
                    s+=npad+str(k)
                #end if
                s+='\n'
            #end for
        #end if
        if all and depth!=nindent:
            for k in qable:
                v = self[k]
                s+=npad+str(k)+'\n'
                s+=v.tree(depth,all,types,nindent+1)
                if isinstance(k,str):
                    s+=npad+'end '+k+'\n'
                #end if
            #end for
        else:
            for k in qable:
                v = self[k]
                if types:
                    s+=npad+'{0:<15} = '.format(k)
                    if hasattr(v,'__class__'):
                        s+='{0:<20}'.format(v.__class__.__name__)
                    else:
                        s+='{0:<20}'.format(type(v))
                    #end if
                else:
                    s+=npad+str(k)
                #end if
                s+='\n'
                s+=v.tree(depth,all,types,nindent+1)
            #end for
        #end if
        return s
    #end def tree

    def data_repr(self,nindent=1,ret_str_keys=False):
        pad = '    '
        npad = nindent*pad
        normal = []
        qable  = []
        str_keys = True
        for k,v in self._items():
            k_str = isinstance(k,str)
            str_keys &= k_str
            if not k_str or k[0]!='_':
                if isinstance(v,object_interface):
                    qable.append(k)
                else:
                    normal.append(k)
                #end if
            #end if
        #end for
        normal = sorted_generic(normal)
        qable  = sorted_generic(qable)
        if str_keys:
            nkmax = 0
            for k in normal:
                nkmax = max(nkmax,len(k))
            #end for
            for k in qable:
                nkmax = max(nkmax,len(k))
            #end for
            k_fmt   = '{0:<'+str(nkmax)+'} = '
            k_delim = '='
            k_func  = str
        else:
            nkmax   = 20
            k_fmt   = '{0:<20} : '
            k_delim = ':'
            k_func  = repr
            o_delim = ''
        #end if
        print(str_keys,list(self.keys()))
        indent = npad+(nkmax+3)*' '
        if nindent==1:
            if str_keys:
                s = 'd = obj(\n'
            else:
                s = 'd = obj({\n'
            #end if
        else:
            s=''
        #end if
        for k in normal:
            v = self[k]
            vstr = (repr(v)+',').replace('\n','\n'+indent)
            s+=npad+k_fmt.format(k_func(k))+vstr+'\n'
        #end for
        for k in qable:
            v = self[k]
            sv,contains_str_keys = v.data_repr(nindent+1,ret_str_keys=True)
            if contains_str_keys:
                o_open  = ''
                o_close = ''
            else:
                o_open  = '{'
                o_close = '}'
            #end if
            s+=npad+k_func(k)+' {} obj({}\n'.format(k_delim,o_open)
            s+=sv
            s+=npad+pad+'{}),\n'.format(o_close)
        #end for
        if nindent==1:
            if str_keys:
                s += pad + ')\n'
            else:
                s += pad + '})\n'
            #end if
        #end if
        if not ret_str_keys:
            return s
        else:
            return s,str_keys
        #end if
    #end def data_repr


    # dict interface
    def keys(self):
        return self.__dict__.keys()
    #end def keys

    def values(self):
        return self.__dict__.values()
    #end def values

    def items(self):
        return self.__dict__.items()
    #end def items

    def copy(self):
        return deepcopy(self)
    #end def copy

    def clear(self):
        self.__dict__.clear()
    #end def clear


    # save/load
    def save(self,fpath=None):
        if fpath is None:
            fpath='./'+self.__class__.__name__+'.p'
        #end if
        fobj = open(fpath,'wb')
        binary = pickle.HIGHEST_PROTOCOL
        pickle.dump(self,fobj,binary)
        fobj.close()
        del fobj
        del binary
        return
    #end def save

    def load(self,fpath=None):
        if fpath is None:
            fpath='./'+self.__class__.__name__+'.p'
        #end if
        fobj = open(fpath,'rb')

        try:
            tmp = pickle.load(fobj)
        except (ImportError,ModuleNotFoundError):
            fobj.seek(0)
            try:
                # Old pickles from before Nexus was packaged (PR #5700, December 20 2025)
                # won't have the correct module path. The custom unpickler will handle this by 
                # prepending "nexus." to the module path
                tmp = NexusUnpickler(fobj).load()
            except UnpicklingError:
                try:
                    # NumPy pickles can use latin1 encoding
                    # They will likely still fail from an underflow since they are not pickle-compliant
                    tmp = NexusUnpickler(fobj).load(encoding='latin1')
                except UnpicklingError:
                    # fallback for files created with protocol 5
                    # in environments that only support up to protocol 4
                    try:
                        import pickle5
                        tmp = pickle5.load(fobj)
                    except ImportError:
                        have_pickle5 = False
                        error("Highest pickle protocol in current python version is {}, but {} is written using a higher protocol. Install pickle5, e.g. via pip, to enable protocol 5 in python <= 3.7.x".format(pickle.HIGHEST_PROTOCOL, fpath))
                #end try
            #end try
        #end try
        fobj.close()
        d = self.__dict__
        d.clear()
        for k,v in tmp.__dict__.items():
            d[k] = v
        #end for
        del fobj
        del tmp
        return
    #end def load


    # log, warning, and error messages
    def open_log(self,filepath):
        self._logfile = open(filepath,'w')
    #end def open_log

    def close_log(self):
        self._logfile.close()
    #end def close_log

    def write(self,s):
        self._logfile.write(s)
    #end def write

    def log(self,*items,**kwargs):
        if 'logfile' not in kwargs:
            kwargs['logfile'] = self._logfile
        #end if
        log(*items,**kwargs)
    #end def log

    def warn(self, msg: str, indent: str = "    "):
        """Warning from inside a Nexus class.

        See ``generic.warn()`` for description.
        """
        warn(msg, indent, warn_type="class", cls=type(self).__qualname__)
    #end def warn

    def error(self,message,header=None,exit=True,trace=-2):
        """Report an error inside a class.

        This is the same as ``generic.error()``, but with ``trace=-2`` as the
        default. This has the benefit of only printing a traceback up to the
        call location, and not inadvertently pointing someone to ``generic.py``.
        """
        if header is None:
            header = self.__class__.__name__
        #end if
        error(message,header,exit,trace,logfile=self._logfile)
    #end def error

    @classmethod
    def class_log(cls,message):
        log(message,logfile=cls._logfile)
    #end def class_log

    @classmethod
    def class_warn(cls,message,header=None,post_header=' Warning:'):
        if header is None:
            header=cls.__name__
        #end if
        warn(message,header,logfile=cls._logfile)
    #end def class_warn

    @classmethod
    def class_error(cls,message,header=None,exit=True,trace=-2,post_header=' Error:'):
        """Report an error relating to a class.

        See Also
        --------
        object_interface.error : Used inside subclasses of ``object_interface``.
        generic.error : Called when you are not reporting an error specific to a class.
        """
        if header is None:
            header = cls.__name__
        #end if
        error(message,header,exit,trace,logfile=cls._logfile)
    #end def class_error

    @classmethod
    def class_has(cls,k):
        return hasattr(cls,k)
    #end def classmethod

    @classmethod
    def class_keys(cls):
        return cls.__dict__.keys()
    #end def class_keys

    @classmethod
    def class_items(cls):
        return cls.__dict__.items()
    #end def class_items

    @classmethod
    def class_get(cls,k):
        return getattr(cls,k)
    #end def class_set

    @classmethod
    def class_set(cls,**kwargs):
        for k,v in kwargs.items():
            setattr(cls,k,v)
        #end for
    #end def class_set

    @classmethod
    def class_set_single(cls,k,v):
        setattr(cls,k,v)
    #end def class_set_single

    @classmethod
    def class_set_optional(cls,**kwargs):
        for k,v in kwargs.items():
            if not hasattr(cls,k):
                setattr(cls,k,v)
            #end if
        #end for
    #end def class_set_optional


    # access preserving functions
    #  dict interface
    def _keys(self,*args,**kwargs):
        return object_interface.keys(self,*args,**kwargs)
    def _values(self,*args,**kwargs):
        object_interface.values(self,*args,**kwargs)
    def _items(self,*args,**kwargs):         
        return object_interface.items(self,*args,**kwargs)         
    def _copy(self,*args,**kwargs):              
        return object_interface.copy(self,*args,**kwargs)
    def _clear(self,*args,**kwargs):
        object_interface.clear(self,*args,**kwargs)
    #  save/load
    def _save(self,*args,**kwargs):
        object_interface.save(self,*args,**kwargs)
    def _load(self,*args,**kwargs):
        object_interface.load(self,*args,**kwargs)
    #  log, warning, and error messages
    def _open_log(self,*args,**kwargs):
        object_interface.open_log(self,*args,**kwargs)
    def _close_log(self,*args,**kwargs):
        object_interface.close_log(self,*args,**kwargs)
    def _write(self,*args,**kwargs):
        object_interface.write(self,*args,**kwargs)
    def _log(self,*args,**kwargs):
        object_interface.log(self,*args,**kwargs)
    def _error(self,*args,**kwargs):
        object_interface.error(self,*args,**kwargs)
    def _warn(self,*args,**kwargs):
        object_interface.warn(self,*args,**kwargs)

#end class object_interface



class obj_deprecated(object_interface):

    def __init__(self,*vars,**kwargs):
        for var in vars:
            if isinstance(var,(dict,object_interface)):
                for k,v in var.items():
                    self[k] = v
                #end for
            else:
                self[var] = None
            #end if
        #end for
        for k,v in kwargs.items():
            self[k] = v
        #end for
    #end def __init__


    # list interface
    def append(self,value):
        self[len(self)] = value
    #end def append


    # return representations
    def list(self,*keys):
        nkeys = len(keys)
        if nkeys==0:
            keys = self._sorted_keys()
        elif nkeys==1 and isinstance(keys[0],(list,tuple)):
            keys = keys[0]
        #end if
        values = []
        for key in keys:
            values.append(self[key])
        #end if
        return values
    #end def list

    def list_optional(self,*keys):
        nkeys = len(keys)
        if nkeys==0:
            keys = self._sorted_keys()
        elif nkeys==1 and isinstance(keys[0],(list,tuple)):
            keys = keys[0]
        #end if
        values = []
        for key in keys:
            if key in self:
                values.append(self[key])
            else:
                values.append(None)
            #end if
        #end if
        return values
    #end def list_optional

    def tuple(self,*keys):
        return tuple(obj_deprecated.list(self,*keys))
    #end def tuple

    def dict(self,*keys):
        nkeys = len(keys)
        if nkeys==0:
            keys = self._keys()
        elif nkeys==1 and isinstance(keys[0],(list,tuple)):
            keys = keys[0]
        #end if
        d = dict()
        for k in keys:
            d[k] = self[k]
        #end for
        return d
    #end def dict

    def to_dict(self):
        d = dict()
        for k,v in self._items():
            if isinstance(v,obj_deprecated):
                d[k] = v._to_dict()
            else:
                d[k] = v
            #end if
        #end for
        return d
    #end def to_dict

    def obj(self,*keys):
        nkeys = len(keys)
        if nkeys==0:
            keys = self._keys()
        elif nkeys==1 and isinstance(keys[0],(list,tuple)):
            keys = keys[0]
        #end if
        o = obj_deprecated()
        for k in keys:
            o[k] = self[k]
        #end for
        return o
    #end def obj

    def to_obj(self):
        o = obj_deprecated()
        for k,v in self._items():
            if isinstance(v,obj_deprecated):
                o[k] = v._to_obj()
            else:
                o[k] = v
            #end if
        #end for
        return o
    #end def to_obj


    # list extensions
    def first(self):
        return self[min(self._keys())]
    #end def first

    def last(self):
        return self[max(self._keys())]
    #end def last

    def select_random(self): 
        return self[randint(0,len(self)-1)]
    #end def select_random


    # dict extensions
    def sorted_keys(self):
        return sorted_generic(self._keys())
    #end def sorted_keys


    def random_key(self):
        key = None
        nkeys = len(self)
        if nkeys>0:
            key = list(self._keys())[randint(0,nkeys-1)]
        #end if
        return key
    #end def random_key


    def set(self,*objs,**kwargs):
        for key,value in kwargs.items():
            self[key]=value
        #end for
        if len(objs)>0:
            for o in objs:
                for k,v in o.items():
                    self[k] = v
                #end for
            #end for
        #end if
        return self
    #end def set

    def set_optional(self,*objs,**kwargs):
        for key,value in kwargs.items():
            if key not in self:
                self[key]=value
            #end if
        #end for
        if len(objs)>0:
            for o in objs:
                for k,v in o.items():
                    if k not in self:
                        self[k] = v
                    #end if
                #end for
            #end for
        #end if
        return self
    #end def set_optional

    def get(self,key,value=None): # follow dict interface, no plural
        if key in self:
            value = self[key]
        #end if
        return value
    #end def get

    def get_optional(self,key,value=None):
        if key in self:
            value = self[key]
        #end if
        return value
    #end def get_optional

    def get_required(self,key):
        if key in self:
            value = self[key]
        else:
            obj_deprecated.error(self,'a required key is not present\nkey required: {0}\nkeys present: {1}'.format(key,self._sorted_keys()))
        #end if
        return value
    #end def get_required

    def delete(self,*keys):
        nkeys = len(keys)
        single = False
        if nkeys==0:
            keys = self._sorted_keys()
        elif nkeys==1 and isinstance(keys[0],(list,tuple)):
            keys = keys[0]
        elif nkeys==1:
            single = True
        #end if
        values = []
        for key in keys:
            values.append(self[key])
            del self[key]
        #end for
        if single:
            return values[0]
        else:
            return values
        #end if
    #end def delete

    def delete_optional(self,key,value=None):
        if key in self:
            value = self[key]
            del self[key]
        #end if
        return value
    #end def delete_optional

    def delete_required(self,key):
        if key in self:
            value = self[key]
            del self[key]
        else:
            obj_deprecated.error(self,'a required key is not present\nkey required: {0}\nkeys present: {1}'.format(key,self._sorted_keys()))
        #end if
        return value
    #end def delete_required

    def add(self,key,value):
        self[key] = value
    #end def add

    def add_optional(self,key,value):
        if key not in self:
            self[key] = value
        #end if
    #end def add_optional

    def transfer_from(self,other,keys=None,copy=False,overwrite=True):
        if keys is None:
            if isinstance(other,object_interface):
                keys = other._keys()
            else:
                keys = other.keys()
            #end if
        #end if
        if copy:
            copier = deepcopy
        else:
            copier = nocopy
        #end if
        if overwrite:
            for k in keys:
                self[k]=copier(other[k])
            #end for
        else:
            for k in keys:
                if k not in self:
                    self[k]=copier(other[k])
                #end if
            #end for            
        #end if
    #end def transfer_from

    def transfer_to(self,other,keys=None,copy=False,overwrite=True):
        if keys is None:
            keys = self._keys()
        #end if
        if copy:
            copier = deepcopy
        else:
            copier = nocopy
        #end if
        if overwrite:
            for k in keys:
                other[k]=copier(self[k])
            #end for
        else:
            for k in keys:
                if k not in self:
                    other[k]=copier(self[k])
                #end if
            #end for            
        #end if
    #end def transfer_to

    def move_from(self,other,keys=None,optional=False):
        if keys is None:
            if isinstance(other,object_interface):
                keys = list(other._keys())
            else:
                keys = list(other.keys())
            #end if
        #end if
        if not optional:
            for k in keys:
                self[k]=other[k]
                del other[k]
            #end for
        else:
            for k in keys:
                if k in other:
                    self[k]=other[k]
                    del other[k]
                #end if
            #end for
        #end if
    #end def move_from

    def move_to(self,other,keys=None,optional=False):
        if keys is None:
            keys = list(self._keys())
        #end if
        if not optional:
            for k in keys:
                other[k]=self[k]
                del self[k]
            #end for
        else:
            for k in keys:
                if k in self:
                    other[k]=self[k]
                    del self[k]
                #end if
            #end for
        #end if
    #end def move_to

    def move_from_optional(self,other,keys=None):
        self.move_from(other,keys,optional=True)
    #end def move_from_optional

    def move_to_optional(self,other,keys=None):
        self.move_to(other,keys,optional=True)
    #end def move_to_optional

    def copy_from(self,other,keys=None,deep=True):
        obj_deprecated.transfer_from(self,other,keys,copy=deep)
    #end def copy_from

    def copy_to(self,other,keys=None,deep=True):
        obj_deprecated.transfer_to(self,other,keys,copy=deep)
    #end def copy_to

    def extract(self,keys=None,optional=False):
        ext = obj_deprecated()
        ext.move_from(self,keys,optional=optional)
        return ext
    #end def extract

    def extract_optional(self,keys=None):
        return self.extract(keys,optional=True)
    #end def extract_optional

    def check_required(self,keys,exit=True):
        if not isinstance(keys,set):
            keys = set(keys)
        #end if
        missing = keys-set(self.keys())
        if exit and len(missing)>0:
            self._error('required keys are missing\nmissing keys: {0}'.format(sorted_generic(missing)))
        #end if
        return missing
    #end def check_required

    def check_types(self,types,optional=False,exit=True):
        kfail = None
        tfail = None
        if not optional:
            for k,t in types.items():
                if not isinstance(self[k],t):
                    kfail = k
                    tfail = t
                    break
                #end if
            #end for
        else:
            for k,t in types.items():
                if k in self and not isinstance(self[k],t):
                    kfail = k
                    tfail = t
                    break
                #end if
            #end for
        #end if
        if exit and kfail is not None:
            self._error('incorrect type encountered for key value\ntype required: {0}\ntype encountered: {1}\ninvalid key: {2}'.format(tfail.__name__,self[kfail].__class__.__name__,kfail))
        #end if
        return kfail,tfail
    #end def check_types

    def check_types_optional(self,types,exit=True):
        return self.check_types(types,exit=exit,optional=True)
    #end def check_types_optional

    def shallow_copy(self):
        new = self.__class__()
        for k,v in self._items():
            new[k] = v
        #end for
        return new
    #end def shallow_copy

    def inverse(self):
        new = self.__class__()
        for k,v in self._items():
            new[v] = k
        #end for
        return new
    #end def inverse

    def path_exists(self,path):
        o = self
        if isinstance(path, str | bytes | Path):
            path = path_string(path)
            path = path.split('/')
        #end if
        for p in path:
            if p not in o:
                return False
            #end if
            o = o[p]
        #end for
        return True
    #end def path_exists

    def set_path(self,path,value=None):
        o = self
        cls = self.__class__
        if isinstance(path, str | bytes | Path):
            path = path_string(path)
            path = path.split('/')
        #end if
        for p in path[0:-1]:
            if p not in o:
                o[p] = cls()
            #end if
            o = o[p]
        #end for
        o[path[-1]] = value
    #end def set_path

    def get_path(self,path,value=None):
        o = self
        if isinstance(path, str | bytes | Path):
            path = path_string(path)
            path = path.split('/')
        #end if
        for p in path[0:-1]:
            if p not in o:
                return value
            #end if
            o = o[p]
        #end for
        lp = path[-1]
        if lp not in o:
            return value
        else:
            return o[lp]
        #end if
    #end def get_path

    def serial(self,s=None,path=None):
        first = s is None
        if first:
            s = obj_deprecated()
            path = ''
        #end if
        for k,v in self._items():
            p = path+str(k)
            if isinstance(v,obj_deprecated):
                if len(v)==0:
                    s[p]=v
                else:
                    v._serial(s,p+'/')
                #end if
            else:
                s[p]=v
            #end if
        #end for
        if first:
            return s
        #end if
    #end def serial


    # access preserving functions
    #  list interface
    def _append(self,*args,**kwargs):
        obj_deprecated.append(self,*args,**kwargs)
    #  return representations
    def _list(self,*args,**kwargs):
        return obj_deprecated.list(self,*args,**kwargs)
    def _list_optional(self,*args,**kwargs):
        return obj_deprecated.list_optional(self,*args,**kwargs)
    def _tuple(self,*args,**kwargs):
        return obj_deprecated.tuple(self,*args,**kwargs)
    def _dict(self,*args,**kwargs):
        return obj_deprecated.dict(self,*args,**kwargs)
    def _to_dict(self,*args,**kwargs):
        return obj_deprecated.to_dict(self,*args,**kwargs)
    def _obj(self,*args,**kwargs):
        return obj_deprecated.obj(self,*args,**kwargs)
    def _to_obj(self,*args,**kwargs):
        return obj_deprecated.to_obj(self,*args,**kwargs)
    #  list extensions
    def _first(self,*args,**kwargs):
        return obj_deprecated.first(self,*args,**kwargs)
    def _last(self,*args,**kwargs):
        return obj_deprecated.last(self,*args,**kwargs)
    def _select_random(self,*args,**kwargs):
        return obj_deprecated.select_random(self,*args,**kwargs)
    #  dict extensions
    def _sorted_keys(self,*args,**kwargs):
        return obj_deprecated.sorted_keys(self,*args,**kwargs)
    def _random_key(self,*args,**kwargs):
        obj_deprecated.random_key(self,*args,**kwargs)
    def _set(self,*args,**kwargs):
        obj_deprecated.set(self,*args,**kwargs)
    def _set_optional(self,*args,**kwargs):
        obj_deprecated.set_optional(self,*args,**kwargs)
    def _get(self,*args,**kwargs):
        obj_deprecated.get(self,*args,**kwargs)
    def _get_optional(self,*args,**kwargs):
        obj_deprecated.get_optional(self,*args,**kwargs)
    def _get_required(self,*args,**kwargs):
        obj_deprecated.get_required(self,*args,**kwargs)
    def _delete(self,*args,**kwargs):
        obj_deprecated.delete(self,*args,**kwargs)
    def _delete_optional(self,*args,**kwargs):
        obj_deprecated.delete_optional(self,*args,**kwargs)
    def _delete_required(self,*args,**kwargs):
        obj_deprecated.delete_required(self,*args,**kwargs)
    def _add(self,*args,**kwargs):
        obj_deprecated.add(self,*args,**kwargs)
    def _add_optional(self,*args,**kwargs):
        obj_deprecated.add_optional(self,*args,**kwargs)
    def _transfer_from(self,*args,**kwargs):
        obj_deprecated.transfer_from(self,*args,**kwargs)
    def _transfer_to(self,*args,**kwargs):
        obj_deprecated.transfer_to(self,*args,**kwargs)
    def _move_from(self,*args,**kwargs):
        obj_deprecated.move_from(self,*args,**kwargs)
    def _move_to(self,*args,**kwargs):
        obj_deprecated.move_to(self,*args,**kwargs)
    def _move_from_optional(self,*args,**kwargs):
        obj_deprecated.move_from_optional(self,*args,**kwargs)
    def _move_to_optional(self,*args,**kwargs):
        obj_deprecated.move_to_optional(self,*args,**kwargs)
    def _copy_from(self,*args,**kwargs):
        obj_deprecated.copy_from(self,*args,**kwargs)
    def _copy_to(self,*args,**kwargs):
        obj_deprecated.copy_to(self,*args,**kwargs)
    def _extract(self,*args,**kwargs):
        obj_deprecated.extract(self,*args,**kwargs)
    def _extract_optional(self,*args,**kwargs):
        obj_deprecated.extract_optional(self,*args,**kwargs)
    def _check_required(self,*args,**kwargs):
        obj_deprecated.check_required(self,*args,**kwargs)
    def _check_types(self,*args,**kwargs):
        obj_deprecated.check_types(self,*args,**kwargs)
    def _check_types_optional(self,*args,**kwargs):
        obj_deprecated.check_types_optional(self,*args,**kwargs)
    def _shallow_copy(self,*args,**kwargs):
        obj_deprecated.shallow_copy(self,*args,**kwargs)
    def _inverse(self,*args,**kwargs):
        return obj_deprecated.inverse(self,*args,**kwargs)
    def _path_exists(self,*args,**kwargs):
        obj_deprecated.path_exists(self,*args,**kwargs)
    def _set_path(self,*args,**kwargs):
        obj_deprecated.set_path(self,*args,**kwargs)
    def _get_path(self,*args,**kwargs):
        obj_deprecated.get_path(self,*args,**kwargs)
    def _serial(self,*args,**kwargs):
        return obj_deprecated.serial(self,*args,**kwargs)


    # compatibility with the standard dict interface
    @classmethod
    def fromkeys(cls,keys,value=None):
        return cls(dict.fromkeys(keys,value))

    def pop(self,*args,**kwargs):
        return self.__dict__.pop(*args,**kwargs)

    def popitem(self,*args,**kwargs):
        return self.__dict__.popitem(*args,**kwargs)

    def setdefault(self,*args,**kwargs):
        return self.__dict__.setdefault(*args,**kwargs)

    def update(self,*args,**kwargs):
        return self.__dict__.update(*args,**kwargs)

#end class obj_deprecated





class DevBaseDeprecated(obj_deprecated):
    def not_implemented(self,name=None):
        if name is None:
            msg = 'a member function has not been implemented for class "{}"'.format(self.__class__.__name__)
        else:
            msg = 'member function "{}" has not been implemented for class "{}"'.format(name,self.__class__.__name__)
        self.error(msg,trace=True)
    #end def not_implemented
#end class DevBaseDeprecated



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
            msg = 'encountered a void item from unavailable module '+str(module)+'  \nthis python module must be installed on your system to use this feature'
        else:
            msg = 'item '+str(item)+' is from unavailable module '+str(module)+'  \nthis python module must be installed on your system to use this feature'
        #end if
        error(msg,'import')
    #end def _unavailable

        
    @classmethod
    def _class_unavailable(cls):
        msg = 'encountered a void item from an unavailable module'
        error(msg,'import')
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
