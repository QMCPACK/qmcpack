##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  developer.py                                                      #
#    Defines developer environment.  Supplies base class for generic #
#    development (Nexus or beyond)                                   #
#                                                                    #
#  Content summary:                                                  #
#    log, warn, error                                                #
#      Function interface to logging and error handling.             #
#                                                                    #
#    DevBase                                                         #
#      Base class inheriting generic abilities for obj, etc.         #
#      Allows for unimplemented functions.                           #
#                                                                    #
#====================================================================#


from .developer_tools import save,load,_pp_repr,_pp_str,dotdict,obj,DevBase  # noqa: F401
from .debug import ci, interact  # noqa: F401


from .generic import NexusError, log, error, warn, message  # noqa: F401
from .generic import unavailable, available, Void  # noqa: F401
from .generic import obj_deprecated, DevBaseDeprecated  # noqa: F401


import traceback
from .generic import generic_settings




def deprecation_error():
    msg = (
        'A now-deprecated member function of obj has been called.\n'
        'Please report this issue to the Nexus developers immediately.\n'
        'To temporarily restore the deprecated implementation, uncomment the\n'
        'two lines at the bottom of:\n'
        '  {}'.format(__file__)
        )
    highlight = '='*79
    stack = ''.join(traceback.format_stack()[:-1])
    report = (
        '\n{0}\n{1}\n{0}\n'
        '{0}\nTraceback (most recent call last):\n{2}{0}\n'
        '{0}\n{1}\n{0}'.format(highlight,msg,stack)
        )
    raise RuntimeError(report)
#end def deprecation_error


class obj_defended:
    """Defend against deprecated method calls in the Nexus codebase"""
    def append(*args,**kwargs): deprecation_error()
    def list(*args,**kwargs): deprecation_error()
    def list_optional(*args,**kwargs): deprecation_error()
    def tuple(*args,**kwargs): deprecation_error()
    def dict(*args,**kwargs): deprecation_error()
    def to_dict(*args,**kwargs): deprecation_error()
    def obj(*args,**kwargs): deprecation_error()
    def to_obj(*args,**kwargs): deprecation_error()
    def first(*args,**kwargs): deprecation_error()
    def last(*args,**kwargs): deprecation_error()
    def select_random(*args,**kwargs): deprecation_error()
    def sorted_keys(*args,**kwargs): deprecation_error()
    def random_key(*args,**kwargs): deprecation_error()
    def set(*args,**kwargs): deprecation_error()
    def set_optional(*args,**kwargs): deprecation_error()
    def get(*args,**kwargs): deprecation_error()
    def get_optional(*args,**kwargs): deprecation_error()
    def get_required(*args,**kwargs): deprecation_error()
    def delete(*args,**kwargs): deprecation_error()
    def delete_optional(*args,**kwargs): deprecation_error()
    def delete_required(*args,**kwargs): deprecation_error()
    def add(*args,**kwargs): deprecation_error()
    def add_optional(*args,**kwargs): deprecation_error()
    def transfer_from(*args,**kwargs): deprecation_error()
    def transfer_to(*args,**kwargs): deprecation_error()
    def move_from(*args,**kwargs): deprecation_error()
    def move_to(*args,**kwargs): deprecation_error()
    def move_from_optional(*args,**kwargs): deprecation_error()
    def move_to_optional(*args,**kwargs): deprecation_error()
    def copy_from(*args,**kwargs): deprecation_error()
    def copy_to(*args,**kwargs): deprecation_error()
    def extract(*args,**kwargs): deprecation_error()
    def extract_optional(*args,**kwargs): deprecation_error()
    def check_required(*args,**kwargs): deprecation_error()
    def check_types(*args,**kwargs): deprecation_error()
    def check_types_optional(*args,**kwargs): deprecation_error()
    def shallow_copy(*args,**kwargs): deprecation_error()
    def inverse(*args,**kwargs): deprecation_error()
    def path_exists(*args,**kwargs): deprecation_error()
    def set_path(*args,**kwargs): deprecation_error()
    def get_path(*args,**kwargs): deprecation_error()
    def serial(*args,**kwargs): deprecation_error()
    def _append(*args,**kwargs): deprecation_error()
    def _list(*args,**kwargs): deprecation_error()
    def _list_optional(*args,**kwargs): deprecation_error()
    def _tuple(*args,**kwargs): deprecation_error()
    def _dict(*args,**kwargs): deprecation_error()
    def _to_dict(*args,**kwargs): deprecation_error()
    def _obj(*args,**kwargs): deprecation_error()
    def _to_obj(*args,**kwargs): deprecation_error()
    def _first(*args,**kwargs): deprecation_error()
    def _last(*args,**kwargs): deprecation_error()
    def _select_random(*args,**kwargs): deprecation_error()
    def _sorted_keys(*args,**kwargs): deprecation_error()
    def _random_key(*args,**kwargs): deprecation_error()
    def _set(*args,**kwargs): deprecation_error()
    def _set_optional(*args,**kwargs): deprecation_error()
    def _get(*args,**kwargs): deprecation_error()
    def _get_optional(*args,**kwargs): deprecation_error()
    def _get_required(*args,**kwargs): deprecation_error()
    def _delete(*args,**kwargs): deprecation_error()
    def _delete_optional(*args,**kwargs): deprecation_error()
    def _delete_required(*args,**kwargs): deprecation_error()
    def _add(*args,**kwargs): deprecation_error()
    def _add_optional(*args,**kwargs): deprecation_error()
    def _transfer_from(*args,**kwargs): deprecation_error()
    def _transfer_to(*args,**kwargs): deprecation_error()
    def _move_from(*args,**kwargs): deprecation_error()
    def _move_to(*args,**kwargs): deprecation_error()
    def _move_from_optional(*args,**kwargs): deprecation_error()
    def _move_to_optional(*args,**kwargs): deprecation_error()
    def _copy_from(*args,**kwargs): deprecation_error()
    def _copy_to(*args,**kwargs): deprecation_error()
    def _extract(*args,**kwargs): deprecation_error()
    def _extract_optional(*args,**kwargs): deprecation_error()
    def _check_required(*args,**kwargs): deprecation_error()
    def _check_types(*args,**kwargs): deprecation_error()
    def _check_types_optional(*args,**kwargs): deprecation_error()
    def _shallow_copy(*args,**kwargs): deprecation_error()
    def _inverse(*args,**kwargs): deprecation_error()
    def _path_exists(*args,**kwargs): deprecation_error()
    def _set_path(*args,**kwargs): deprecation_error()
    def _get_path(*args,**kwargs): deprecation_error()
    def _serial(*args,**kwargs): deprecation_error()
#end class obj_defended


class obj_nexus(obj,obj_defended):
    # change from default iteration over values to keys, blow up
    def __iter__(self): deprecation_error()

    # change from deepcopy to shallow copy, blow up
    def copy(self): deprecation_error()
#end class obj_nexus



class DevBaseNexus(DevBase,obj_defended):

    # change from default iteration over values to keys, blow up
    def __iter__(self): deprecation_error()

    # change from deepcopy to shallow copy, blow up
    def copy(self): deprecation_error()

    # logging - unique to Nexus-style DevBase (future refactor)
    @property
    def _logfile(self):
        return generic_settings.devlog

    def log(self,*a,**kw):
        kw.setdefault('logfile',self._logfile)
        log(*a,**kw)

    def warn(self,msg,indent='    '):
        warn(
            msg,
            indent,
            warn_type = 'class',
            cls       = type(self).__qualname__,
            )

    def error(self,msg,*,header=None,exit=True,trace=-2):
        if header is None:
            header = self.__class__.__name__
        error(msg,header,exit,trace,logfile=self._logfile)

#end class DevBaseNexus



def to_obj(d):
    o = obj()
    for k,v in d.items():
        if hasattr(v,'__dict__'):
            o[k] = to_obj(v)
        else:
            o[k] = v
    return o
#end def to_obj


obj     = obj_nexus
DevBase = DevBaseNexus


# restore original/old obj and DevBase classes
#DevBase = DevBaseDeprecated
#obj     = obj_deprecated
