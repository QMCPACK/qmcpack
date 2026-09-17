##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  developer.py                                                      #
#    Defines developer environment.  Supplies base class for generic #
#    development (Nexus or beyond)                                   #
#                                                                    #
#  Content summary:                                                  #
#    log, warn,                                                      #
#      Function interface to logging and warning handling.           #
#                                                                    #
#    DevBase                                                         #
#      Base class inheriting generic abilities for obj, etc.         #
#      Allows for unimplemented functions.                           #
#                                                                    #
#====================================================================#


from typing import NoReturn

from .developer_tools import save,load,_pp_repr,_pp_str,dotdict,obj,DevBase  # noqa: F401
from .developer_tools import sorted_py2 as sorted_generic
from .debug import ci, interact  # noqa: F401


from .generic import NexusError, FileFormatError, NotAnElementError  # noqa: F401
from .generic import error, nxs_print, warn, message  # noqa: F401


import traceback


def deprecation_error():
    msg = (
        'A now-deprecated member function of obj has been called.\n'
        'Please report this issue to the Nexus developers immediately.\n'
        )
    highlight = '='*79
    stack = ''.join(traceback.format_stack()[:-1])
    report = (
        f'\n{highlight}\n{msg}\n{highlight}\n'
        f'{highlight}\nTraceback (most recent call last):\n{stack}{highlight}\n'
        f'{highlight}\n{msg}\n{highlight}'
        )
    raise RuntimeError(report)
#end def deprecation_error


class DevBaseNexus(DevBase):

    # change from default iteration over values to keys, blow up
    def __iter__(self): deprecation_error()

    # change from deepcopy to shallow copy, blow up
    def copy(self): deprecation_error()


    def nxs_print(self,*a,**kw):
        nxs_print(*a,**kw)

    def warn(self,msg,indent='    '):
        warn(
            msg,
            indent,
            warn_type = 'class',
            cls       = type(self).__qualname__,
            )

    def error(self, msg: str, *, header: str | None = None) -> NoReturn:
        if header is None:
            header = type(self).__name__
        error(msg, header)

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

DevBase = DevBaseNexus
