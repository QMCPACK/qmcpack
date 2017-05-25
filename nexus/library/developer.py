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
#    Void                                                            #
#      Class instances used to represent missing elements.           #
#      Execution stops when any action is performed on a Void object.#
#                                                                    #
#    unavailable                                                     #
#      Function to create named void objects.                        #
#      Used when imported entities do not exist on the local machine.#
#      Allows execution tp proceed normally so long as none of these #
#        non-existent entities are used during runtime execution.    #
#      This enables the maximum amount of Nexus functionality to be  #
#        accessed given the available modules.                       #
#                                                                    #
#====================================================================#


from generic import obj,object_interface,log,error,warn,devlog
from debug import ci,interact


class DevBase(obj):
    def not_implemented(self):
        self.error('a base class function has not been implemented',trace=True)
    #end def not_implemented
#end class DevBase


class enum(object_interface):
    def __init__(self,*keys):
        if len(keys)==1 and isinstance(keys[0],(list,tuple)):
            keys = keys[0]
        #end if
        n=0
        for key in keys:
            self[key] = n
            n+=1
        #end for
    #end def __init__
        
    # override some object interface methods
    # specifically forbid modification
    def __setitem__(self,name,value):
        self._error('attempted to modify immutable enum object')
    #end def __setitem__

    def __delitem__(self,name):
        self._error('attempted to modify immutable enum object')
    #end def __delitem__

    def clear(self):
        self._error('attempted to modify immutable enum object')
    #end def clear

    def _clear(self):
        enum.clear(self)
    #end def _clear
#end class enum


class Void:
    void_items = dict()

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
        DevBase.class_error(msg,'Void')
    #end def _unavailable

        
    @classmethod
    def _class_unavailable(cls):
        msg = 'encountered a void item from an unavailable module'
        DevBase.class_error(msg,'Void')
    #end def _class_unavailable
        


    def __init__(self,module=None,item=None):
        Void.void_items[id(self)] = module,item
    #end def __init__


    #list of magic functions taken from the following sources
    #  http://web.archive.org/web/20110131211638/http://diveintopython3.org/special-method-names.html
    #  http://www.ironpythoninaction.com/magic-methods.html


    #class methods
    @classmethod
    def __instancecheck__(cls,*args,**kwargs):
        Void._class_unavailable()
    @classmethod
    def __subclasscheck__(cls,*args,**kwargs):
        Void._class_unavailable()
    @classmethod
    def __subclasshook__(cls,*args,**kwargs):
        Void._class_unavailable()
    

    #member methods
    def __new__(self,*args,**kwargs):
        Void._unavailable(self)
    def __eq__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ne__(self,*args,**kwargs):
        Void._unavailable(self)
    def __lt__(self,*args,**kwargs):
        Void._unavailable(self)
    def __le__(self,*args,**kwargs):
        Void._unavailable(self)
    def __gt__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ge__(self,*args,**kwargs):
        Void._unavailable(self)
    def __nonzero__(self,*args,**kwargs):
        Void._unavailable(self)
    def __subclasses__(self,*args,**kwargs):
        Void._unavailable(self)
    def __call__(self,*args,**kwargs):
        Void._unavailable(self)
    def __hash__(self,*args,**kwargs):
        Void._unavailable(self)
    #def __del__(self,*args,**kwargs):
    #    Void._unavailable(self)
    def __dir__(self,*args,**kwargs):
        Void._unavailable(self)
    def __getitem__(self,*args,**kwargs):
        Void._unavailable(self)
    def __setitem__(self,*args,**kwargs):
        Void._unavailable(self)
    def __delitem__(self,*args,**kwargs):
        Void._unavailable(self)
    def __len__(self,*args,**kwargs):
        Void._unavailable(self)
    def __contains__(self,*args,**kwargs):
        Void._unavailable(self)
    def __iter__(self,*args,**kwargs):
        Void._unavailable(self)
    def __reversed__(self,*args,**kwargs):
        Void._unavailable(self)
    def __missing__(self,*args,**kwargs):
        Void._unavailable(self)
    def __length_hint__(self,*args,**kwargs):
        Void._unavailable(self)
    def __repr__(self,*args,**kwargs):
        Void._unavailable(self)
    def __str__(self,*args,**kwargs):
        Void._unavailable(self)
    def __unicode__(self,*args,**kwargs):
        Void._unavailable(self)
    def __getattr__(self,*args,**kwargs):
        Void._unavailable(self)
    def __setattr__(self,*args,**kwargs):
        Void._unavailable(self)
    def __delattr__(self,*args,**kwargs):
        Void._unavailable(self)
    def __getattribute__(self,*args,**kwargs):
        Void._unavailable(self)
    def __add__(self,*args,**kwargs):
        Void._unavailable(self)
    def __sub__(self,*args,**kwargs):
        Void._unavailable(self)
    def __mul__(self,*args,**kwargs):
        Void._unavailable(self)
    def __floordiv__(self,*args,**kwargs):
        Void._unavailable(self)
    def __div__(self,*args,**kwargs):
        Void._unavailable(self)
    def __truediv__(self,*args,**kwargs):
        Void._unavailable(self)
    def __mod__(self,*args,**kwargs):
        Void._unavailable(self)
    def __divmod__(self,*args,**kwargs):
        Void._unavailable(self)
    def __pow__(self,*args,**kwargs):
        Void._unavailable(self)
    def __lshift__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rshift__(self,*args,**kwargs):
        Void._unavailable(self)
    def __and__(self,*args,**kwargs):
        Void._unavailable(self)
    def __xor__(self,*args,**kwargs):
        Void._unavailable(self)
    def __or__(self,*args,**kwargs):
        Void._unavailable(self)
    def __neg__(self,*args,**kwargs):
        Void._unavailable(self)
    def __pos__(self,*args,**kwargs):
        Void._unavailable(self)
    def __abs__(self,*args,**kwargs):
        Void._unavailable(self)
    def __invert__(self,*args,**kwargs):
        Void._unavailable(self)
    def __complex__(self,*args,**kwargs):
        Void._unavailable(self)
    def __int__(self,*args,**kwargs):
        Void._unavailable(self)
    def __long__(self,*args,**kwargs):
        Void._unavailable(self)
    def __float__(self,*args,**kwargs):
        Void._unavailable(self)
    def __oct__(self,*args,**kwargs):
        Void._unavailable(self)
    def __hex__(self,*args,**kwargs):
        Void._unavailable(self)
    def __index__(self,*args,**kwargs):
        Void._unavailable(self)
    def __enter__(self,*args,**kwargs):
        Void._unavailable(self)
    def __exit__(self,*args,**kwargs):
        Void._unavailable(self)
    def __get__(self,*args,**kwargs):
        Void._unavailable(self)
    def __set__(self,*args,**kwargs):
        Void._unavailable(self)
    def __delete__(self,*args,**kwargs):
        Void._unavailable(self)
    def __doc__(self,*args,**kwargs):
        Void._unavailable(self)
    def __dict__(self,*args,**kwargs):
        Void._unavailable(self)
    def __slots__(self,*args,**kwargs):
        Void._unavailable(self)
    def __class__(self,*args,**kwargs):
        Void._unavailable(self)
    def __bases__(self,*args,**kwargs):
        Void._unavailable(self)
    def __name__(self,*args,**kwargs):
        Void._unavailable(self)
    def __all__(self,*args,**kwargs):
        Void._unavailable(self)
    def __file__(self,*args,**kwargs):
        Void._unavailable(self)
    def __module__(self,*args,**kwargs):
        Void._unavailable(self)
    #def __metaclass__(self,*args,**kwargs):
    #    Void._unavailable(self)
    def __import__(self,*args,**kwargs):
        Void._unavailable(self)
    def __radd__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rsub__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rmul__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rtruediv__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rfloordiv__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rmod__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rdivmod__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rpow__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rlshift__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rrshift__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rand__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rxor__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ror__(self,*args,**kwargs):
        Void._unavailable(self)
    def __iadd__(self,*args,**kwargs):
        Void._unavailable(self)
    def __isub__(self,*args,**kwargs):
        Void._unavailable(self)
    def __imul__(self,*args,**kwargs):
        Void._unavailable(self)
    def __itruediv__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ifloordiv__(self,*args,**kwargs):
        Void._unavailable(self)
    def __imod__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ipow__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ilshift__(self,*args,**kwargs):
        Void._unavailable(self)
    def __irshift__(self,*args,**kwargs):
        Void._unavailable(self)
    def __iand__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ixor__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ior__(self,*args,**kwargs):
        Void._unavailable(self)
    def __round__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ceil__(self,*args,**kwargs):
        Void._unavailable(self)
    def __floor__(self,*args,**kwargs):
        Void._unavailable(self)
    def __trunc__(self,*args,**kwargs):
        Void._unavailable(self)
    def __bool__(self,*args,**kwargs):
        Void._unavailable(self)
    def __copy__(self,*args,**kwargs):
        Void._unavailable(self)
    def __deepcopy__(self,*args,**kwargs):
        Void._unavailable(self)
    def __getstate__(self,*args,**kwargs):
        Void._unavailable(self)
    def __reduce__(self,*args,**kwargs):
        Void._unavailable(self)
    def __reduce_ex__(self,*args,**kwargs):
        Void._unavailable(self)
    def __getnewargs__(self,*args,**kwargs):
        Void._unavailable(self)
    def __setstate__(self,*args,**kwargs):
        Void._unavailable(self)
    def __bytes__(self,*args,**kwargs):
        Void._unavailable(self)
    def __format__(self,*args,**kwargs):
        Void._unavailable(self)
    def __next__(self,*args,**kwargs):
        Void._unavailable(self)
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
