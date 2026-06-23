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


from .utilities import to_str
from .generic import obj, object_interface, hidden, NexusError, log, error, warn, message
from .debug import ci, interact


class DevBase(obj):
    def not_implemented(self,name=None):
        if name is None:
            msg = 'a member function has not been implemented for class "{}"'.format(self.__class__.__name__)
        else:
            msg = 'member function "{}" has not been implemented for class "{}"'.format(name,self.__class__.__name__)
        self.error(msg,trace=True)
    #end def not_implemented
#end class DevBase



class Void(object):
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
    #def __new__(self,*args,**kwargs):
    #    Void._unavailable(self)
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
    #def __slots__(self,*args,**kwargs):
    #    Void._unavailable(self)
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


def valid_variable_name(s):
    """Check if a variable name contains invalid characters"""
    if not any([i in ('!"#$%&\'()*+,-./:;<=>?@[\\]^`{|}-\n\t ') for i in s]):
        return True
    else:
        return False
    #end if
#end def valid_variable_name









import sys
import copy
import pickle
from .generic import sorted_generic, generic_settings



def save(o,filepath):
    with open(filepath,'wb') as f:
        binary = pickle.HIGHEST_PROTOCOL
        pickle.dump(o,f,binary)

def load(filepath):
    with open(filepath,'rb') as f:
        dl = pickle.load(f)
    return dl



def _pp_repr(self):
    s=''
    for k in sorted_generic(self.keys()):
        v = self.__dict__[k]
        if hasattr(v,'__class__'):
            s+='  {0:<20}  {1:<20}\n'.format(str(k),v.__class__.__name__)
        else:
            s+='  {0:<20}  {1:<20}\n'.format(str(k),type(v))
    return s


def _pp_str(self,nindent=1):
    pad = '  '
    npad = nindent*pad
    s=''
    normal = []
    qable  = []
    for k,v in self.items():
        if isinstance(v,(obj2,DevBase2)):
            qable.append(k)
        else:
            normal.append(k)
    normal = sorted_generic(normal)
    qable  = sorted_generic(qable)
    indent = npad+18*' '
    for k in normal:
        v = self[k]
        vstr = str(v).replace('\n','\n'+indent)
        s+=npad+'{0:<15} = '.format(str(k))+vstr+'\n'
    for k in qable:
        v = self[k]
        s+=npad+str(k)+'\n'
        s+=v.__str__(nindent+1)
        if isinstance(k,str):
            s+=npad+'end '+k+'\n'
    return s



class obj2:
    # dict interface
    def __init__(self,*args,**kwargs):   self.__dict__.update(dict(*args,**kwargs))
    def items(self):              return self.__dict__.items()
    def clear(self):              return self.__dict__.clear()
    
    # change from deep copy to shallow is pernicious
    #   nuke it until purged from all code
    #   probably should fully remove DevBase2.copy in a second pass
    #def copy(self):               return self.__class__(self.__dict__)
    def copy(self):
        raise RuntimeError('shallow copy called by obj!!!')
        return self.__class__(self.__dict__)

    def fromkeys(self,*a,**kw):   return self.__class__(self.__dict__.fromkeys(*a,**kw))
    def get(self,*a,**kw):        return self.__dict__.get(*a,**kw)
    def keys(self):               return self.__dict__.keys()
    def pop(self,*a,**kw):        return self.__dict__.pop(*a,**kw)
    def values(self):             return self.__dict__.values()
    def popitem(self,*a,**kw):    return self.__dict__.popitem(*a,**kw)
    def setdefault(self,*a,**kw): return self.__dict__.setdefault(*a,**kw)
    def update(self,*a,**kw):     return self.__dict__.update(*a,**kw)

    # protected dict interface
    def _items(self):              return self.__dict__.items()
    def _clear(self):              return self.__dict__.clear()
    def _copy(self):               return self.__class__(self.__dict__)
    def _fromkeys(self,*a,**kw):   return self.__class__(self.__dict__.fromkeys(*a,**kw))
    def _get(self,*a,**kw):        return self.__dict__.get(*a,**kw)
    def _keys(self):               return self.__dict__.keys()
    def _pop(self,*a,**kw):        return self.__dict__.pop(*a,**kw)
    def _values(self):             return self.__dict__.values()
    def _popitem(self,*a,**kw):    return self.__dict__.popitem(*a,**kw)
    def _setdefault(self,*a,**kw): return self.__dict__.setdefault(*a,**kw)
    def _update(self,*a,**kw):     return self.__dict__.update(*a,**kw)

    # basic functions, including dot access
    def __len__(self):               return len(self.__dict__)
    def __contains__(self,key):      return key in self.__dict__
    def __getitem__(self,key):       return self.__dict__[key]
    def __setitem__(self,key,value): self.__dict__[key]=value
    def __delitem__(self,key):       del self.__dict__[key]

    # iter also diverges, blow up
    def __iter__(self):
        raise RuntimeError('obj iteration called!!!')
        for item in self.__dict__: 
            yield item

    # pretty print
    __repr__ = _pp_repr
    __str__  = _pp_str



class DevBase2:
    # similar to/same as dict
    def __len__(self):          return len(self.__dict__)
    def __contains__(self,key): return key in self.__dict__
    def __getitem__(self,key):  return self.__dict__[key]
    def __setitem__(self,key,value): self.__dict__[key]=value
    def __delitem__(self,key):       del self.__dict__[key]
    def keys(self):             return self.__dict__.keys()
    def values(self):           return self.__dict__.values()
    def items(self):            return self.__dict__.items()
    def update(self,*a,**kw):   return self.__dict__.update(*a,**kw)
    def clear(self):            return self.__dict__.clear()

    # correctly iterate over values, not keys
    def __iter__(self):
        for item in self.__dict__.values():
            yield item

    # (deep) copy
    def copy(self):
        return copy.deepcopy(self)

    # pretty print
    __repr__ = _pp_repr
    __str__  = _pp_str

    # save and load
    def save(self,filepath=None):
        if filepath is None:
            filepath='./'+self.__class__.__name__+'.p'
        save(self,filepath)

    def load(self,filepath=None):
        if filepath is None:
            filepath='./'+self.__class__.__name__+'.p'
        tmp = load(filepath)
        d = self.__dict__
        d.clear()
        for k,v in tmp.__dict__.items():
            d[k] = v

    # logging
    @property
    def _logfile(self):
        return generic_settings.devlog

    def log(self,*a,**kw):
        kw.setdefault('logfile',self._logfile)
        log(*a,**kw)

    def warn(self,message,header=None):
        if header is None:
            header=self.__class__.__name__
        warn(message,header,logfile=self._logfile)

    def error(self,message,header=None,exit=True,trace=-2):
        if header is None:
            header = self.__class__.__name__
        error(message,header,exit,trace,self._logfile)

    # general dev, ditch?
    def not_implemented(self,name=None):
        if name is None:
            msg = 'a member function has not been implemented for class "{}"'.format(self.__class__.__name__)
        else:
            msg = 'member function "{}" has not been implemented for class "{}"'.format(name,self.__class__.__name__)
        self.error(msg,trace=True)




def to_obj(d):
    o = obj()
    for k,v in d.items():
        if hasattr(v,'__dict__'):
            o[k] = to_obj(v)
        else:
            o[k] = v
    return o
