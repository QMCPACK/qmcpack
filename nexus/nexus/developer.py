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


from .debug import ci, interact

from .generic import NexusError, log, error, warn, message 
from .generic import unavailable, available, Void
from .generic import obj_deprecated, DevBaseDeprecated


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
        if isinstance(v,(obj,DevBase)):
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



class dotdict(dict):
    '''A dictionary with dot-access for keys'''
    def __getattr__(self, item):
        return self[item]
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    def copy(self): return self.__class__(self)
    def __deepcopy__(self, memo):
        result = self.__class__.__new__(self.__class__)
        memo[id(self)] = result
        for key, value in self.items():
            result[copy.deepcopy(key, memo)] = copy.deepcopy(value, memo)
        return result



class obj:
    # dict interface
    @classmethod
    def fromkeys(cls, keys, value=None):
        return cls(dict.fromkeys(keys, value))

    def __init__(self,*args,**kwargs):   self.__dict__.update(dict(*args,**kwargs))
    def items(self):              return self.__dict__.items()
    def clear(self):              return self.__dict__.clear()
    
    # change from deep copy to shallow is pernicious
    #   blow up until purged from all code
    #def copy(self):               return self.__class__(self.__dict__)
    def copy(self):
        raise RuntimeError('shallow copy called by obj!!!')
        return self.__class__(self.__dict__)

    def get(self,*a,**kw):        return self.__dict__.get(*a,**kw)
    def keys(self):               return self.__dict__.keys()
    def pop(self,*a,**kw):        return self.__dict__.pop(*a,**kw)
    def values(self):             return self.__dict__.values()
    def popitem(self,*a,**kw):    return self.__dict__.popitem(*a,**kw)
    def setdefault(self,*a,**kw): return self.__dict__.setdefault(*a,**kw)
    def update(self,*a,**kw):     return self.__dict__.update(*a,**kw)

    # basic functions, including dot access
    def __len__(self):               return len(self.__dict__)
    def __contains__(self,key):      return key in self.__dict__
    def __getitem__(self,key):       return self.__dict__[key]
    def __setitem__(self,key,value): self.__dict__[key]=value
    def __delitem__(self,key):       del self.__dict__[key]
    def __eq__(self,other):          return self.__dict__==other

    # change from default iteration over values to keys, blow up
    #def __iter__(self):
    #    return iter(self.__dict__)
    def __iter__(self):
        raise RuntimeError('obj iteration called!!!')
        for item in self.__dict__: 
            yield item

    # pretty print
    __repr__ = _pp_repr
    __str__  = _pp_str
#end class obj



class DevBase:
    # similar to/same as dict
    def __len__(self):               return len(self.__dict__)
    def __contains__(self,key):      return key in self.__dict__
    def __getitem__(self,key):       return self.__dict__[key]
    def __setitem__(self,key,value): self.__dict__[key]=value
    def __delitem__(self,key):       del self.__dict__[key]
    def keys(self):                  return self.__dict__.keys()
    def values(self):                return self.__dict__.values()
    def items(self):                 return self.__dict__.items()
    def update(self,*a,**kw):        return self.__dict__.update(*a,**kw)
    def clear(self):                 return self.__dict__.clear()

    # correctly iterate over values, not keys
    def __iter__(self):
        raise RuntimeError('DevBase bare iteration!')
        for item in self.__dict__.values():
            yield item

    # (deep) copy
    def copy(self):
        raise RuntimeError('copy called by DevBase!!!')
        return copy.deepcopy(self)

    # pretty print
    __repr__ = _pp_repr
    __str__  = _pp_str


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

#end class DevBase



def to_obj(d):
    o = obj()
    for k,v in d.items():
        if hasattr(v,'__dict__'):
            o[k] = to_obj(v)
        else:
            o[k] = v
    return o
#end def to_obj



# restore original/old obj and DevBase classes
#DevBase = DevBaseDeprecated
#obj     = obj_deprecated
