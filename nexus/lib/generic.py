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
#    generic                                                         #
#      More efficient implementation of AllAbilities+obj interface.  #
#      Intended to allow for method namespace infringement without   #
#        loss of access to functionality.                            #
#      Limited use so far.                                           #
#                                                                    #
#    hidden                                                          #
#      Like generic, but allows for hidden storage.                  #
#      Can be used, e.g., to make an ordered object class.           #
#      See use in qmcpack_input.py.                                  #
#                                                                    #
#====================================================================#


import sys
import traceback
from copy import deepcopy
import pickle
from random import randint

from utilities import sorted_py2


class generic_settings:
    devlog      = sys.stdout
    raise_error = False
#end class generic_settings


class NexusError(Exception):
    None
#end class NexusError


exit_call = sys.exit


def nocopy(value):
    return value
#end def nocopy



sorted_generic = sorted_py2



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


def warn(msg,header=None,indent='    ',logfile=None):
    if logfile is None:
        logfile = generic_settings.devlog
    #end if
    post_header=' warning:'
    message(msg,header,post_header,indent,logfile)
#end def warn


def error(msg,header=None,exit=True,trace=True,indent='    ',logfile=None):
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
        if trace:
            traceback.print_stack()
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
            eqval = svar==ovar
            if isinstance(eqval,bool):
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
        except:
            try:
                tmp = pickle.load(fobj,encoding='latin1')
            except:
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

    def warn(self,message,header=None):
        if header is None:
            header=self.__class__.__name__
        #end if
        warn(message,header,logfile=self._logfile)
    #end def warn

    def error(self,message,header=None,exit=True,trace=True):
        if header==None:
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
        if header==None:
            header=cls.__name__
        #end if
        warn(message,header,logfile=cls._logfile)
    #end def class_warn

    @classmethod
    def class_error(cls,message,header=None,exit=True,trace=True,post_header=' Error:'):
        if header==None:
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



class obj(object_interface):

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
        return tuple(obj.list(self,*keys))
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
            if isinstance(v,obj):
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
        o = obj()
        for k in keys:
            o[k] = self[k]
        #end for
        return o
    #end def obj

    def to_obj(self):
        o = obj()
        for k,v in self._items():
            if isinstance(v,obj):
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
            obj.error(self,'a required key is not present\nkey required: {0}\nkeys present: {1}'.format(key,self._sorted_keys()))
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
            obj.error(self,'a required key is not present\nkey required: {0}\nkeys present: {1}'.format(key,self._sorted_keys()))
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
        obj.transfer_from(self,other,keys,copy=deep)
    #end def copy_from

    def copy_to(self,other,keys=None,deep=True):
        obj.transfer_to(self,other,keys,copy=deep)
    #end def copy_to

    def extract(self,keys=None,optional=False):
        ext = obj()
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
        if isinstance(path,str):
            path = path.split('/')
        #end if
        for p in path:
            if not p in o:
                return False
            #end if
            o = o[p]
        #end for
        return True
    #end def path_exists

    def set_path(self,path,value=None):
        o = self
        cls = self.__class__
        if isinstance(path,str):
            path = path.split('/')
        #end if
        for p in path[0:-1]:
            if not p in o:
                o[p] = cls()
            #end if
            o = o[p]
        #end for
        o[path[-1]] = value
    #end def set_path

    def get_path(self,path,value=None):
        o = self
        if isinstance(path,str):
            path = path.split('/')
        #end if
        for p in path[0:-1]:
            if not p in o:
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
            s = obj()
            path = ''
        #end if
        for k,v in self._items():
            p = path+str(k)
            if isinstance(v,obj):
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
        obj.append(self,*args,**kwargs)
    #  return representations
    def _list(self,*args,**kwargs):
        return obj.list(self,*args,**kwargs)
    def _list_optional(self,*args,**kwargs):
        return obj.list_optional(self,*args,**kwargs)
    def _tuple(self,*args,**kwargs):
        return obj.tuple(self,*args,**kwargs)
    def _dict(self,*args,**kwargs):
        return obj.dict(self,*args,**kwargs)
    def _to_dict(self,*args,**kwargs):
        return obj.to_dict(self,*args,**kwargs)
    def _obj(self,*args,**kwargs):
        return obj.obj(self,*args,**kwargs)
    def _to_obj(self,*args,**kwargs):
        return obj.to_obj(self,*args,**kwargs)
    #  list extensions
    def _first(self,*args,**kwargs):
        return obj.first(self,*args,**kwargs)
    def _last(self,*args,**kwargs):
        return obj.last(self,*args,**kwargs)
    def _select_random(self,*args,**kwargs):
        return obj.select_random(self,*args,**kwargs)
    #  dict extensions
    def _sorted_keys(self,*args,**kwargs):
        return obj.sorted_keys(self,*args,**kwargs)
    def _random_key(self,*args,**kwargs):
        obj.random_key(self,*args,**kwargs)
    def _set(self,*args,**kwargs):
        obj.set(self,*args,**kwargs)
    def _set_optional(self,*args,**kwargs):
        obj.set_optional(self,*args,**kwargs)
    def _get(self,*args,**kwargs):
        obj.get(self,*args,**kwargs)
    def _get_optional(self,*args,**kwargs):
        obj.get_optional(self,*args,**kwargs)
    def _get_required(self,*args,**kwargs):
        obj.get_required(self,*args,**kwargs)
    def _delete(self,*args,**kwargs):
        obj.delete(self,*args,**kwargs)
    def _delete_optional(self,*args,**kwargs):
        obj.delete_optional(self,*args,**kwargs)
    def _delete_required(self,*args,**kwargs):
        obj.delete_required(self,*args,**kwargs)
    def _add(self,*args,**kwargs):
        obj.add(self,*args,**kwargs)
    def _add_optional(self,*args,**kwargs):
        obj.add_optional(self,*args,**kwargs)
    def _transfer_from(self,*args,**kwargs):
        obj.transfer_from(self,*args,**kwargs)
    def _transfer_to(self,*args,**kwargs):
        obj.transfer_to(self,*args,**kwargs)
    def _move_from(self,*args,**kwargs):
        obj.move_from(self,*args,**kwargs)
    def _move_to(self,*args,**kwargs):
        obj.move_to(self,*args,**kwargs)
    def _move_from_optional(self,*args,**kwargs):
        obj.move_from_optional(self,*args,**kwargs)
    def _move_to_optional(self,*args,**kwargs):
        obj.move_to_optional(self,*args,**kwargs)
    def _copy_from(self,*args,**kwargs):
        obj.copy_from(self,*args,**kwargs)
    def _copy_to(self,*args,**kwargs):
        obj.copy_to(self,*args,**kwargs)
    def _extract(self,*args,**kwargs):
        obj.extract(self,*args,**kwargs)
    def _extract_optional(self,*args,**kwargs):
        obj.extract_optional(self,*args,**kwargs)
    def _check_required(self,*args,**kwargs):
        obj.check_required(self,*args,**kwargs)
    def _check_types(self,*args,**kwargs):
        obj.check_types(self,*args,**kwargs)
    def _check_types_optional(self,*args,**kwargs):
        obj.check_types_optional(self,*args,**kwargs)
    def _shallow_copy(self,*args,**kwargs):
        obj.shallow_copy(self,*args,**kwargs)
    def _inverse(self,*args,**kwargs):
        return obj.inverse(self,*args,**kwargs)
    def _path_exists(self,*args,**kwargs):
        obj.path_exists(self,*args,**kwargs)
    def _set_path(self,*args,**kwargs):
        obj.set_path(self,*args,**kwargs)
    def _get_path(self,*args,**kwargs):
        obj.get_path(self,*args,**kwargs)
    def _serial(self,*args,**kwargs):
        return obj.serial(self,*args,**kwargs)

#end class obj





class hobj(obj):
    def __init__(self,*args,**kwargs):
        obj.__init__(self,*args,**kwargs)
    #end def __init__

    @property
    def _dict(self):
        return self.__dict__
    #end def _dict

    @property
    def _alt(self):
        return self.__dict__
    #end def _alt

    def __len__(self):
        return len(self._dict)
    #end def __len__

    def __contains__(self,name):
        return name in self._dict
    #end def __contains__

    def __getitem__(self,name):
        return self._dict[name]
    #end def __getitem__

    def __setitem__(self,name,value):
        self._dict[name] = value
    #end def __setitem__

    def __delitem__(self,name):
        del self._dict[name]
    #end def __delitem__

    def __iter__(self):
        d = self._dict
        for item in d.__dict__:
            yield d[item]
        #end for
    #end def __iter__

    def keys(self):
        return self._dict.keys()
    #end def keys

    def values(self):
        return self._dict.values()
    #end def keys

    def items(self):
        return self._dict.items()
    #end def items

    def clear(self):
        self._dict.clear()
    #end def clear

    # access preserving functions
    #  dict interface
    def _keys(self,*args,**kwargs):
        return hobj.keys(self,*args,**kwargs)
    def _values(self,*args,**kwargs):
        hobj.values(self,*args,**kwargs)
    def _items(self,*args,**kwargs):         
        return hobj.items(self,*args,**kwargs)         
    def _clear(self,*args,**kwargs):
        hobj.clear(self,*args,**kwargs)
#end class hobj



class hidden(hobj):
    def __init__(self,*vals,**kwargs):
        d = object.__getattribute__(self,'__dict__')
        d['_hidden_'] = hobj()
        d['_public_'] = hobj()
        hobj.__init__(self,*vals,**kwargs)
    #end def __init__

    @property
    def _dict(self):
        return self.__dict__['_public_']
    #end def __get_dict

    @property
    def _alt(self):
        return self.__dict__['_hidden_']
    #end def __alt

    def __getattribute__(self,name):
        d = object.__getattribute__(self,'__dict__')
        if '_public_' in d:
            p = d['_public_']
            if name in p:
                return p[name]
            else:
                return object.__getattribute__(self,name)
            #end if
        else:
            return object.__getattribute__(self,name)
        #end if
    #end def __getattribute__

    def __setattr__(self,name,value):
        self._dict[name] = value
    #end def __setattr__

    def __delattr__(self,name):
        del self._dict[name]
    #end def __delattr__

    def hidden(self):
        return self.__dict__['_hidden_']
    #end def hidden

    def public(self):
        return self.__dict__['_public_']
    #end def public

    def _hidden(self):
        return hidden.hidden(self)
    #end def _hidden

    def _public(self):
        return hidden.public(self)
    #end def _public

    def open_log(self,filepath):
        self._alt._open_log(filepath)
    #end def open_log

    def close_log(self):
        self._alt._close_log()
    #end def close_log

    def write(self,s):
        self._alt._write(s)
    #end def write

    def log(self,*items,**kwargs):
        self._alt._log(*items,**kwargs)
    #end def log

    def __repr__(self):
        s=''
        for k in self._sorted_keys():
            if not isinstance(k,str) or k[0]!='_':
                v=self._dict[k]
                if hasattr(v,'__class__'):
                    s+='  {0:<20}  {1:<20}\n'.format(k,v.__class__.__name__)
                else:
                    s+='  {0:<20}  {1:<20}\n'.format(k,type(v))
                #end if
            #end if
        #end for
        return s
    #end def __repr__

    #  log, warning, and error messages
    def _open_log(self,*args,**kwargs):
        hidden.open_log(self,*args,**kwargs)
    def _close_log(self,*args,**kwargs):
        hidden.close_log(self,*args,**kwargs)
    def _write(self,*args,**kwargs):
        hidden.write(self,*args,**kwargs)
    def _log(self,*args,**kwargs):
        hidden.log(self,*args,**kwargs)

#end class hidden
