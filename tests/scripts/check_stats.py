#! /usr/bin/env python3

from __future__ import print_function, division

# Statical error checking code for use by testing framework (stat.h5 files)
# Jaron Krogel/ORNL


# check_stats.py packages obj and HDFreader classes from Nexus.
#   Note that h5py is required (which depends on numpy).



######################################################################
# from generic.py
######################################################################

import sys
import traceback
from copy import deepcopy
import pickle
from random import randint


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



def log(*items,**kwargs):
    indent=None
    logfile=generic_settings.devlog
    if len(kwargs)>0:
        n=0
        if 'indent' in kwargs:
            indent = kwargs['indent']
            n+=1
        #end if
        if 'logfile' in kwargs:
            logfile = kwargs['logfile']
            n+=1
        #end if
        if n!=len(kwargs):
            valid = 'indent logfile'.split()
            invalid = set(kwargs.keys())-set(valid)
            error('invalid keyword arguments provided\ninvalid keywords: {0}\nvalid options are: {1}'.format(sorted(invalid),valid))
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
        for k in sorted(self._keys()):
            if not isinstance(k,str) or k[0]!='_':
                v=self.__dict__[k]
                if hasattr(v,'__class__'):
                    s+='  {0:<20}  {1:<20}\n'.format(k,v.__class__.__name__)
                else:
                    s+='  {0:<20}  {1:<20}\n'.format(k,type(v))
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
        normal.sort()
        qable.sort()
        indent = npad+18*' '
        for k in normal:
            v = self[k]
            vstr = str(v).replace('\n','\n'+indent)
            s+=npad+'{0:<15} = '.format(k)+vstr+'\n'
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
        if fpath==None:
            fpath='./'+self.__class__.__name__+'.p'
        #end if
        fobj = open(fpath,'w')
        binary = pickle.HIGHEST_PROTOCOL
        pickle.dump(self,fobj,binary)
        fobj.close()
        del fobj
        del binary
        return
    #end def save

    def load(self,fpath=None):
        if fpath==None:
            fpath='./'+self.__class__.__name__+'.p'
        #end if
        fobj = open(fpath,'r')
        tmp = pickle.load(fobj)
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
            keys = sorted(self._keys())
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
            keys = sorted(self._keys())
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
            keys = sorted(self._keys())
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
            keys = sorted(self._keys())
        elif nkeys==1 and isinstance(keys[0],(list,tuple)):
            keys = keys[0]
        #end if
        o = obj()
        for k in keys:
            o[k] = self[k]
        #end for
        return o
    #end def obj


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
    def random_key(self):
        key = None
        nkeys = len(self)
        if nkeys>0:
            key = self._keys()[randint(0,nkeys-1)]
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
            obj.error(self,'a required key is not present\nkey required: {0}\nkeys present: {1}'.format(key,sorted(self._keys())))
        #end if
        return value
    #end def get_required

    def delete(self,*keys):
        nkeys = len(keys)
        single = False
        if nkeys==0:
            keys = sorted(self._keys())
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
            obj.error(self,'a required key is not present\nkey required: {0}\nkeys present: {1}'.format(key,sorted(self._keys())))
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
        if keys==None:
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
        if keys==None:
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

    def move_from(self,other,keys=None):
        if keys==None:
            if isinstance(other,object_interface):
                keys = other._keys()
            else:
                keys = other.keys()
            #end if
        #end if
        for k in keys:
            self[k]=other[k]
            del other[k]
        #end for
    #end def move_from

    def move_to(self,other,keys=None):
        if keys==None:
            keys = self._keys()
        #end if
        for k in keys:
            other[k]=self[k]
            del self[k]
        #end for
    #end def move_to

    def copy_from(self,other,keys=None,deep=True):
        obj.transfer_from(self,other,keys,copy=deep)
    #end def copy_from

    def copy_to(self,other,keys=None,deep=True):
        obj.transfer_to(self,other,keys,copy=deep)
    #end def copy_to

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
    #  list extensions
    def _first(self,*args,**kwargs):
        return obj.first(self,*args,**kwargs)
    def _last(self,*args,**kwargs):
        return obj.last(self,*args,**kwargs)
    def _select_random(self,*args,**kwargs):
        return obj.select_random(self,*args,**kwargs)
    #  dict extensions
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
    def _copy_from(self,*args,**kwargs):
        obj.copy_from(self,*args,**kwargs)
    def _copy_to(self,*args,**kwargs):
        obj.copy_to(self,*args,**kwargs)
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

######################################################################
# end from generic.py
######################################################################


######################################################################
# from superstring.py
######################################################################

import string

def contains_any(str, set):
    for c in set:
        if c in str: return 1;
    return 0;
#end def contains_any

invalid_variable_name_chars=set('!"#$%&\'()*+,-./:;<=>?@[\\]^`{|}-\n\t ')
def valid_variable_name(s):
    return not contains_any(s,invalid_variable_name_chars)
#end def valid_variable_name

######################################################################
# end from superstring.py
######################################################################


######################################################################
# from debug.py
######################################################################

import code
import inspect

def ci(locs=None,globs=None):
    if locs is None or globs is None:
        cur_frame = inspect.currentframe()
        caller_frame = cur_frame.f_back
        locs  = caller_frame.f_locals
        globs = caller_frame.f_globals
    #end if
    code.interact(local=dict(globs,**locs))
#end def ci

ls = locals
gs = globals
interact = ci

######################################################################
# end from debug.py
######################################################################


######################################################################
# from developer.py
######################################################################

class DevBase(obj):
    def not_implemented(self):
        self.error('a base class function has not been implemented',trace=True)
    #end def not_implemented
#end class DevBase

######################################################################
# end from developer.py
######################################################################


######################################################################
# from hdfreader.py
######################################################################
from numpy import array,ndarray,minimum,abs,ix_,resize
import sys
import keyword
from inspect import getmembers
import h5py


class HDFglobals(DevBase):
    view = False
#end class HDFglobals


class HDFgroup(DevBase):
    def _escape_name(self,name):
        if name in self._escape_names:
            name=name+'_'
        #end if
        return name
    #end def escape_name

    def _set_parent(self,parent):
        self._parent=parent
        return
    #end def set_parent

    def _add_dataset(self,name,dataset):
        self._datasets[name]=dataset
        return 
    #end def add_dataset

    def _add_group(self,name,group):
        group._name=name
        self._groups[name]=group
        return 
    #end def add_group

    def _contains_group(self,name):
        return name in self._groups.keys()
    #end def _contains_group

    def _contains_dataset(self,name):
        return name in self._datasets.keys()
    #end def _contains_dataset

    def _to_string(self):
        s=''
        if len(self._datasets)>0:
            s+='  datasets:\n'
            for k,v in self._datasets.items():
                s+= '    '+k+'\n'
            #end for
        #end if
        if len(self._groups)>0:
            s+= '  groups:\n'
            for k,v in self._groups.items():
                s+= '    '+k+'\n'
            #end for
        #end if
        return s
    #end def list

#    def __str__(self):
#        return self._to_string()
#    #end def __str__
#
#    def __repr__(self):
#        return self._to_string()
#    #end def __repr__

    def __init__(self):
        self._name=''
        self._parent=None
        self._groups={};
        self._datasets={};
        self._group_counts={}

        self._escape_names=None
        self._escape_names=set(dict(getmembers(self)).keys()) | set(keyword.kwlist)
        return
    #end def __init__


    def _remove_hidden(self,deep=True):
        if '_parent' in self:
            del self._parent
        #end if
        if deep:
            for name,value in self.items():
                if isinstance(value,HDFgroup):
                    value._remove_hidden()
                #end if
            #end for
        #end if
        for name in list(self.keys()):
            if name[0]=='_':
                del self[name]
            #end if
        #end for
    #end def _remove_hidden


    # read in all data views (h5py datasets) into arrays
    #   useful for converting a single group read in view form to full arrays
    def read_arrays(self):
        self._remove_hidden()
        for k,v in self.items():
            if isinstance(v,HDFgroup):
                v.read_arrays()
            else:
                self[k] = array(v)
            #end if
        #end for
    #end def read_arrays


    def get_keys(self):
        if '_groups' in self:
            keys = list(self._groups.keys())
        else:
            keys = list(self.keys())
        #end if
        return keys
    #end def get_keys
#end class HDFgroup




class HDFreader(DevBase):
    datasets = set(["<class 'h5py.highlevel.Dataset'>","<class 'h5py._hl.dataset.Dataset'>"])
    groups   = set(["<class 'h5py.highlevel.Group'>","<class 'h5py._hl.group.Group'>"])
    
    def __init__(self,fpath,verbose=False,view=False):
        
        HDFglobals.view = view

        if verbose:
            print('  Initializing HDFreader')

        self.fpath=fpath
        if verbose:
            print('    loading h5 file')

        try:
            self.hdf = h5py.File(fpath,'r')
        except IOError:
            self._success = False
            self.hdf = obj(obj=obj())
        else:
            self._success = True
        #end if

        if verbose:
            print('    converting h5 file to dynamic object')
        #convert the hdf 'dict' into a dynamic object
        self.nlevels=1
        self.ilevel=0
        #  Set the current hdf group
        self.obj = HDFgroup()
        self.cur=[self.obj]
        self.hcur=[self.hdf]

        if self._success:
            cur   = self.cur[self.ilevel]
            hcur  = self.hcur[self.ilevel]
            for kr,v in hcur.items():
                k=cur._escape_name(kr)
                if valid_variable_name(k):
                    vtype = str(type(v))
                    if vtype in HDFreader.datasets:
                        self.add_dataset(cur,k,v)
                    elif vtype in HDFreader.groups:
                        self.add_group(hcur,cur,k,v)
                    else:
                        print('hdfreader error: encountered invalid type: '+vtype)
                        sys.exit()
                    #end if
                else:
                    print('hdfreader warning: attribute '+k+' is not a valid variable name and has been ignored')
                #end if
            #end for
        #end if

        if verbose:
            print('  end HDFreader Initialization')

        return
    #end def __init__


    def increment_level(self):
        self.ilevel+=1
        self.nlevels = max(self.ilevel+1,self.nlevels)
        if self.ilevel+1==self.nlevels:
            self.cur.append(None)
            self.hcur.append(None)
        #end if
        self.pad = self.ilevel*'  '
        return
    #end def increment_level

    def decrement_level(self):
        self.ilevel-=1
        self.pad = self.ilevel*'  '
        return
    #end def decrement_level

    def add_dataset(self,cur,k,v):
        if not HDFglobals.view:
            cur[k]=array(v)
        else:
            cur[k] = v
        #end if
        cur._add_dataset(k,cur[k])
        return
    #end def add_dataset

    def add_group(self,hcur,cur,k,v):
        cur[k] = HDFgroup()
        cur._add_group(k,cur[k])
        cur._groups[k]._parent = cur
        self.increment_level()
        self.cur[self.ilevel]  = cur._groups[k]
        self.hcur[self.ilevel] = hcur[k]

        cur   = self.cur[self.ilevel]
        hcur  = self.hcur[self.ilevel]
        for kr,v in hcur.items():
            k=cur._escape_name(kr)
            if valid_variable_name(k):
                vtype = str(type(v))
                if vtype in HDFreader.datasets:
                    self.add_dataset(cur,k,v)
                elif vtype in HDFreader.groups:
                    self.add_group(hcur,cur,k,v)
                #end if
            else:
                print('hdfreader warning: attribute '+k+' is not a valid variable name and has been ignored')
            #end if
        #end for

        return
    #end def add_group
#end class HDFreader



def read_hdf(fpath,verbose=False,view=False):
    return HDFreader(fpath=fpath,verbose=verbose,view=view).obj
#end def read_hdf

######################################################################
# end from hdfreader.py
######################################################################









import os
import sys
from optparse import OptionParser
from numpy import zeros,sqrt,longdouble,loadtxt
can_plot = False
try:
    import matplotlib
    gui_envs = ['GTKAgg','TKAgg','Qt4Agg','WXAgg']
    for gui in gui_envs:
        try:
            matplotlib.use(gui,warn=False, force=True)
            from matplotlib import pyplot
            can_plot = True
            break
        except:
            continue
        #end try
    #end for
    from matplotlib.pyplot import figure,plot,xlabel,ylabel,title,show,ylim,legend,xlim,rcParams,savefig,bar,xticks,subplot,grid,setp,errorbar,loglog,semilogx,semilogy,text

    params = {'legend.fontsize':14,'figure.facecolor':'white','figure.subplot.hspace':0.,
          'axes.labelsize':16,'xtick.labelsize':14,'ytick.labelsize':14}
    rcParams.update(params)
except (ImportError,RuntimeError):
    can_plot = False
#end try



class ColorWheel(DevBase):
    def __init__(self):
        colors  = 'Black Maroon DarkOrange Green DarkBlue Purple Gray Firebrick Orange MediumSeaGreen DodgerBlue MediumOrchid'.split()
        lines  = '- -- -. :'.split()
        markers = '. v s o ^ d p'.split() 
        ls = []
        for line in lines:
            for color in colors:
                ls.append((color,line))
            #end for
        #end for
        ms = []
        for i in range(len(markers)):
            ms.append((colors[i],markers[i]))
        #end for
        mls = []
        ic=-1
        for line in lines:
            for marker in markers:
                ic = (ic+1)%len(colors)
                mls.append((colors[ic],marker+line))
            #end for
        #end for
        self.line_styles   = ls
        self.marker_styles = ms
        self.marker_line_styles = mls
        self.reset()
    #end def __init__

    def next_line(self):
        self.iline = (self.iline+1)%len(self.line_styles)
        return self.line_styles[self.iline]
    #end def next_line

    def next_marker(self):
        self.imarker = (self.imarker+1)%len(self.marker_styles)
        return self.marker_styles[self.imarker]
    #end def next_marker

    def next_marker_line(self):
        self.imarker_line = (self.imarker_line+1)%len(self.marker_line_styles)
        return self.marker_line_styles[self.imarker_line]
    #end def next_marker_line

    def reset(self):
        self.iline        = -1
        self.imarker      = -1
        self.imarker_line = -1
    #end def reset

    def reset_line(self):
        self.iline        = -1
    #end def reset_line

    def reset_marker(self):
        self.imarker      = -1
    #end def reset_marker

    def reset_marker_line(self):
        self.imarker_line = -1
    #end def reset_marker_line
#end class ColorWheel
color_wheel = ColorWheel()



checkstats_settings = obj(
    verbose = False,
    )

def vlog(*args,**kwargs):
    if checkstats_settings.verbose:
        n = kwargs.get('n',0)
        if n==0:
            log(*args)
        else:
            log(*args,indent=n*'  ')
        #end if
    #end if
#end def vlog





# standalone definition of error function from Abramowitz & Stegun
# credit: http://www.johndcook.com/blog/2009/01/19/stand-alone-error-function-erf/
# consider also: https://math.stackexchange.com/questions/42920/efficient-and-accurate-approximation-of-error-function
def erf(x):
    # constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911

    # Save the sign of x
    sign = 1
    if x < 0:
        sign = -1
    x = abs(x)

    # A & S 7.1.26
    t = 1.0/(1.0 + p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)

    return sign*y
#end def erf


# standalone inverse error function 
# credit: https://stackoverflow.com/questions/42381244/pure-python-inverse-error-function
import math
def polevl(x, coefs, N):
    ans = 0
    power = len(coefs) - 1
    for coef in coefs:
        ans += coef * x**power
        power -= 1
    return ans
#end def polevl

def p1evl(x, coefs, N):
    return polevl(x, [1] + coefs, N)
#end def p1evl

def erfinv(z):
    if z < -1 or z > 1:
        raise ValueError("'z' must be between -1 and 1 inclusive")

    if z == 0:
        return 0
    if z == 1:
        return float('inf')
    if z == -1:
        return -float('inf')

    # From scipy special/cephes/ndrti.c
    def ndtri(y):
        # approximation for 0 <= abs(z - 0.5) <= 3/8
        P0 = [
            -5.99633501014107895267E1,
             9.80010754185999661536E1,
             -5.66762857469070293439E1,
             1.39312609387279679503E1,
             -1.23916583867381258016E0,
             ]

        Q0 = [
            1.95448858338141759834E0,
            4.67627912898881538453E0,
            8.63602421390890590575E1,
            -2.25462687854119370527E2,
            2.00260212380060660359E2,
            -8.20372256168333339912E1,
            1.59056225126211695515E1,
            -1.18331621121330003142E0,
            ]

        # Approximation for interval z = sqrt(-2 log y ) between 2 and 8
        # i.e., y between exp(-2) = .135 and exp(-32) = 1.27e-14.
        P1 = [
            4.05544892305962419923E0,
            3.15251094599893866154E1,
            5.71628192246421288162E1,
            4.40805073893200834700E1,
            1.46849561928858024014E1,
            2.18663306850790267539E0,
            -1.40256079171354495875E-1,
            -3.50424626827848203418E-2,
            -8.57456785154685413611E-4,
            ]

        Q1 = [
            1.57799883256466749731E1,
            4.53907635128879210584E1,
            4.13172038254672030440E1,
            1.50425385692907503408E1,
            2.50464946208309415979E0,
            -1.42182922854787788574E-1,
            -3.80806407691578277194E-2,
            -9.33259480895457427372E-4,
            ]

        # Approximation for interval z = sqrt(-2 log y ) between 8 and 64
        # i.e., y between exp(-32) = 1.27e-14 and exp(-2048) = 3.67e-890.
        P2 = [
            3.23774891776946035970E0,
            6.91522889068984211695E0,
            3.93881025292474443415E0,
            1.33303460815807542389E0,
            2.01485389549179081538E-1,
            1.23716634817820021358E-2,
            3.01581553508235416007E-4,
            2.65806974686737550832E-6,
            6.23974539184983293730E-9,
            ]

        Q2 = [
            6.02427039364742014255E0,
            3.67983563856160859403E0,
            1.37702099489081330271E0,
            2.16236993594496635890E-1,
            1.34204006088543189037E-2,
            3.28014464682127739104E-4,
            2.89247864745380683936E-6,
            6.79019408009981274425E-9,
            ]

        s2pi = 2.50662827463100050242
        code = 1

        if y > (1.0 - 0.13533528323661269189):      # 0.135... = exp(-2)
            y = 1.0 - y
            code = 0

        if y > 0.13533528323661269189:
            y = y - 0.5
            y2 = y * y
            x = y + y * (y2 * polevl(y2, P0, 4) / p1evl(y2, Q0, 8))
            x = x * s2pi
            return x

        x = math.sqrt(-2.0 * math.log(y))
        x0 = x - math.log(x) / x

        z = 1.0 / x
        if x < 8.0:                 # y > exp(-32) = 1.2664165549e-14
            x1 = z * polevl(z, P1, 8) / p1evl(z, Q1, 8)
        else:
            x1 = z * polevl(z, P2, 8) / p1evl(z, Q2, 8)

        x = x0 - x1
        if code != 0:
            x = -x

        return x

    result = ndtri((z + 1) / 2.0) / math.sqrt(2)

    return result
#end def inv_erf





# Returns failure error code to OS.
# Explicitly prints 'fail' after an optional message.
def exit_fail(msg=None):
    if msg!=None:
        print(msg)
    #end if
    print('Test status: fail')
    exit(1)
#end def exit_fail


# Returns success error code to OS.
# Explicitly prints 'pass' after an optional message.
def exit_pass(msg=None):
    if msg!=None:
        print(msg)
    #end if
    print('Test status: pass')
    exit(0)
#end def exit_pass


# Calculates the mean, variance, errorbar, and autocorrelation time
# for a N-d array of statistical data values.
# If 'exclude' is provided, the first 'exclude' values will be 
# excluded from the analysis.
def simstats(x,dim=None,exclude=None):
    if exclude!=None:
        x = x[exclude:]
    #end if
    shape = x.shape
    ndim  = len(shape)
    if dim==None:
        dim=ndim-1
    #end if
    permute = dim!=ndim-1
    reshape = ndim>2
    nblocks = shape[dim]
    if permute:
        r = list(range(ndim))
        r.pop(dim)
        r.append(dim)
        permutation = tuple(r)
        r = list(range(ndim))
        r.pop(ndim-1)
        r.insert(dim,ndim-1)
        invperm     = tuple(r)
        x=x.transpose(permutation)
        shape = tuple(array(shape)[array(permutation)])
        dim = ndim-1
    #end if
    if reshape:        
        nvars = prod(shape[0:dim])
        x=x.reshape(nvars,nblocks)
        rdim=dim
        dim=1
    else:
        nvars = shape[0]
    #end if

    mean  = x.mean(dim)
    var   = x.var(dim)

    N=nblocks

    if ndim==1:
        i=0          
        tempC=0.5
        kappa=0.0
        mtmp=mean
        if abs(var)<1e-15:
            kappa = 1.0
        else:
            ovar=1.0/var
            while (tempC>0 and i<(N-1)):
                kappa=kappa+2.0*tempC
                i=i+1
                #tempC=corr(i,x,mean,var)
                tempC = ovar/(N-i)*sum((x[0:N-i]-mtmp)*(x[i:N]-mtmp))
            #end while
            if kappa == 0.0:
                kappa = 1.0
            #end if
        #end if
        Neff=(N+0.0)/(kappa+0.0)
        if (Neff == 0.0):
            Neff = 1.0
        #end if
        error=sqrt(var/Neff)
    else:
        error = zeros(mean.shape)
        kappa = zeros(mean.shape)
        for v in range(nvars):
            i=0          
            tempC=0.5
            kap=0.0
            vtmp = var[v]
            mtmp = mean[v]
            if abs(vtmp)<1e-15:
                kap = 1.0
            else:
                ovar   = 1.0/vtmp
                while (tempC>0 and i<(N-1)):
                    i += 1
                    kap += 2.0*tempC
                    tempC = ovar/(N-i)*sum((x[v,0:N-i]-mtmp)*(x[v,i:N]-mtmp))
                #end while
                if kap == 0.0:
                    kap = 1.0
                #end if
            #end if
            Neff=(N+0.0)/(kap+0.0)
            if (Neff == 0.0):
                Neff = 1.0
            #end if
            kappa[v]=kap
            error[v]=sqrt(vtmp/Neff)
        #end for    
    #end if

    if reshape:
        x     =     x.reshape(shape)
        mean  =  mean.reshape(shape[0:rdim])
        var   =   var.reshape(shape[0:rdim])
        error = error.reshape(shape[0:rdim])
        kappa = kappa.reshape(shape[0:rdim])
    #end if
    if permute:
        x=x.transpose(invperm)
    #end if

    return (mean,var,error,kappa)
#end def simstats



def load_scalar_file(options,selector):
    output_files = options.output_files
    if selector=='auto':
        if 'dmc' in output_files:
            selector = 'dmc'
        elif 'scalar' in output_files:
            selector = 'scalar'
        else:
            exit_fail('could not load scalar file, no files present')
        #end if
    elif selector not in ('scalar','dmc'):
        exit_fail('could not load scalar file, invalid selector\ninvalid selector: {0}\nvalid options: scalar, dmc'.format(selector))
    #end if
    if selector not in output_files:
        exit_fail('could not load scalar file, file is not present\nfile type requested: {0}'.format(selector))
    #end if
    filepath = os.path.join(options.path,output_files[selector])
    lt = loadtxt(filepath)
    if len(lt.shape)==1:
        lt.shape = (1,len(lt))
    #end if
    data = lt[:,1:].transpose()
    fobj = open(filepath,'r')
    variables = fobj.readline().split()[2:]
    fobj.close()
    scalars = obj(
        file_type = selector,
        data      = obj(),
        )
    for i,var in enumerate(variables):
        scalars.data[var] = data[i,:]
    #end for
    return scalars
#end def load_scalar_file



# Reads command line options.
def read_command_line():
    try:

        parser = OptionParser(
            usage='usage: %prog [options]',
            add_help_option=False,
            version='%prog 0.1'
            )

        parser.add_option('-h','--help',dest='help',
                          action='store_true',default=False,
                          help='Print help information and exit (default=%default).'
                          )
        parser.add_option('-p','--prefix',dest='prefix',
                          default='qmc',
                          help='Prefix for output files (default=%default).  Can be a path including the file prefix.'
                          )
        parser.add_option('-s','--series',dest='series',
                          default='0',
                          help='Output series to analyze (default=%default).'
                          )
        parser.add_option('-e','--equilibration',dest='equilibration',
                          default='0',
                          help='Equilibration length in blocks (default=%default).'
                          )
        parser.add_option('-n','--nsigma',dest='nsigma',
                          default='3',
                          help='Sigma requirement for pass/fail (default=%default).'
                          )
        parser.add_option('-q','--quantity',dest='quantity',
                          default='none',
                          help = 'Quantity to check (required).  If a non-default name for the quantity is used, pass in the quantity and name as a pair.'
                          )
        parser.add_option('-c','--npartial_sums',dest='npartial_sums',
                          default='none',
                          help = 'Partial sum count for the reference data (required)'
                          )
        parser.add_option('-r','--ref','--reference',dest='reference_file',
                          default='none',
                          help = 'Path to reference file containing full and partial sum reference information.  The test fails if any full or partial sum exceeds nsigma deviation from the reference values.  For cases like the density or spin density, the -f option should additionally be used (see below).  For the energy density, a block by block check against relevant summed energy terms in scalar.dat or dmc.dat files is additionally made.'
                          )
        parser.add_option('-f','--fixed','--fixed_sum',dest='fixed_sum',
                          action='store_true',default=False,
                          help = 'Full sum of data takes on a fixed, non-stochastic value.  In this case, when checking against reference data, check that each block satisfies the fixed sum condition.  This is appropriate, e.g. for the electron density where the full sum of each block must equal the number of electrons.  Typically the appropriate value is inferred automatically and applied by default (in other cases default=%default).'
                          )
        parser.add_option('-m','--make_ref','--make_reference',dest='make_reference',
                          default='none',
                          help='Used during test construction phase.  Pass an integer list via -m corresponding to the number of partial sums to perform on the reference stat data followed by a series of MC step factors.  The number of partial means must divide evenly into the number of stat field values for the quantity in question.  The step factors relate the length of the test run (shorter) to the reference run (longer): #MC_test*factor=#MC_reference.  Files containing the reference data will be produced, one for each step factor.  For the partial sums, the reference sigma is increased so that the test fails with the expected probability specified by the inputted nsigma.'
                          )
        parser.add_option('-t','--plot_trace',dest='plot_trace',
                          action='store_true',default=False,
                          help='Plot traces of full and partial sums (default=%default).'
                          )
        parser.add_option('-v','--verbose',
                          action='store_true',default=False,
                          help='Print detailed information (default=%default).'
                          )

        allowed_quantities = [
            'density',
            'spindensity',
            'energydensity',
            '1rdm',
            '1redm',
            'obdm',
            'momentum',
            ]

        opt,files_in = parser.parse_args()
        options = obj()
        options.transfer_from(opt.__dict__)

        if options.help:
            print('\n'+parser.format_help().strip())
            print('\n\nExample usage:')
            print('\n  Making reference data to create a test:')
            print("    check_stats.py -p qmc -s 0 -q spindensity -e 10 -c 8 -v -m '0 10 100'")
            print('\n  Using reference data to perform a test:')
            print('    check_stats.py -p qmc -s 0 -q spindensity -e 10 -c 8 -n 3 -r qmc.s000.stat_ref_spindensity_10.dat')
            print()
            exit()
        #end if

        if len(files_in)>0:
            exit_fail('check_stats does not accept file as input, only command line arguments\nfiles provided: {0}'.format(files_in))
        #end if

        checkstats_settings.verbose = options.verbose

        vlog('\nreading command line inputs')

        options.series        = int(options.series)
        options.equilibration = int(options.equilibration)
        options.nsigma        = float(options.nsigma)
        options.path,options.prefix = os.path.split(options.prefix)

        if options.plot_trace and not can_plot:
            vlog('trace plots requested, but plotting libraries are not available\ndisabling plots',n=1)
            options.plot_trace = False
        #end if

        if options.path=='':
            options.path = './'
        #end if

        options.qlabel = None
        if ' ' in options.quantity or ',' in options.quantity:
            qlist = options.quantity.strip('"').strip("'").replace(',',' ').split()
            if len(qlist)!=2:
                exit_fail('quantity can accept only one or two values\nyou provided {0}: {1}'.format(len(qlist),qlist))
            #end if
            options.quantity,options.qlabel = qlist
        #end if
        if options.qlabel is None:
            default_label = obj({
                'density'       : 'Density'        ,
                'spindensity'   : 'SpinDensity'    ,
                'energydensity' : 'EnergyDensity'  ,
                '1rdm'          : 'DensityMatrices',
                '1redm'         : 'DensityMatrices',
                'obdm'          : 'OneBodyDensityMatrices' ,
                'momentum'      : 'nofk'           ,
                })
            options.qlabel = default_label[options.quantity]
        #end if
        if options.quantity=='none':
            exit_fail('must provide quantity')
        elif options.quantity not in allowed_quantities:
            exit_fail('unrecognized quantity provided\nallowed quantities: {0}\nquantity provided: {1}'.format(allowed_quantities,options.quantity))
        #end if
            
        if options.npartial_sums=='none':
            exit_fail('-c option is required')
        #end if
        options.npartial_sums = int(options.npartial_sums)

        if options.reference_file!='none':
            if not os.path.exists(options.reference_file):
                exit_fail('reference file does not exist\nreference file provided: {0}'.format(options.reference_file))
            #end if
            options.make_reference = False
        elif options.make_reference!='none':
            try:
                mr = array(options.make_reference.split(),dtype=int)
            except:
                exit_fail('make_reference must be a list of integers\nyou provided: {0}'.format(options.make_reference))
            #end try
            if len(mr)<1:
                exit_fail('make_reference must contain at least one MC length factor')
            #end if
            options.mc_factors     = mr
            options.make_reference = True
        else:
            exit_fail('must provide either reference_file or make_reference')
        #end if

        fixed_sum_quants = set(['density','spindensity','energydensity'])
        if options.quantity in fixed_sum_quants:
            options.fixed_sum = True
        #end if

        vlog('inputted options:\n'+str(options),n=1)

    except Exception as e:
        exit_fail('error during command line read:\n'+str(e))
    #end try

    return options
#end def read_command_line




def process_stat_file(options):
    vlog('processing stat.h5 file')

    values = obj()
    
    try:
        # find all output files matching prefix
        vlog('searching for qmcpack output files',n=1)
        vlog('search path:\n  '+options.path,n=2)
        prefix = options.prefix+'.s'+str(options.series).zfill(3)
        files = os.listdir(options.path)
        output_files = obj()
        for file in files:
            if file.startswith(prefix):
                if file.endswith('.stat.h5'):
                    output_files.stat = file
                elif file.endswith('.scalar.dat'):
                    output_files.scalar = file
                elif file.endswith('.dmc.dat'):
                    output_files.dmc = file
                #end if
            #end if
        #end for
        options.output_files = output_files
        vlog('files found:\n'+str(output_files).rstrip(),n=2)
        if 'stat' not in output_files:
            exit_fail('stat.h5 file matching prefix {0} was not found\nsearch path: {1}'.format(prefix,options.path))
        #end if

        # read data from the stat file
        vlog('opening stat.h5 file',n=1)
        stat = read_hdf(os.path.join(options.path,output_files.stat),view=True)
        vlog('file contents:\n'+repr(stat).rstrip(),n=2)
        vlog('extracting {0} data'.format(options.quantity),n=1)
        vlog('searching for {0} with label {1}'.format(options.quantity,options.qlabel),n=2)
        if options.qlabel in stat:
            qstat = stat[options.qlabel]
            vlog('{0} data contents:\n{1}'.format(options.quantity,repr(qstat).rstrip()),n=2)
        else:
            exit_fail('could not find {0} data with label {1}'.format(options.quantity,options.qlabel))
        #end if
        quantity_paths = obj({
            'density'       : obj(tot='value'),
            'spindensity'   : obj(u='u/value',
                                  d='d/value'),
            '1rdm'          : obj(u='number_matrix/u/value',
                                  d='number_matrix/d/value'),
            '1redm'         : obj(u='energy_matrix/u/value',
                                  d='energy_matrix/d/value'),
            'obdm'          : obj(u='number_matrix/u/value',
                                  d='number_matrix/d/value'),
            'energydensity' : obj(W=('spacegrid1/value',0,3),
                                  T=('spacegrid1/value',1,3),
                                  V=('spacegrid1/value',2,3)),
            'momentum'      : obj(tot='value'),
            })
        qpaths = quantity_paths[options.quantity]
        vlog('search paths:\n{0}'.format(str(qpaths).rstrip()),n=2)
        qdata = obj()
        dfull = None
        for dname,dpath in qpaths.items():
            packed = isinstance(dpath,tuple)
            if packed:
                dpath,dindex,dcount = dpath
            #end if
            if not qstat.path_exists(dpath):
                exit_fail('{0} data not found in file {1}\npath searched: {2}'.format(options.quantity,output_files.stat,dpath))
            #end if
            if not packed:
                d = array(qstat.get_path(dpath),dtype=float)
            else:
                if dfull is None:
                    dfull = array(qstat.get_path(dpath),dtype=float)
                    dfull.shape = dfull.shape[0],dfull.shape[1]//dcount,dcount
                #end if
                d = dfull[:,:,dindex]
                d.shape = dfull.shape[0],dfull.shape[1]
            #end if
            qdata[dname] = d
            vlog('{0} data found with shape {1}'.format(dname,d.shape),n=2)
            if len(d.shape)>2:
                d.shape = d.shape[0],d.size//d.shape[0]
                vlog('reshaped {0} data to {1}'.format(dname,d.shape),n=2)
            #end if
            options.nblocks = d.shape[0]
        #end for

        # process the data, taking full and partial sums
        vlog('processing {0} data'.format(options.quantity),n=1)
        for dname,d in qdata.items():
            vlog('processing {0} data'.format(dname),n=2)
            if d.shape[1]%options.npartial_sums!=0:
                exit_fail('cannot make partial sums\nnumber of requested partial sums does not divide evenly into the number of values available\nrequested partial sums: {0}\nnumber of values present: {1}\nnvalue/npartial_sums: {2}'.format(options.npartial_sums,d.shape[1],float(d.shape[1])/options.npartial_sums))
            #end if
            data = obj()
            data.full_sum = d.sum(1)
            vlog('full sum data shape: {0}'.format(data.full_sum.shape),n=3)
            data.partial_sums = zeros((d.shape[0],options.npartial_sums))
            psize = d.shape[1]//options.npartial_sums
            for p in range(options.npartial_sums):
                data.partial_sums[:,p] = d[:,p*psize:(p+1)*psize].sum(1)
            #end for
            vlog('partial sum data shape: {0}'.format(data.partial_sums.shape),n=3)
            fmean,var,ferror,kappa = simstats(data.full_sum,exclude=options.equilibration)
            vlog('full sum mean : {0}'.format(fmean),n=3)
            vlog('full sum error: {0}'.format(ferror),n=3)
            pmean,var,perror,kappa = simstats(data.partial_sums,dim=0,exclude=options.equilibration)
            vlog('partial sum mean : {0}'.format(pmean),n=3)
            vlog('partial sum error: {0}'.format(perror),n=3)
            values[dname] = obj(
                full_mean     = fmean,
                full_error    = ferror,
                partial_mean  = pmean,
                partial_error = perror,
                data          = data,
                )
        #end for

        # check that all values have been processed
        missing = set(qpaths.keys())-set(values.keys())
        if len(missing)>0:
            exit_fail('some values not processed\nvalues missing: {0}'.format(sorted(missing)))
        #end if

        # plot quantity traces, if requested
        if options.plot_trace:
            vlog('creating trace plots of full and partial sums',n=1)
            for dname,dvalues in values.items():
                label = options.quantity
                if len(values)>1:
                    label+=' '+dname
                #end if
                data = dvalues.data
                figure()
                plot(data.full_sum)
                title('Trace of {0} full sum'.format(label))
                xlabel('Block index')
                figure()
                plot(data.partial_sums)
                title('Trace of {0} partial sums'.format(label))
                xlabel('Block index')
            #end for
            show()
        #end if
    except Exception as e:
        exit_fail('error during stat file processing:\n'+str(e))
    #end try

    return values
#end def process_stat_file



def make_reference_files(options,values):
    vlog('\nmaking reference files')

    # create a reference file for each Monte Carlo sample factor
    for mcfac in options.mc_factors:
        errfac = sqrt(1.0+mcfac)
        filename = '{0}.s{1}.stat_ref_{2}_{3}.dat'.format(options.prefix,str(options.series).zfill(3),options.quantity,mcfac)
        filepath = os.path.join(options.path,filename)
        vlog('writing reference file for {0}x shorter test runs'.format(mcfac),n=1)
        vlog('reference file location: '+filepath,n=2)
        f = open(filepath,'w')
        # write descriptive header line
        line = '# '
        for dname in sorted(values.keys()):
            line += '  {0:<16}  {1:<16}'.format(dname,dname+'_err')
        #end for
        f.write(line+'\n')
        # write means and errors of full sum
        line = ''
        for dname in sorted(values.keys()):
            dvalues = values[dname]
            fmean   = dvalues.full_mean
            ferror  = dvalues.full_error
            line += '  {0: 16.8e}  {1: 16.8e}'.format(fmean,errfac*ferror)
        #end for
        f.write(line+'\n')
        # write means and errors of partial sums
        for p in range(options.npartial_sums):
            line = ''
            for dname in sorted(values.keys()):
                dvalues = values[dname]
                pmean   = dvalues.partial_mean
                perror  = dvalues.partial_error
                line += '  {0: 16.8e}  {1: 16.8e}'.format(pmean[p],errfac*perror[p])
            #end for
            f.write(line+'\n')
        #end for
        f.close()
    #end for

    # create a trace file containing full and partial sum data per block
    filename = '{0}.s{1}.stat_trace_{2}.dat'.format(options.prefix,str(options.series).zfill(3),options.quantity)
    filepath = os.path.join(options.path,filename)
    vlog('writing trace file containing full and partial sums per block',n=1)
    vlog('trace file location: '+filepath,n=2)
    f = open(filepath,'w')
    # write descriptive header line
    line = '# index  '
    for dname in sorted(values.keys()):
        line += '  {0:<16}'.format(dname+'_full')
        for p in range(options.npartial_sums):
            line += '  {0:<16}'.format(dname+'_partial_'+str(p))
        #end for
    #end for
    f.write(line+'\n')
    # write full and partial sum data per block
    for b in range(options.nblocks):
        line = ' {0:>6}'.format(b)
        for dname in sorted(values.keys()):
            dvalues = values[dname].data
            fsum  = dvalues.full_sum
            psums = dvalues.partial_sums[b]
            line += '  {0: 16.8e}'.format(fsum[b])
            for psum in psums:
                line += '  {0: 16.8e}'.format(psum)
            #end for
        #end for
        f.write(line+'\n')
    #end for
    f.close()
    vlog('\n')
#end def make_reference_files



# Checks computed values from stat.h5 files against specified reference values.
passfail = {True:'pass',False:'fail'}
def check_values(options,values):
    vlog('\nchecking against reference values')

    success = True
    msg     = ''

    try:
        msg += '\nTests for series {0} quantity "{1}"\n'.format(options.series,options.quantity)

        # find nsigma for each partial sum
        #   overall probability of partial sum failure is according to original nsigma
        vlog('adjusting nsigma to account for partial sum count',n=1)
        x = longdouble(options.nsigma/sqrt(2.))
        N = options.npartial_sums
        nsigma_partial = sqrt(2.)*erfinv(erf(x)**(1./N))
        vlog('overall full/partial test nsigma: {0}'.format(options.nsigma),n=2)
        vlog('adjusted per partial sum nsigma : {0}'.format(nsigma_partial),n=2)

        # read in the reference file
        vlog('reading reference file',n=1)
        vlog('reference file location: {0}'.format(options.reference_file),n=2)
        f = open(options.reference_file,'r')
        dnames = f.readline().split()[1::2]
        vlog('sub-quantities found: {0}'.format(dnames),n=2)
        if set(dnames)!=set(values.keys()):
            missing = set(values.keys())-set(dnames)
            extra   = set(dnames)-set(values.keys())
            if missing>0:
                exit_fail('some sub-quantities are missing\npresent in test files: {0}\npresent in reference files: {1}\nmissing: {2}'.format(sorted(values.keys()),sorted(dnames),sorted(missing)))
            elif extra>0:
                exit_fail('some sub-quantities are extra\npresent in test files: {0}\npresent in reference files: {1}\nextra: {2}'.format(sorted(values.keys()),sorted(dnames),sorted(extra)))
            else:
                exit_fail('developer error, this point should be impossible to reach')
            #end if
        #end if
        ref = array(f.read().split(),dtype=float)
        ref.shape = len(ref)//(2*len(dnames)),2*len(dnames)
        full    = ref[0,:].ravel()
        partial = ref[1:,:].T
        if len(ref)-1!=options.npartial_sums:
            exit_fail('test and reference partial sum counts do not match\ntest partial sum count: {0}\nreference partial sum count: {1}'.format(options.npartial_sums,len(ref)-1))
        #end if
        vlog('partial sum count found: {0}'.format(len(ref)-1),n=2)
        ref_values = obj()
        for dname in dnames:
            ref_values[dname] = obj()
        #end for
        n=0
        for dname in dnames:
            ref_values[dname].set(
                full_mean     = full[n],
                full_error    = full[n+1],
                partial_mean  = partial[n,:].ravel(),
                partial_error = partial[n+1,:].ravel(),
                )
            n+=2
        #end for
        npartial = len(ref)-1
        f.close()
        dnames = sorted(dnames)
        vlog('reference file read successfully',n=2)

        # for cases with fixed full sum, check the per block sum
        if options.fixed_sum:
            vlog('checking per block fixed sums',n=1)
            msg+='\n  Fixed sum per block tests:\n'
            dnames_fixed = dnames
            if options.quantity=='energydensity':
                dnames_fixed = ['W']
            #end if
            fixed_sum_success = True
            ftol = 1e-8
            eq = options.equilibration
            for dname in dnames_fixed:
                ref_vals  = ref_values[dname]
                ref_mean  = ref_vals.full_mean
                ref_error = ref_vals.full_error
                if abs(ref_error/ref_mean)>ftol:
                    exit_fail('reference fixed sum is not fixed as asserted\ncannot check per block fixed sums\nplease check reference data')
                #end if
                test_vals = values[dname].data.full_sum
                for i,v in enumerate(test_vals[eq:]):
                    if abs((v-ref_mean)/ref_mean)>ftol:
                        fixed_sum_success = False
                        msg += '    {0} {1} {2}!={3}\n'.format(dname,i,v,ref_mean)
                    #end if
                #end for
            #end for
            if fixed_sum_success:
                fmsg = 'all per block sums match the reference'
            else:
                fmsg = 'some per block sums do not match the reference'
            #end if
            vlog(fmsg,n=2)
            msg += '    '+fmsg+'\n'
            msg += '    status of this test: {0}\n'.format(passfail[fixed_sum_success])
            success &= fixed_sum_success
        #end if

        # for the energy density, check per block sums against the scalar file
        if options.quantity=='energydensity':
            vlog('checking energy density terms per block',n=1)
            msg+='\n  Energy density sums vs. scalar file per block tests:\n'
            scalars = load_scalar_file(options,'auto')
            ed_success = True
            ftol = 1e-8
            ed_values = obj(
                T = values.T.data.full_sum,
                V = values.V.data.full_sum,
                )
            ed_values.E = ed_values.T + ed_values.V
            if scalars.file_type=='scalar':
                comparisons = obj(
                    E='LocalEnergy',
                    T='Kinetic',
                    V='LocalPotential',
                    )
            elif scalars.file_type=='dmc':
                comparisons = obj(E='LocalEnergy')
            else:
                exit_fail('unrecognized scalar file type: {0}'.format(scalars.file_type))
            #end if
            for k in sorted(comparisons.keys()):
                ed_vals = ed_values[k]
                sc_vals = scalars.data[comparisons[k]]
                if scalars.file_type=='dmc':
                    if len(sc_vals)%len(ed_vals)==0 and len(sc_vals)>len(ed_vals):
                        steps = len(sc_vals)//len(ed_vals)
                        sc_vals.shape = len(sc_vals)//steps,steps
                        sc_vals = sc_vals.mean(1)
                    #end if
                #end if
                if len(ed_vals)!=len(sc_vals):
                    exit_fail('energy density per block test cannot be completed\nnumber of energy density and scalar blocks do not match\nenergy density blocks: {0}\nscalar file blocks: {1}'.format(len(ed_vals),len(sc_vals)))
                #end if
                for i,(edv,scv) in enumerate(zip(ed_vals,sc_vals)):
                    if abs((edv-scv)/scv)>ftol:
                        ed_success = False
                        msg += '    {0} {1} {2}!={3}\n'.format(k,i,edv,scv)
                    #end if
                #end for
            #end for
            if ed_success:
                fmsg = 'all per block sums match the scalar file'
            else:
                fmsg = 'some per block sums do not match the scalar file'
            #end if
            vlog(fmsg,n=2)
            msg += '    '+fmsg+'\n'
            msg += '    status of this test: {0}\n'.format(passfail[ed_success])
            success &= ed_success
        #end if


        # function used immediately below to test a mean value vs reference
        def check_mean(label,mean_comp,error_comp,mean_ref,error_ref,nsigma):
            msg='\n  Testing quantity: {0}\n'.format(label)

            # ensure error_ref is large enough for non-statistical quantities
            ctol = 1e-12
            if abs(error_ref/mean_ref)<ctol:
                error_ref = ctol*mean_ref
            #end if

            quant_success = abs(mean_comp-mean_ref) <= nsigma*error_ref

            delta = mean_comp-mean_ref
            delta_err = sqrt(error_comp**2+error_ref**2)

            msg+='    reference mean value     : {0: 12.8f}\n'.format(mean_ref)
            msg+='    reference error bar      : {0: 12.8f}\n'.format(error_ref)
            msg+='    computed  mean value     : {0: 12.8f}\n'.format(mean_comp)
            msg+='    computed  error bar      : {0: 12.8f}\n'.format(error_comp)
            msg+='    pass tolerance           : {0: 12.8f}  ({1: 12.8f} sigma)\n'.format(nsigma*error_ref,nsigma)
            if error_ref > 0.0:
                msg+='    deviation from reference : {0: 12.8f}  ({1: 12.8f} sigma)\n'.format(delta,delta/error_ref)
            #end if
            msg+='    error bar of deviation   : {0: 12.8f}\n'.format(delta_err)
            if error_ref > 0.0:
                msg+='    significance probability : {0: 12.8f}  (gaussian statistics)\n'.format(erf(abs(delta/error_ref)/math.sqrt(2.0)))
            #end if
            msg+='    status of this test      :   {0}\n'.format(passfail[quant_success])

            return quant_success,msg
        #end def check_mean


        # check full and partial sums vs the reference
        vlog('checking full and partial sums',n=1)
        for dname in dnames:
            vals = values[dname]
            ref_vals = ref_values[dname]
            # check full sum
            vlog('checking full sum mean for "{0}"'.format(dname),n=2)
            qsuccess,qmsg = check_mean(
                label      = '{0} full sum'.format(dname),
                mean_comp  = vals.full_mean,
                error_comp = vals.full_error,
                mean_ref   = ref_vals.full_mean,
                error_ref  = ref_vals.full_error,
                nsigma     = options.nsigma,
                )
            vlog('status for full sum: {0}'.format(passfail[qsuccess]),n=3)
            msg     += qmsg
            success &= qsuccess
            # check partial sums
            vlog('checking partial sum means for "{0}"'.format(dname),n=2)
            for p in range(len(ref_vals.partial_mean)):
                qsuccess,qmsg = check_mean(
                    label      = '{0} partial sum {1}'.format(dname,p),
                    mean_comp  = vals.partial_mean[p],
                    error_comp = vals.partial_error[p],
                    mean_ref   = ref_vals.partial_mean[p],
                    error_ref  = ref_vals.partial_error[p],
                    nsigma     = nsigma_partial,
                    )
                vlog('status for partial sum {0}: {1}'.format(p,passfail[qsuccess]),n=3)
                msg     += qmsg
                success &= qsuccess
            #end for
        #end for
    except Exception as e:
        exit_fail('error during value check:\n'+str(e))
    #end try

    return success,msg
#end def check_values




# Main execution
if __name__=='__main__':
    # Read and interpret command line options.
    options = read_command_line()

    # Compute means of desired quantities from a stat.h5 file.
    values = process_stat_file(options)

    if options.make_reference:
        # Make reference files to create tests from.
        make_reference_files(options,values)
        exit()
    else:
        # Check computed means against reference solutions.
        success,msg = check_values(options,values)

        # Pass success/failure exit codes and strings to the OS.
        if success:
            exit_pass(msg)
        else:
            exit_fail(msg)
        #end if
    #end if
#end if
