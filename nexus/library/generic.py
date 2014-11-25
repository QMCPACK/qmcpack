
import sys
import traceback
from copy import deepcopy
from abilities import AllAbilities,genbase

exit_call = exit

class obj(AllAbilities):
    logfile = sys.stdout

    def copy(self):
        return self._copy()
    def copies(self,count):
        return self._copies(count)
    def iteritems(self):
        return self._iteritems()
    def keys(self):
        return self._keys()
    def values(self):
        return self._values()
    def inverse(self):
        return self._inverse()
    def set(self,**variables):
        return self._set(**variables)
    def clear(self):
        self._clear()
    def add_attribute(self,name,value=None):
        self._add_attribute(name,value)
    def add_attributes(self,*names,**attributes):
        self._add_attributes(*names,**attributes)
    #def transfer_from(self,other,copy=False):
    #    self._transfer_from(other,copy)
    #def transfer_to(self,other,copy=False):
    #    self._transfer_to(other,copy)
    def append(self,value):
        self._append(value)
    def save(self,fpath=None):
        self._save(fpath)
    def load(self,fpath):
        self._load(fpath)


    def list(self,*names):
        if len(names)==0:
            names = list(self.keys())
            names.sort()
        #end if
        values = []
        for name in names:
            values.append(self[name])
        #end if
        return values
    #end def list

    def tuple(self,*names):
        return tuple(self.list(*names))
    #end def tuple

    def obj(self,*names):
        o = obj()
        o.transfer_from(self,keys=names,copy=False)
        return o
    #end def obj

    def first(self):
        return self[min(self.keys())]
    #end def first

    def last(self):
        return self[max(self.keys())]
    #end def last
    
    def open_log(self,filepath):
        self.logfile = open(filepath,'w')
    #end def open_log

    def close_log(self):
        self.logfile.close()
    #end def close_log

    def _write(self,s):
        self.logfile.write(s)
    #end def _write

    def write(self,s):
        self._write(s)
    #end def write

    def logna(self,*items):
        s=''
        for item in items:
            s+=str(item)+' '
        #end for
        self.logfile.write(s)
    #end def logna

    def log(self,*items):
        s=''
        for item in items:
            s+=str(item)+' '
        #end for
        s+='\n'
        self.logfile.write(s)
    #end def log

    def error(self,message,header=None,exit=True,trace=True,post_header=' Error:'):
        pad = 4*' '
        if header==None:
            header = self.__class__.__name__
        #end if
        self.log(header+post_header)
        self.log(pad+message.replace('\n','\n'+pad))
        if exit:
            self.log('  exiting.\n')
            if trace:
                traceback.print_stack()
            #end if
            exit_call()
        #end if
    #end def error

    def warn(self,message,header=None,post_header=' Warning:'):
        pad = 4*' '
        if header==None:
            header=self.__class__.__name__
        #end if
        self.log(header+post_header)
        self.log(pad+message.replace('\n','\n'+pad))
    #end def error

    @classmethod
    def class_error(cls,message,header=None,exit=True,trace=True,post_header=' Error:'):
        pad = 4*' '
        if header==None:
            header = cls.__name__
        #end if
        cls.logfile.write(header+post_header+'\n')
        cls.logfile.write(('\n'+message).replace('\n','\n'+pad)+'\n')
        if exit:
            cls.logfile.write('  exiting.\n\n')
            if trace:
                traceback.print_stack()
            #end if
            exit_call()
        #end if
    #end def class_error

    @classmethod
    def class_warn(cls,message,header=None,post_header=' Warning:'):
        pad = 4*' '
        if header==None:
            header=cls.__name__
        #end if
        cls.logfile.write(header+post_header+'\n')
        cls.logfile.write(('\n'+message).replace('\n','\n'+pad)+'\n')
    #end def error

    def transfer_from(self,other,keys=None,copy=False):
        if keys==None:
            keys = other.keys()
        #end if
        if not copy:
            for k in keys:
                self[k]=other[k]
            #end for
        else:
            for k in keys:
                self[k]=deepcopy(other[k])
            #end for
        #end if
    #end def transfer_from

    def transfer_to(self,other,keys=None,copy=False):
        if keys==None:
            keys = self.keys()
        #end if
        if not copy:
            for k in keys:
                other[k]=self[k]
            #end for
        else:
            for k in keys:
                other[k]=deepcopy(self[k])
            #end for
        #end if
    #end def transfer_to

    def move_from(self,other,keys=None):
        if keys==None:
            keys = other.keys()
        #end if
        for k in keys:
            self[k]=other[k]
            del other[k]
        #end for
    #end def move_from

    def move_to(self,other,keys=None):
        other.move_from(self,keys)
    #end def move_to

    def copy_from(self,other,keys=None,deep=False):
        self.transfer_from(other,keys,copy=deep)
    #end def copy_from

    def copy_to(self,other,keys=None,deep=False):
        self.transfer_to(other,keys,copy=deep)
    #end def copy_to

    def delete(self,*names):
        if len(names)==1:
            name = names[0]
            value = self[name]
            del self[name]
            return value
        else:
            if len(names)==0:
                names = sorted(obj.keys(self))
            #end if
            values = []
            for name in names:
                values.append(self[name])
                del self[name]
            #end for
            return values
        #end if
    #end def delete

#end class obj




import copy
import cPickle

class generic(genbase):
    logfile = sys.stdout

    def __init__(self,*vals,**kwargs):
        if len(vals)==1 and isinstance(vals[0],(dict,generic)):
            self.add_attributes(**vals[0])
            self.add_attributes(**kwargs)
        else:
            self.add_attributes(*vals,**kwargs)
        #end for
    #end def __init__

    def _dict(self):
        return self.__dict__
    #end def __get_dict

    def _alt(self):
        return self.__dict__
    #end def __alt

    def __len__(self):
        return len(self._dict())
    #end def __len__

    def __contains__(self,name):
        return name in self._dict()
    #end def __contains__

    def __getitem__(self,name):
        return self._dict()[name]
    #end def __getitem__

    def __setitem__(self,name,value):
        self._dict()[name] = value
    #end def __setitem__

    def __delitem__(self,name):
        del self._dict()[name]
    #end def __delitem__

    def __repr__(self):                                                      
        s='' 
        stype = type(s)
        d = self._dict()
        mem = list(d.keys())
        mem.sort()
        for m in mem:
            v=d[m]
            if hasattr(v,'__class__'):
                s+='  {0:<20}  {1:<20}\n'.format(m,v.__class__.__name__)
            else:
                s+='  {0:<20}  {1:<20}\n'.format(m,type(v))
            #end if
        #end for
        return s
    #end def __repr__            

    def __iter__(self):
        d = self._dict()
        for item in d.__dict__:
            yield d[item]
        #end for
    #end def __iter__

    def __str__(self,nindent=1):
        pad = '  '
        npad = nindent*pad
        s=''
        stype = type(s)
        normal = []
        qable  = []
        for k,v in self._dict().iteritems():
            if type(k)!=stype or k[0]!='_':
                if isinstance(v,(generic,obj)):
                    qable.append(k)
                else:
                    normal.append(k)
                #end if
            #end if
        #end for
        normal.sort()
        qable.sort()
        for k in normal:
            v = self[k]
            indent = npad+18*' '
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



    def copy(self):
        return copy.deepcopy(self)
    #end def copy

    def iteritems(self):
        return self._dict().iteritems()
    #end def iteritems

    def keys(self):
        return self._dict().keys()
    #end def keys

    def values(self):
        return self._dict().values()
    #end def keys

    def inverse(self):
        new = self.__class__()
        d = dict((v,k) for k, v in self.iteritems())
        new.add_attributes(**d)
        return new
    #end def inverse

    def set(self,**variables):
        for name,value in variables.iteritems():
            self[name]=value
        #end for
        return self
    #end def set

    def clear(self):
        self._dict().clear()
    #end def clear

    def add_attribute(self,name,value=None):
        self[name] = value
    #end def add_attribute

    def add_attributes(self,*names,**attributes):
        for name in names:
            self[name] = None
        #end for
        for name,value in attributes.iteritems():
            self[name]=value
        #end for
    #end def add_attributes

    def append(self,value):
        self[len(self)] = value
    #end def append

    def save(self,fpath=None):
        if fpath==None:
            fpath='./'+self.__class__.__name__+'.p'
        #end if
        fobj = open(fpath,'w')
        binary = cPickle.HIGHEST_PROTOCOL
        cPickle.dump(self,fobj,binary)
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
        tmp = cPickle.load(fobj)
        fobj.close()
        d = self.__dict__
        d.clear()
        for k,v in tmp.__dict__.iteritems():
            d[k] = v
        #end for
        del fobj
        del tmp
        return
    #end def load

    def transfer_from(self,other,keys=None,copy=False):
        if keys==None:
            keys = other.keys()
        #end if
        if not copy:
            for k in keys:
                self[k]=other[k]
            #end for
        else:
            for k in keys:
                self[k]=deepcopy(other[k])
            #end for
        #end if
    #end def transfer_from

    def transfer_to(self,other,keys=None,copy=False):
        if keys==None:
            keys = self.keys()
        #end if
        if not copy:
            for k in keys:
                other[k]=self[k]
            #end for
        else:
            for k in keys:
                other[k]=deepcopy(self[k])
            #end for
        #end if
    #end def transfer_to

    def move_from(self,other,keys=None):
        if keys==None:
            keys = other.keys()
        #end if
        for k in keys:
            self[k]=other[k]
            del other[k]
        #end for
    #end def move_from

    def move_to(self,other,keys=None):
        other.move_from(self,keys)
    #end def move_to

    def copy_from(self,other,keys=None,deep=False):
        self.transfer_from(other,keys,copy=deep)
    #end def copy_from

    def copy_to(self,other,keys=None,deep=False):
        self.transfer_to(other,keys,copy=deep)
    #end def copy_to

    def delete(self,*names):
        if len(names)==1:
            name = names[0]
            value = self[name]
            del self[name]
            return value
        else:
            if len(names)==0:
                names = sorted(generic.keys(self))
            #end if
            values = []
            for name in names:
                values.append(self[name])
                del self[name]
            #end for
            return values
        #end if
    #end def delete

    def list(self,*names):
        if len(names)==0:
            names = list(generic.keys(self))
            names.sort()
        #end if
        values = []
        for name in names:
            values.append(self[name])
        #end if
        return values
    #end def list

    def tuple(self,*names):
        return tuple(self.list(*names))
    #end def tuple

    def obj(self,*names):
        o = obj()
        o.transfer_from(self,keys=names,copy=False)
        return o
    #end def obj

    def first(self):
        return self[min(self.keys())]
    #end def first

    def last(self):
        return self[max(self.keys())]
    #end def last
    
    def open_log(self,filepath):
        self._alt().logfile = open(filepath,'w')
    #end def open_log

    def close_log(self):
        self._alt().logfile.close()
    #end def close_log

    def write(self,s):
        self._alt().logfile.write(s)
    #end def write

    def logna(self,*items):
        s=''
        for item in items:
            s+=str(item)+' '
        #end for
        self._alt().logfile.write(s)
    #end def logna

    def log(self,*items):
        s=''
        for item in items:
            s+=str(item)+' '
        #end for
        s+='\n'
        self._alt().logfile.write(s)
    #end def log

    def error(self,message,header=None,exit=True,trace=True,post_header=' Error:'):
        pad = 4*' '
        if header==None:
            header = self.__class__.__name__
        #end if
        self.log(header+post_header)
        self.log(pad+message.replace('\n','\n'+pad))
        if exit:
            self.log('  exiting.\n')
            if trace:
                traceback.print_stack()
            #end if
            exit_call()
        #end if
    #end def error

    def warn(self,message,header=None,post_header=' Warning:'):
        pad = 4*' '
        if header==None:
            header=self.__class__.__name__
        #end if
        self.log(header+post_header)
        self.log(pad+message.replace('\n','\n'+pad))
    #end def error

    @classmethod
    def class_error(cls,message,header=None,exit=True,trace=True,post_header=' Error:'):
        pad = 4*' '
        if header==None:
            header = cls.__name__
        #end if
        cls.logfile.write(header+post_header)
        cls.logfile.write(pad+message.replace('\n','\n'+pad)+'\n')
        if exit:
            cls.logfile.write('  exiting.\n\n')
            if trace:
                traceback.print_stack()
            #end if
            exit_call()
        #end if
    #end def class_error



    def _copy(self,*args,**kwargs):              
        return generic.copy(self,*args,**kwargs)
    def _iteritems(self,*args,**kwargs):         
        return generic.iteritems(self,*args,**kwargs)         
    def _keys(self,*args,**kwargs):
        return generic.keys(self,*args,**kwargs)
    def _values(self,*args,**kwargs):
        generic.values(self,*args,**kwargs)
    def _inverse(self,*args,**kwargs):
        return generic.inverse(self,*args,**kwargs)
    def _set(self,*args,**kwargs):
        generic.set(self,*args,**kwargs)
    def _clear(self,*args,**kwargs):
        generic.clear(self,*args,**kwargs)
    def _add_attribute(self,*args,**kwargs):
        generic.add_attribute(self,*args,**kwargs)
    def _add_attributes(self,*args,**kwargs):
        generic.add_attributes(self,*args,**kwargs)
    def _append(self,*args,**kwargs):
        generic.append(self,*args,**kwargs)
    def _save(self,*args,**kwargs):
        generic.save(self,*args,**kwargs)
    def _load(self,*args,**kwargs):
        generic.load(self,*args,**kwargs)
    def _transfer_from(self,*args,**kwargs):
        generic.transfer_from(self,*args,**kwargs)
    def _transfer_to(self,*args,**kwargs):
        generic.transfer_to(self,*args,**kwargs)
    def _move_from(self,*args,**kwargs):
        generic.move_from(self,*args,**kwargs)
    def _move_to(self,*args,**kwargs):
        generic.move_to(self,*args,**kwargs)
    def _copy_from(self,*args,**kwargs):
        generic.copy_from(self,*args,**kwargs)
    def _copy_to(self,*args,**kwargs):
        generic.copy_to(self,*args,**kwargs)
    def _delete(self,*args,**kwargs):
        generic.delete(self,*args,**kwargs)
    def _list(self,*args,**kwargs):
        return generic.list(self,*args,**kwargs)
    def _tuple(self,*args,**kwargs):
        return generic.tuple(self,*args,**kwargs)
    def _obj(self,*args,**kwargs):
        return generic.obj(self,*args,**kwargs)
    def _first(self,*args,**kwargs):
        return generic.first(self,*args,**kwargs)
    def _last(self,*args,**kwargs):
        return generic.last(self,*args,**kwargs)
    def _open_log(self,*args,**kwargs):
        generic.open_log(self,*args,**kwargs)
    def _close_log(self,*args,**kwargs):
        generic.close_log(self,*args,**kwargs)
    def _write(self,*args,**kwargs):
        generic.write(self,*args,**kwargs)
    def _logna(self,*args,**kwargs):
        generic.logna(self,*args,**kwargs)
    def _log(self,*args,**kwargs):
        generic.log(self,*args,**kwargs)
    def _error(self,*args,**kwargs):
        generic.error(self,*args,**kwargs)
    def _warn(self,*args,**kwargs):
        generic.warn(self,*args,**kwargs)
#end class generic



class hidden(generic):
    def __init__(self,*vals,**kwargs):
        d = object.__getattribute__(self,'__dict__')
        d['_hidden_'] = generic()
        d['_public_'] = generic()
        d = self._dict()
        generic.__init__(self,*vals,**kwargs)
    #end def __init__

    def _dict(self):
        return self.__dict__['_public_']
    #end def __get_dict

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
        self._dict()[name] = value
    #end def __setattr__

    def __delattr__(self,name):
        del self._dict()[name]
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
#end class hidden
