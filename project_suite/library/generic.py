
import sys
import traceback
from copy import deepcopy
from abilities import AllAbilities

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
#end class obj
