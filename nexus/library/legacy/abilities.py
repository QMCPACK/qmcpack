##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  abilities.py                                                      # 
#    Base classes for all Nexus classes.  Defines abilities such as  #
#    querying, printing, etc.  Lays the foundation for developer     #
#    and user environments.                                          #
#                                                                    #
#  Content summary:                                                  #
#    genbase                                                         #
#      Essentially a forward declaration of the generic class.       #
#      See generic.py.                                               #
#                                                                    #
#    Callable                                                        #
#      Helper class to create 'static' methods.  Seldom used.        #
#                                                                    #
#    Copyable                                                        #
#      Base class for copy abilities.                                #
#                                                                    #
#    Comparable                                                      #
#      Base class for comparison abilities (=)                       #
#                                                                    #
#    Iterable                                                        #
#      Base class for iteration abilities.                           #
#      Allows derived classes to behave somewhat like dictionaries.  #
#                                                                    #
#    Queryable                                                       #
#      Base class for query abilities.                               #
#      Covers len, in, del, clear, etc.                              #
#                                                                    #
#    Logable                                                         #
#      Base class for generic logging abilities.                     #
#      Covers log, error, warn, etc.                                 #
#      Functionality superceded by obj class in generic.py.          #
#                                                                    #
#    Saveable                                                        #
#      Base class for generic save/load abilities.                   #
#                                                                    #
#    Validatable                                                     #
#      Base class for validation abilities.  Seldom used.            #
#                                                                    #
#    AllAbilities                                                    #
#      Summary base class of multiple abilities.                     #
#      Starting base class for most Nexus classes.                   #
#      Confluence of Copyable, Comparable, Queryable, Savable, and   #
#        Validatable classes.                                        #
#                                                                    # 
#====================================================================#



class genbase(object):
    None
#end class genbase


class Callable:                        #helper class to create 'static' methods
    def __init__(self, anycallable):
        self.__call__ = anycallable
    #end def init
#end class Callable


import copy
class Copyable:
    def _copy(self):
        return copy.deepcopy(self)
        #return copy.copy(self)
    #end def copy

    def _copies(self,count):
        copies = []
        for i in range(count):
            copies.append(self._copy())
        #end for
        return copies
    #end def copies
#end class Copyable


from numpy import array
class Comparable:
    array_type = type(array([1]))
    def __eq__(self,other): 
        if type(other)!=type(self):#must both be class objects
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
            if stype == Comparable.array_type:
                eq = eq and (svar==ovar).all()
            else:
                eq = eq and svar==ovar
            #end if
        #end for
        return eq
    #end def __eq__
#end class Comparable


class Iterable:
    def __iter__(self):
        for item in self.__dict__:
            yield self.__dict__[item]
        #end for
    #end def __iter__

    def _iteritems(self):
        return self.__dict__.iteritems()
    #end def iteritems

    def _keys(self):
        return self.__dict__.keys()
    #end def _keys

    def _values(self):
        return self.__dict__.values()
    #end def _keys

    def _inverse(self):
        new = self.__class__()
        new.__dict__ = dict((v,k) for k, v in self._iteritems())
        return new
    #end def _inverse
#end class Iterable


import inspect
class Queryable(Iterable):
    def __str__(self,nindent=1):
        pad = '  '
        npad = nindent*pad
        s=''
        stype = type(s)
        normal = []
        qable  = []
        for k,v in self._iteritems():
            if type(k)!=stype or k[0]!='_':
                if isinstance(v,(Queryable,genbase)):
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

    def __repr__(self):
        s=''
        stype = type(s)
        mem = list(self.__dict__.keys())
        mem.sort()
        for m in mem:
            if type(m)!=stype or m[0]!='_':
#                s+='  '+m+'   '+str(type(self.__dict__[m]))+'\n'
#                s+='  '+m+'   '+self.__dict__[m].__class__.__name__+'\n'
                v=self.__dict__[m]
                if hasattr(v,'__class__'):
                    s+='  {0:<20}  {1:<20}\n'.format(m,v.__class__.__name__)
                else:
                    s+='  {0:<20}  {1:<20}\n'.format(m,type(v))
                #end if
            #end if
        #end for
        return s
    #end def __repr__

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
    #end def __setitem__

    def _set(self,**variables):
        for name,value in variables.iteritems():
            self[name]=value
        #end for
        return self
    #end def _set

    def _clear(self):
        self.__dict__.clear()
    #end def _clear

    def _add_attribute(self,name,value=None):
        self[name] = value
    #end def _add_attribute

    def _add_attributes(self,*names,**attributes):
        for name in names:
            self[name] = None
        #end for
        for name,value in attributes.iteritems():
            self[name]=value
        #end for
    #end def _add_attributes

    def _transfer_from(self,other,copy=False):
        if not copy:
            for k,v in other._iteritems():
                self[k]=v
            #end for
        else:
            for k,v in other._iteritems():
                self[k]=deepcopy(v)
            #end for
        #end if
    #end def _transfer_from

    def _transfer_to(self,other,copy=False):
        other._transfer_from(self,copy)
    #end def _transfer_to

    def _append(self,value):
        self[len(self)] = value
    #end def _append
#end class Queryable


import sys
class Logable:
    def _open_log(self,filepath,increment=0,prefix=''):
        self._log_file=open(filepath,'w')
        self._log_pad_prefix=prefix
        self._log_pad_increment=increment
        self._log_pad_shift=0
        self._log_pad_update()
    #end def _open_log

    def _close_log(self):
        if self._log_file!=sys.stdout:
            self._log_file.close()
        #end if
    #end def _close_log

    def _set_log(self,logfile=sys.stdout,increment=0,prefix=''):
        self._log_file=logfile
        self._log_pad_prefix=prefix
        self._log_pad_increment=increment
        self._log_pad_shift=0
        self._log_pad_update()
    #end def _set_log

    def _get_log(self):
        L = Logable()
        L._log_file          = self._log_file               
        L._log_pad_prefix    = str(self._log_pad_prefix)
        L._log_pad_increment = int(self._log_pad_increment)
        L._log_pad_shift     = int(self._log_pad_shift)
        return L
    #end def _get_log

    def _transfer_log(self,L):
        self._log_file          = L._log_file               
        self._log_pad_prefix    = L._log_pad_prefix   
        self._log_pad_increment = L._log_pad_increment
        self._log_pad_shift     = L._log_pad_shift    
        self._log_pad_update()
    #end def _transfer_log

    def _increment_pad(self):
        self._log_pad_shift+=self._log_pad_increment
        self._log_pad_update()
    #end def _increment_pad

    def _decrement_pad(self):
        self._log_pad_shift-=self._log_pad_increment
        self._log_pad_update()
    #end def _decrement_pad

    def _set_pad(self,increment=0,prefix=''):
        self._log_pad_prefix=prefix
        self._log_pad_increment=increment
        self._log_pad_update()
    #end def _set_pad

    def _log(self,s=''):
        self._log_file.write(self._log_pad+str(s)+'\n')
    #end def _log

    def _exit(self,s=''):
        self._log(self.__class__.__name__+' exiting')
        self._log(s)
        sys.exit()
    #end def _exit

    def _error(self,s='',exit_now=False):
        self._log(self.__class__.__name__+' Error:  '+s)
        if exit_now:
            self._exit()
        #end if
    #end def _error

    def _warn(self,s):
        self._log(self.__class__.__name__+' Warning:  '+s)
        return
    #end def _warn


    #user should not need to use these functions
    def _log_pad_update(self):
        self._log_pad = self._log_pad_prefix+self._log_pad_shift*' '
    #end def _log_pad_update
#end class Logable


import cPickle
import pickle
class Savable:
    def _save(self,fpath=None,fast=True):
        #self._unlink_dynamic_methods(self)
        if fpath==None:
            fpath='./'+self.__class__.__name__+'.p'
        #end if
        #self._savepath=fpath
        fobj = open(fpath,'w')
        if fast:
            binary = cPickle.HIGHEST_PROTOCOL
            cPickle.dump(self,fobj,binary)
        else:
            binary = pickle.HIGHEST_PROTOCOL
            pickle.dump(self,fobj,binary)
        #end if
        fobj.close()
        del fobj
        del binary
        return
    #end def _save

    def _load(self,fpath=None,fast=True):
        if fpath==None:
            fpath='./'+self.__class__.__name__+'.p'
        #end if
        fobj = open(fpath,'r')
        if fast:
            tmp = cPickle.load(fobj)
        else:
            tmp = pickle.load(fobj)
        #end if
        fobj.close()
        self.__dict__.clear()
        #self.__dict__=tmp.__dict__
        for k,v in tmp.__dict__.iteritems():
            self.__dict__[k] = v
        #end for
        del fobj
        del tmp
        #self._relink_dynamic_methods()
        return
    #end def _load

    #def _unlink_dynamic_methods(self,obj):
    #    if hasattr(obj,'__dict__'):
    #        for k,v in obj.__dict__.iteritems():
    #            if type(v)==type(self._unlink_dynamic_methods):
    #                obj.__dict__[k]=None
    #            else:
    #                self._unlink_dynamic_methods(v) 
    #            #end if
    #        #end for
    #    elif hasattr(obj,'__iter__'):
    #        for a in obj:
    #            self._unlink_dynamic_methods(a)
    #        #end for
    #    #end if
    #    return
    ##end def _unlink_dynamic_methods
    #
    #def _relink_dynamic_methods(self,obj):
    #    if hasattr(obj,'_set_dynamic_methods'):
    #        obj._set_dynamic_methods()
    #    #end if
    #    if hasattr(obj,'__dict__'):
    #        for a in obj.__dict__:
    #            self._relink_dynamic_methods(a)
    #        #end for
    #    elif hasattr(obj,'__iter__'):
    #        for a in obj:
    #            self._relink_dynamic_methods(a)
    #        #end for
    #    #end if
    #    return
    ##end def _relink_dynamic_methods
    #
    #def _set_dynamic_methods(self):
    #    return
    ##end def _set_dynamic_methods
#end class Savable


class Validatable:
    def _has_type_bindings(self):
        return '_type_bindings' in self.__class__.__dict__
    #end def _has_type_bindings

    def _has_rules(self):
        return '_rules' in self.__class__.__dict__
    #end def _has_rules

    def _bind_type(self,vname,vtype):
        if not self._has_type_bindings():
            self.__class__.__dict__['_type_bindings']=dict()
        #end if
        self.__class__.__dict__['_type_bindings'][vname]=vtype
        return
    #end def _bind_type

    #rules must have an executable form, eg. vrule='self.var > 1.0'
    def _add_rule(self,vrule):
        if not self._has_rules():
            self.__class__.__dict__['_rules']=list()
        #end if
        self.__class__.__dict__['_rules'].append(vrule)
        return
    #end def

    def _is_complete(self):
        complete = True
        for k,v in self.__dict__.iteritems():
            if k[0]!='_':
                complete = complete and v!=None
            #end if
        #end for
        return complete
    #end def _is_complete

    def _is_type_valid(self,return_exceptions=False):
        type_valid=True
        exceptions = dict()
        if self._has_type_bindings():
            for k,v in self.__class__.__dict__['_type_bindings'].iteritems():
                var = self.__dict__[k]
                if hasattr(v,'__class__'):
                    vtype = var.__class__.__name__
                else:
                    vtype = type(var)
                #end if
                var_valid = vtype==v
                type_valid = type_valid and var_valid
                if not var_valid:
                    exceptions[k]=v,vtype
                #end if
            #end for
        #end if
        if not return_exceptions:
            return type_valid
        else:
            return type_valid,exceptions
        #end if
    #end def _is_type_valid

    def _is_legal(self):
        legal=True
        if self._has_rules():
            for rule in self.__class__.__dict__['_rules']:
                exec 'obeys_rule='+rule
                legal = legal and obeys_rule
            #end for
        #end if
        return legal
    #end def _is_legal

    def _is_valid(self,exit_on_fail=False,return_exceptions=False):
        complete   = self._is_complete()
        type_valid,type_exceptions = self._is_type_valid(return_exceptions=True)
        legal      = self._is_legal()
        valid = complete and type_valid and legal
        if not valid and exit_on_fail:
            name = self.__class__.__name__
            if not complete:
                print '       '+name+' has None attributes'
            #end if
            if not type_valid:
                print '       '+name+' has invalid attribute types'
                self._write_type_exceptions(type_exceptions)
            #end if
            if not legal:
                print '       '+name+' violates user-defined rules'
            #end if
            sys.exit('Error: '+name+' failed _is_valid, exiting')
        #end if
        return valid
    #end def _is_valid

    def _write_type_exceptions(self,type_exceptions,pad='         '):
        for k,v in type_exceptions.iteritems():
            print pad+'variable '+k+' should be '+v[0]+' type, is '+v[1]+' type'
        #end for
    #end def _write_type_exceptions        
#end class Validatable


class AllAbilities(Copyable,Comparable,Queryable,Savable,Validatable):
    def __init__(self,*vals,**attributes):
        names = []
        for val in vals:
            if isinstance(val,str):
                names.append(val)
            elif isinstance(val,Iterable):
                for name,value in val._iteritems():
                    self[name]=value
                #end for
            elif isinstance(val,dict):
                for name,value in val.iteritems():
                    self[name]=value
                #end for
            #end for
        #end for
        self._add_attributes(*names,**attributes)
    #end def __init__
#end class AllAbilities




class inTest(Comparable):
    count=0
    def __init__(self):
        inTest.count+=1
        count = inTest.count
        self.a=array([count,2,3,4,5.0])
        self.l=['list',3,1.3,count]
        self.d={'dict':count,'b':2.0,'c':'cval'}
        self.s=set(['set','apple','banana','pear',count])
        self.t=('tup',1,2.0,count)
        self.i = 4+count
        self.r = 3.14+count
        self.st= 'string'+str(count)
        return
    #end def __init__
#end class inTest

class Test(AllAbilities):
    def __init__(self):
        self.a=array([1,2,3,4,5.0])
        self.l=['list',3,1.3,inTest()]
        self.d={'dict':1,'b':2.0,'c':'cval','o':inTest()}
        self.s=set(['set','apple','banana','pear'])
        self.t=('tup',1,2.0,inTest())
        self.i = 4
        self.r = 3.14
        self.st= 'string'
        self.o = inTest()
        return
    #end def __init__
#end class Test

def test():
    s=True
    it = inTest()
    t = Test()
    t2= Test()

    print
    print 'copying'
    tc = t._copy()

    print
    print 'comparing'
    istrue = t==tc
    isfalse= t==t2
    print '  test == copy ',istrue
    print '  test == test2',isfalse
    s=s and istrue and not isfalse

    print
    print 'iterating'
    for v in t:
        print v
    #end for
    for k,v in t._iteritems():
        print k,'=',v
    #end for

    print
    print 'querying'
    qs=True
    for k,v in t._iteritems():
        qs=qs and k in t
    #end for
    s=s and qs
    print '  in ',qs
    print '  __str__'
    print t
    print '  __repr__'
    print t.__repr__
    print '  adding attribute myatt=5'
    t._add_attribute('myatt',5)
    print t
    print '  adding attribute i=2342'
    t._add_attribute('myatt',2342)
    print t


    print
    print 'saving'
    t._save()
    tt=Test()
    tt._load(t._savepath)
    print '  save==load',t==tt
    print


    return s
#end def test
#del inTest
#del Test

