
from numpy import array,ndarray,minimum,abs,ix_,resize
import sys
import keyword
from inspect import getmembers

from superstring import valid_variable_name
from generic import obj
from developer import DevBase,unavailable
try:
    import h5py
except ImportError:
    h5py = unavailable('h5py')
#end try
from debug import *



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
            for k,v in self._datasets.iteritems():
                s+= '    '+k+'\n'
            #end for
        #end if
        if len(self._groups)>0:
            s+= '  groups:\n'
            for k,v in self._groups.iteritems():
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
            for name,value in self.iteritems():
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


    #project interface methods

    def zero(self,*names):
        for name in names:
            if name in self and isinstance(self[name],ndarray):
                self[name][:] = 0
            #end if
        #end for
        for name,value in self._groups.iteritems():
            if isinstance(value,str):
                ci(ls(),gs())
            #end if
            value.zero(*names)
        #end for
        #self.sum(*names)
    #end def zero


    def minsize(self,other,*names):
        name_set = set(names)
        snames = set(self.keys()) & name_set
        onames = set(other.keys()) & name_set
        if len(snames)>0 and len(snames)<len(names):
            self.error('only some names found in self\n  names_requested:'+str(names)+'\n  names found:  '+str(snames))
        elif len(onames)>0 and len(onames)<len(names):
            self.error('only some names found in other\n  names_requested:'+str(names)+'\n  names found:  '+str(onames))
        #end if
        if len(snames)>0:
            for name in names:
                svalue = self[name]
                ovalue = other[name]
                if not isinstance(svalue,ndarray) or not isinstance(ovalue,ndarray):
                    self.error(name+' is not an array')
                #end if
                shape  = minimum(svalue.shape,ovalue.shape)
                self[name] = resize(svalue,shape)
            #end for
        #end if
        for name,value in self._groups.iteritems():
            if name in other and isinstance(other[name],HDFgroup):
                value.minsize(other[name])
            else:
                self.error(name+' not found in minsize partner')
            #end if
        #end for
        #self.sum(*names)
    #end def minsize


    def accumulate(self,other,*names):
        name_set = set(names)
        snames = set(self.keys()) & name_set
        onames = set(other.keys()) & name_set
        if len(snames)>0 and len(snames)<len(names):
            self.error('only some names found in self\n  names_requested:'+str(names)+'\n  names found:  '+str(snames))
        elif len(onames)>0 and len(onames)<len(names):
            self.error('only some names found in other\n  names_requested:'+str(names)+'\n  names found:  '+str(onames))
        #end if
        if len(snames)>0:
            for name in names:
                svalue = self[name]
                ovalue = other[name]
                if not isinstance(svalue,ndarray) or not isinstance(ovalue,ndarray):
                    self.error(name+' is not an array')
                #end if
                shape  = minimum(svalue.shape,ovalue.shape)
                if abs(shape-array(svalue.shape)).sum() > 0:
                    self.error(name+' in partner is too large')
                #end if
                ranges = []
                for s in shape:
                    ranges.append(range(s))
                #end for
                #add the part of the other data that fits into own data
                svalue += ovalue[ix_(*ranges)]
            #end for
        #end if
        for name,value in self._groups.iteritems():
            if name in other and isinstance(other[name],HDFgroup):
                value.accumulate(other[name])
            else:
                self.error(name+' not found in accumulate partner')
            #end if
        #end for
        #self.sum(*names)
    #end def accumulate

    
    def normalize(self,normalization,*names):
        for name in names:
            if name in self and isinstance(self[name],ndarray):
                self[name] /= normalization
            #end if
        #end for
        for name,value in self._groups.iteritems():
            value.normalize(normalization,*names)
        #end for
        #self.sum(*names)
    #end def normalize

        
    def sum(self,*names):
        for name in names:
            if name in self and isinstance(self[name],ndarray) and name=='value':
                s = self[name].mean(0).sum()
                print '                sum = {0}'.format(s)
            #end if
        #end for
    #end def sum

#end class HDFgroup




class HDFreader(DevBase):
    datasets = set(["<class 'h5py.highlevel.Dataset'>","<class 'h5py._hl.dataset.Dataset'>"])
    groups   = set(["<class 'h5py.highlevel.Group'>","<class 'h5py._hl.group.Group'>"])
    
    def __init__(self,fpath,verbose=False):
        
        if verbose:
            print '  Initializing HDFreader'

        self.fpath=fpath
        if verbose:
            print '    loading h5 file'

        try:
            self.hdf = h5py.File(fpath,'r')
        except IOError:
            self._success = False
            self.hdf = obj(obj=obj())
        else:
            self._success = True
        #end if

        if verbose:
            print '    converting h5 file to dynamic object'
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
            for kr,v in hcur.iteritems():
                k=cur._escape_name(kr)
                if valid_variable_name(k):
                    vtype = str(type(v))
                    if vtype in HDFreader.datasets:
                        self.add_dataset(cur,k,v)
                    elif vtype in HDFreader.groups:
                        self.add_group(hcur,cur,k,v)
                    else:
                        print 'hdfreader error: encountered invalid type: '+vtype
                        sys.exit()
                    #end if
                else:
                    print 'hdfreader warning: attribute '+k+' is not a valid variable name and has been ignored'                    
                #end if
            #end for
        #end if

        if verbose:
            print '  end HDFreader Initialization'

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
        exec 'cur.'+k+'=array(v)'
        exec 'cur._add_dataset("'+k+'",cur.'+k+')'
        return
    #end def add_dataset

    def add_group(self,hcur,cur,k,v):
        exec 'cur.'+k+'=HDFgroup()'
        exec 'cur._add_group("'+k+'",cur.'+k+')'
        cur._groups[k]._parent = cur
        self.increment_level()
        self.cur[self.ilevel]  = cur._groups[k]
        self.hcur[self.ilevel] = hcur[k]

        cur   = self.cur[self.ilevel]
        hcur  = self.hcur[self.ilevel]
        for kr,v in hcur.iteritems():
            k=cur._escape_name(kr)
            if valid_variable_name(k):
                vtype = str(type(v))
                if vtype in HDFreader.datasets:
                    self.add_dataset(cur,k,v)
                elif vtype in HDFreader.groups:
                    self.add_group(hcur,cur,k,v)
                #end if
            else:
                print 'hdfreader warning: attribute '+k+' is not a valid variable name and has been ignored'                    
            #end if
        #end for

        return
    #end def add_group
#end class HDFreader




