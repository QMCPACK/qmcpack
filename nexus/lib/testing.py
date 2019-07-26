
import numpy as np


# determine if two values differ
def value_diff(v1,v2,tol=1e-6,int_as_float=False):
    diff = False
    if int_as_float and isinstance(v1,(int,float)) and isinstance(v2,(int,float)):
        diff = np.abs(float(v1)-float(v2))>tol
    elif not isinstance(v1,type(v2)):
        diff = True
    elif isinstance(v1,(bool,int,str)):
        diff = v1!=v2
    elif isinstance(v1,float):
        diff = np.abs(v1-v2)>tol
    elif isinstance(v1,(list,tuple)):
        v1 = np.array(v1,dtype=object).ravel()
        v2 = np.array(v2,dtype=object).ravel()
        for vv1,vv2 in zip(v1,v2):
            diff |= value_diff(vv1,vv2,tol,int_as_float)
        #end for
    elif isinstance(v1,np.ndarray):
        v1 = v1.ravel()
        v2 = v2.ravel()
        for vv1,vv2 in zip(v1,v2):
            diff |= value_diff(vv1,vv2,tol,int_as_float)
        #end for
    elif isinstance(v1,dict):
        diff = v1!=v2
    elif isinstance(v1,set):
        diff = v1!=v2
    elif v1 is None and v2 is None:
        diff = False
    elif hasattr(v1,'__len__') and hasattr(v2,'__len__') and len(v1)==0 and len(v2)==0:
        None
    else:
        diff = True # unsupported types
    #end if
    return diff
#end def value_diff


# determine if two objects differ
def object_diff(o1,o2,tol=1e-6,full=False,int_as_float=False):
    diff1 = dict()
    diff2 = dict()
    o1    = o1._serial()
    o2    = o2._serial()
    keys1 = set(o1._keys())
    keys2 = set(o2._keys())
    ku1   = keys1 - keys2
    ku2   = keys2 - keys1
    km    = keys1 & keys2
    for k in ku1:
        diff1[k] = o1[k]
    #end for
    for k in ku2:
        diff2[k] = o2[k]
    #end for
    for k in km:
        v1 = o1[k]
        v2 = o2[k]
        if value_diff(v1,v2,tol,int_as_float):
            diff1[k] = v1
            diff2[k] = v2
        #end if
    #end for
    diff = len(diff1)!=0 or len(diff2)!=0
    if not full:
        return diff
    else:
        return diff,diff1,diff2
    #end if
#end def object_diff


# print the difference between two objects
def print_diff(o1,o2): # used in debugging, not actual tests
    from nexus import obj
    print 20*'='
    print o1
    print 20*'='
    print o2
    diff,diff1,diff2 = object_diff(o1,o2,full=True)
    d1 = obj(diff1)
    d2 = obj(diff2)
    print 20*'='
    print d1
    print 20*'='
    print d2
#end def print_diff



# additional convenience functions to use value_diff and object_diff
value_neq = value_diff
def value_eq(*args,**kwargs):
    return not value_neq(*args,**kwargs)
#end def value_eq

object_neq = object_diff
def object_eq(*args,**kwargs):
    return not object_neq(*args,**kwargs)
#end def object_eq




