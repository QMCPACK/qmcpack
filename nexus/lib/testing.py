
import numpy as np


# determine if two values differ
def value_diff(v1,v2,tol=1e-6,int_as_float=False):
    diff = False
    if int_as_float and isinstance(v1,(int,float)) and isinstance(v2,(int,float)):
        diff = np.abs(float(v1)-float(v2))>tol
    elif not isinstance(v1,type(v2)) or not isinstance(v2,type(v1)):
        diff = True
    elif isinstance(v1,(bool,int,str)):
        diff = v1!=v2
    elif isinstance(v1,float):
        diff = np.abs(v1-v2)>tol
    elif isinstance(v1,(list,tuple)):
        v1 = np.array(v1,dtype=object).ravel()
        v2 = np.array(v2,dtype=object).ravel()
        if len(v1)!=len(v2):
            diff = True
        else:
            for vv1,vv2 in zip(v1,v2):
                diff |= value_diff(vv1,vv2,tol,int_as_float)
            #end for
        #end if
    elif isinstance(v1,np.ndarray):
        v1 = v1.ravel()
        v2 = v2.ravel()
        if len(v1)!=len(v2):
            diff = True
        else:
            for vv1,vv2 in zip(v1,v2):
                diff |= value_diff(vv1,vv2,tol,int_as_float)
            #end for
        #end if
    elif isinstance(v1,dict):
        k1 = v1.keys()
        k2 = v2.keys()
        if set(k1)!=set(k2):
            diff = True
        else:
            for k in k1:
                diff |= value_diff(v1[k],v2[k],tol,int_as_float)
            #end for
        #end if
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
def object_diff(o1,o2,tol=1e-6,full=False,int_as_float=False,bypass=False):
    diff1 = dict()
    diff2 = dict()
    if not bypass:
        o1 = o1._serial().__dict__
        o2 = o2._serial().__dict__
    #end if
    keys1 = set(o1.keys())
    keys2 = set(o2.keys())
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
    from generic import obj
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



# find the path to the Nexus directory and other internal paths
def nexus_path(append=None,location=None):
    import os
    testing_path = os.path.realpath(__file__)

    assert(isinstance(testing_path,str))
    assert(len(testing_path)>0)
    assert('/' in testing_path)

    tokens = testing_path.split('/')

    assert(len(tokens)>=3)
    assert(tokens[-1].startswith('testing.py'))
    assert(tokens[-2]=='lib')
    assert(tokens[-3]=='nexus')

    path = os.path.dirname(testing_path)
    path = os.path.dirname(path)

    assert(path.endswith('/nexus'))

    if location is not None:
        if location=='unit':
            append = 'tests/unit'
        else:
            print('nexus location "{}" is unknown'.format(location))
            raise ValueError
        #end if
    #end if
    if append is not None:
        path = os.path.join(path,append)
    #end if

    assert(os.path.exists(path))

    return path
#end def nexus_path



# find the path to a file associated with a unit test
def unit_test_file_path(test,file=None):
    import os
    unit_path  = nexus_path(location='unit')
    files_dir  = 'test_{}_files'.format(test)
    path = os.path.join(unit_path,files_dir)
    if file is not None:
        path = os.path.join(path,file)
    #end if
    assert(os.path.exists(path))
    return path
#end def unit_test_file_path



# collect paths to all files associated with a unit test
def collect_unit_test_file_paths(test,storage):
    import os
    if len(storage)==0:
        test_files_dir = unit_test_file_path(test)
        files = os.listdir(test_files_dir)
        for file in files:
            filepath = os.path.join(test_files_dir,file)
            assert(os.path.exists(filepath))
            storage[file] = filepath
        #end for
    #end if
    return storage
#end def collect_unit_test_file_paths



# find the output path for a test
def unit_test_output_path(test,subtest=None):
    import os
    unit_path  = nexus_path(location='unit')
    files_dir  = 'test_{}_output'.format(test)
    path = os.path.join(unit_path,files_dir)
    if subtest is not None:
        path = os.path.join(path,subtest)
    #end if
    return path
#end def unit_test_output_path



# setup the output directory for a test
def setup_unit_test_output_directory(test,subtest):
    import os
    import shutil
    path = unit_test_output_path(test,subtest)
    assert('nexus' in path)
    assert('unit' in path)
    assert(os.path.basename(path).startswith('test_'))
    assert(path.endswith('/'+subtest))
    if os.path.exists(path):
        shutil.rmtree(path)
    #end if
    os.makedirs(path)
    assert(os.path.exists(path))
    return path
#end def setup_unit_test_output_directory
