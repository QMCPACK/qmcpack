
import numpy as np


# determine if two values differ
def value_diff(v1,v2,tol=1e-6,int_as_float=False):
    diff = False
    v1_bool  = isinstance(v1,(bool,np.bool_))
    v2_bool  = isinstance(v2,(bool,np.bool_))
    v1_int   = isinstance(v1,(int,np.int_)) and not v1_bool
    v2_int   = isinstance(v2,(int,np.int_)) and not v2_bool
    v1_float = isinstance(v1,(float,np.float_))
    v2_float = isinstance(v2,(float,np.float_))
    v1_str   = isinstance(v1,(str,np.string_))
    v2_str   = isinstance(v2,(str,np.string_))
    if id(v1)==id(v2):
        None
    elif int_as_float and (v1_int or v1_float) and (v2_int or v2_float):
        diff = np.abs(float(v1)-float(v2))>tol
    elif v1_float and v2_float:
        diff = np.abs(v1-v2)>tol
    elif v1_int and v2_int:
        diff = v1!=v2
    elif (v1_bool or v1_str) and (v2_bool or v2_str):
        diff = v1!=v2
    elif not isinstance(v1,type(v2)) or not isinstance(v2,type(v1)):
        diff = True
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


# determine if two text blocks differ
def read_text_value(s):
    v = s
    try:
        vi = int(s)
    except:
        try:
            v = float(s)
        except:
            None
        #end try
    #end try
    return v
#end def read_text_value

def read_text_tokens(t):
    tokens = []
    for v in t.split():
        tokens.append(read_text_value(v))
    #end for
    return tokens
#end def read_text_tokens

def text_diff(t1,t2,tol=1e-6,full=False,by_line=False):
    t1 = t1.replace(',',' , ')
    t2 = t2.replace(',',' , ')
    tokens1 = read_text_tokens(t1)
    tokens2 = read_text_tokens(t2)
    diff = value_diff(tokens1,tokens2)
    if not full:
        return diff
    elif not by_line:
        diff1 = dict()
        diff2 = dict()
        nmin = min(len(tokens1),len(tokens2))
        for n,(v1,v2) in enumerate(zip(tokens1[:nmin],tokens2[:nmin])):
            if value_diff(v1,v2):
                diff1[n] = v1
                diff2[n] = v2
            #end if
        #end for
        if len(tokens1)>len(tokens2):
            for n,v in enumerate(tokens1[nmin:]):
                diff1[nmin+n] = v
                diff2[nmin+n] = None
            #end for
        elif len(tokens2)>len(tokens1):
            for n,v in enumerate(tokens2[nmin:]):
                diff1[nmin+n] = None
                diff2[nmin+n] = v
            #end for
        #end if
        return diff,diff1,diff2
    else:
        lines1 = t1.splitlines()
        lines2 = t2.splitlines()
        nmin = min(len(lines1),len(lines2))
        for n,(l1,l2) in enumerate(zip(lines1[:nmin],lines2[:nmin])):
            tokens1 = read_text_tokens(l1)
            tokens2 = read_text_tokens(l2)
            if value_diff(tokens1,tokens2):
                diff1[n] = l1
                diff2[n] = l2
            #end if
        #end for
        if len(lines1)>len(lines2):
            for n,l in enumerate(lines1[nmin:]):
                diff1[nmin+n] = l
                diff2[nmin+n] = None
            #end for
        elif len(lines2)>len(lines1):
            for n,l in enumerate(lines2[nmin:]):
                diff1[nmin+n] = None
                diff2[nmin+n] = l
            #end for
        #end if
        return diff,diff1,diff2
    #end if
#end def text_diff


# print the difference between two objects
def print_diff(o1,o2,text=False,by_line=False): # used in debugging, not actual tests
    from generic import obj
    print(20*'=')
    print(o1)
    print(20*'=')
    print(o2)
    if not text:
        diff,diff1,diff2 = object_diff(o1,o2,full=True)
    else:
        diff,diff1,diff2 = text_diff(o1,o2,full=True,by_line=by_line)
    #end if
    d1 = obj(diff1)
    d2 = obj(diff2)
    print(20*'=')
    print(d1)
    print(20*'=')
    print(d2)
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

text_neq = text_diff
def text_eq(*args,**kwargs):
    return not text_neq(*args,**kwargs)
#end def text_eq


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
            if not file.startswith('.'):
                filepath = os.path.join(test_files_dir,file)
                assert(os.path.exists(filepath))
                storage[file] = filepath
            #end if
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
def setup_unit_test_output_directory(test,subtest,divert=False,file_sets=None):
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

    # divert nexus paths and output, if requested
    if divert:
        from nexus_base import nexus_core
        divert_nexus()
        nexus_core.local_directory  = path
        nexus_core.remote_directory = path
        nexus_core.file_locations = nexus_core.file_locations + [path]
    #end if

    # transfer files into output directory, if requested
    if file_sets is not None:
        assert(isinstance(file_sets,dict))
        filepaths = dict()
        collect_unit_test_file_paths(test,filepaths)
        for fpath,filenames in file_sets.items():
            assert(len(set(filenames)-set(filepaths.keys()))==0)
            dest_path = path
            if fpath is not None:
                dest_path = os.path.join(dest_path,fpath)
                if not os.path.exists(dest_path):
                    os.makedirs(dest_path)
                #end if
            #end if
            assert(os.path.exists(dest_path))
            for filename in filenames:
                source_filepath = filepaths[filename]
                shutil.copy2(source_filepath,dest_path)
                assert(os.path.exists(dest_path))
            #end for
        #end for
    #end if
    return path
#end def setup_unit_test_output_directory



# class used to divert log output when desired
class FakeLog:
    def __init__(self):
        self.reset()
    #end def __init__

    def reset(self):
        self.s = ''
    #end def reset

    def write(self,s):
        self.s+=s
    #end def write

    def close(self):
        None
    #end def close

    def contents(self):
        return self.s
    #end def contents
#end class FakeLog


# dict to temporarily store logger when log output is diverted
logging_storage = dict()

# dict to temporarily store nexus core attributes when diverted
nexus_core_storage    = dict()
nexus_noncore_storage = dict()


# divert nexus log output
def divert_nexus_log():
    from generic import generic_settings,object_interface
    assert(len(logging_storage)==0)
    logging_storage['devlog'] = generic_settings.devlog
    logging_storage['objlog'] = object_interface._logfile 
    logfile = FakeLog()
    generic_settings.devlog   = logfile
    object_interface._logfile = logfile
    return logfile
#end def divert_nexus_log


# restore nexus log output
def restore_nexus_log():
    from generic import generic_settings,object_interface
    assert(set(logging_storage.keys())==set(['devlog','objlog']))
    generic_settings.devlog   = logging_storage.pop('devlog')
    object_interface._logfile = logging_storage.pop('objlog')
    assert(len(logging_storage)==0)
#end def restore_nexus_log


core_keys = [
    'local_directory',
    'remote_directory',
    'mode',
    'stages',
    'stages_set',
    'status',
    'sleep',
    'file_locations',
    'pseudo_dir',
    'pseudopotentials',
    ]
noncore_keys = [
    'pseudo_dir',
    'pseudopotentials',
    ]

# divert nexus core attributes
def divert_nexus_core():
    from nexus_base import nexus_core,nexus_noncore
    assert(len(nexus_core_storage)==0)
    for key in core_keys:
        nexus_core_storage[key] = nexus_core[key]
    #end for
    assert(len(nexus_noncore_storage)==0)
    for key in noncore_keys:
        if key in nexus_noncore:
            nexus_noncore_storage[key] = nexus_noncore[key]
        #end if
    #end for
#end def divert_nexus_core


# restore nexus core attributes
def restore_nexus_core():
    from nexus_base import nexus_core,nexus_noncore
    for key in core_keys:
        nexus_core[key] = nexus_core_storage.pop(key)
    #end for
    assert(len(nexus_core_storage)==0)
    for key in noncore_keys:
        if key in nexus_noncore_storage:
            nexus_noncore[key] = nexus_noncore_storage.pop(key)
        elif key in nexus_noncore:
            del nexus_noncore[key]
        #end if
    #end for
    assert(len(nexus_noncore_storage)==0)
#end def restore_nexus_core


def divert_nexus():
    divert_nexus_log()
    divert_nexus_core()
#end def divert_nexus


def restore_nexus():
    restore_nexus_log()
    restore_nexus_core()
#end def restore_nexus



# declare test failure
#   useful inside try/except blocks
def failed(msg='Test failed.'):
    assert False,msg
#end def failed


class FailedTest(Exception):
    None
#end class FailedTest


global_data = dict(
    job_ref_table = False,
    )


def divert_nexus_errors():
    from generic import generic_settings
    generic_settings.raise_error = True
#end def divert_nexus_errors
