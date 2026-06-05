try:
    import numpy as np
    if np.lib.NumpyVersion(np.__version__) >= '2.0.0b1':
        np.set_printoptions(legacy="1.25")
    numpy_available = True
except:
    numpy_available = False
#end try


def_atol =  0.0
def_rtol = 1e-6


# determine if two floats differ
def float_diff(v1,v2,atol=def_atol,rtol=def_rtol):
    return np.abs(v1-v2)>atol+rtol*np.abs(v2)
#end def float_diff


# determine if two values differ
def value_diff(v1,v2,atol=def_atol,rtol=def_rtol,int_as_float=False):
    diff = False
    v1_bool  = isinstance(v1,(bool,np.bool_))
    v2_bool  = isinstance(v2,(bool,np.bool_))
    v1_int   = isinstance(v1,(int,np.int_)) and not v1_bool
    v2_int   = isinstance(v2,(int,np.int_)) and not v2_bool
    v1_float = isinstance(v1,(float,np.float64))
    v2_float = isinstance(v2,(float,np.float64))
    v1_str   = isinstance(v1,(str,np.bytes_))
    v2_str   = isinstance(v2,(str,np.bytes_))
    if id(v1)==id(v2):
        None
    elif int_as_float and (v1_int or v1_float) and (v2_int or v2_float):
        diff = float_diff(v1,v2,atol=atol,rtol=rtol)
    elif v1_float and v2_float:
        diff = float_diff(v1,v2,atol=atol,rtol=rtol)
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
                diff |= value_diff(vv1,vv2,atol,rtol,int_as_float)
            #end for
        #end if
    elif isinstance(v1,np.ndarray):
        v1 = v1.ravel()
        v2 = v2.ravel()
        if len(v1)!=len(v2):
            diff = True
        else:
            for vv1,vv2 in zip(v1,v2):
                diff |= value_diff(vv1,vv2,atol,rtol,int_as_float)
            #end for
        #end if
    elif isinstance(v1,dict):
        k1 = v1.keys()
        k2 = v2.keys()
        if set(k1)!=set(k2):
            diff = True
        else:
            for k in k1:
                diff |= value_diff(v1[k],v2[k],atol,rtol,int_as_float)
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
def object_diff(o1,o2,atol=def_atol,rtol=def_rtol,int_as_float=False,full=False,bypass=False):
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
        if value_diff(v1,v2,atol,rtol,int_as_float):
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
            pass
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

def text_diff(t1,t2,atol=def_atol,rtol=def_rtol,int_as_float=False,full=False,by_line=False):
    t1 = t1.replace(',',' , ')
    t2 = t2.replace(',',' , ')
    tokens1 = read_text_tokens(t1)
    tokens2 = read_text_tokens(t2)
    diff = value_diff(tokens1,tokens2,atol,rtol,int_as_float)
    if not full:
        return diff
    elif not by_line:
        diff1 = dict()
        diff2 = dict()
        nmin = min(len(tokens1),len(tokens2))
        for n,(v1,v2) in enumerate(zip(tokens1[:nmin],tokens2[:nmin])):
            if value_diff(v1,v2,atol,rtol,int_as_float):
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
        diff1 = dict()
        diff2 = dict()
        lines1 = t1.splitlines()
        lines2 = t2.splitlines()
        nmin = min(len(lines1),len(lines2))
        for n,(l1,l2) in enumerate(zip(lines1[:nmin],lines2[:nmin])):
            tokens1 = read_text_tokens(l1)
            tokens2 = read_text_tokens(l2)
            if value_diff(tokens1,tokens2,atol,rtol,int_as_float):
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
def print_diff(o1,o2,atol=def_atol,rtol=def_rtol,int_as_float=False,text=False,by_line=False): # used in debugging, not actual tests
    from .generic import obj
    hline = '========== {} =========='
    print(hline.format('left object'))
    print(o1)
    print(hline.format('right object'))
    print(o2)
    if not text:
        diff,diff1,diff2 = object_diff(o1,o2,atol,rtol,int_as_float,full=True)
    else:
        diff,diff1,diff2 = text_diff(o1,o2,atol,rtol,int_as_float,full=True,by_line=by_line)
    #end if
    d1 = obj(diff1)
    d2 = obj(diff2)
    print(hline.format('left diff'))
    print(list(sorted(d1.keys())))
    print(d1)
    print(hline.format('right diff'))
    print(list(sorted(d2.keys())))
    print(d2)
#end def print_diff


# check for value equality and if different, print the difference
def check_value_eq(v1,v2,**kwargs):
    verbose = kwargs.pop('verbose',False)
    same = value_eq(v1,v2,**kwargs)
    if not same and (verbose or global_data['verbose']):
        print('\nValues differ, please see below for details')
        hline = '========== {} =========='
        print()
        print(hline.format('left value'))
        print(v1)
        print()
        print(hline.format('right value'))
        print(v2)
        print()
    #end if
    return same
#end def check_value_eq


# check for object equality and if different, print the difference
def check_object_eq(o1,o2,**kwargs):
    verbose = kwargs.pop('verbose',False)
    same = object_eq(o1,o2,**kwargs)
    if not same and (verbose or global_data['verbose']):
        print('\nObjects differ, please see below for details')
        print_diff(o1,o2)
    #end if
    return same
#end def check_object_eq



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


# declare test failure
#   useful inside try/except blocks
def failed(msg='Test failed.'):
    assert False,msg
#end def failed


class FailedTest(Exception):
    None
#end class FailedTest


global_data = dict(
    verbose       = False,
    job_ref_table = False,
    )


def clear_all_sims():
    from .simulation import Simulation
    Simulation.clear_all_sims()
#end def clear_all_sims



def check_final_state():
    from .nexus_base import nexus_core,nexus_core_defaults
    from .nexus_base import nexus_noncore,nexus_noncore_defaults
    from .nexus_base import nexus_core_noncore,nexus_core_noncore_defaults
    
    assert('runs' in nexus_core_defaults)
    assert('basis_dir' in nexus_noncore_defaults)
    assert('pseudo_dir' in nexus_core_noncore_defaults)

    assert(object_eq(nexus_core,nexus_core_defaults))
    assert(object_eq(nexus_noncore,nexus_noncore_defaults))
    assert(object_eq(nexus_core_noncore,nexus_core_noncore_defaults))

    from .simulation import Simulation

    assert(Simulation.sim_count==0)
    assert(len(Simulation.all_sims)==0)
    assert(len(Simulation.sim_directories)==0)
#end def check_final_state



def execute(command):
    from .execute import execute as nexus_execute
    import os
    # for python exe's, restrict pythonpath to this nexus repo
    ct = command.split()
    if len(ct)>0 and ct[0].strip('3').endswith('python'):
        pypath = os.path.abspath(os.path.join(__file__,'..','..'))
        command = 'PYTHONPATH='+pypath+' '+command
    out,err,rc = nexus_execute(command)
    if rc!=0:
        msg = '''Executed system command failed.

Command:
========
{}

stdout:
=======
{}

stderr:
=======
{}

Return code:
============
{}

'''.format(command,out,err,rc)
        failed(msg)
    #end if
    return out,err,rc
#end def execute
