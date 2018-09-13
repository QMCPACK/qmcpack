
import os
import sys
import traceback
from numpy import ndarray,array,abs


def value_diff(v1,v2,tol=1e-6):
    diff = False
    if not isinstance(v1,type(v2)):
        diff = True
    elif isinstance(v1,(bool,int,str)):
        diff = v1!=v2
    elif isinstance(v1,float):
        diff = abs(v1-v2)>tol
    elif isinstance(v1,(list,tuple)):
        v1 = array(v1,dtype=object).ravel()
        v2 = array(v2,dtype=object).ravel()
        for vv1,vv2 in zip(v1,v2):
            diff |= value_diff(vv1,vv2)
        #end for
    elif isinstance(v1,ndarray):
        v1 = v1.ravel()
        v2 = v2.ravel()
        for vv1,vv2 in zip(v1,v2):
            diff |= value_diff(vv1,vv2)
        #end for
    elif isinstance(v1,dict):
        diff = v1!=v2
    elif isinstance(v1,set):
        diff = v1!=v2
    elif v1 is None and v2 is None:
        diff = False
    else:
        diff = True # unsupported types
    #end if
    return diff
#end def value_diff


def object_diff(o1,o2,tol=1e-6,full=False):
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
        if value_diff(v1,v2,tol):
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



value_neq = value_diff
def value_eq(*args,**kwargs):
    return not value_neq(*args,**kwargs)
#end def value_eq

object_neq = object_diff
def object_eq(*args,**kwargs):
    return not object_neq(*args,**kwargs)
#end def object_eq



class NexusTestFail(Exception):
    None
#end class NexusTestFail

class NexusTestMisconstructed(Exception):
    None
#end class NexusTestMisconstructed

class NexusTestTripped(Exception):
    None
#end class NexusTestTripped


class NexusTestBase(object):
    nexus_test_dir = '.nexus_test' # ntest directory

    launch_path    = None # path from which ntest exe was launched
    current_test   = ''   # current NexusTest.name
    current_label  = ''   # current nlabel()
    test_count     =  0   # current test count
    label_count    =  0   # current label count in NexusTest.operation()
    current_assert =  0   # current assert count

    assert_trip    = -1   # developer tool to trip assert's one by one

    @staticmethod
    def assert_called():
        NexusTestBase.current_assert+=1
        ca = NexusTestBase.current_assert
        if ca==NexusTestBase.assert_trip:
            raise NexusTestTripped
        #end if
    #end def assert_called
        
#end class NexusTestBase




def nlabel(label):
    os.chdir(NexusTestBase.launch_path)
    NexusTestBase.current_label = label
    NexusTestBase.label_count  += 1
#end def nlabel


def nenter(path):
    test   = NexusTestBase.current_test
    label  = NexusTestBase.current_label
    tcount = str(NexusTestBase.test_count).zfill(2)
    lcount = str(NexusTestBase.label_count).zfill(2)
    test_dir  = '{0}_{1}'.format(tcount,test)
    label_dir = '{0}_{1}'.format(lcount,label)
    ntdir  = NexusTestBase.nexus_test_dir
    nlpath = NexusTestBase.launch_path
    if len(label)==0:
        enter_path = os.path.join(nlpath,ntdir,test_dir)
    else:
        enter_path = os.path.join(nlpath,ntdir,test_dir,label_dir)
    #end if
    os.makedirs(enter_path)
#end def nenter


def nleave():
    os.chdir(NexusTestBase.launch_path)
#end def nleave


def npass():
    None
#end def npass


def nfail(exception=NexusTestFail('Nexus test failed.')):
    raise exception
#end def nfail


def nassert(result):
    if not isinstance(result,bool):
        raise NexusTestMisconstructed
    elif not result:
        nfail()
    else:
        npass()
    #end if
    NexusTestBase.assert_called()
#end def nassert



class NexusTest(NexusTestBase):

    status_options = dict(
        unused = 0,
        passed = 1,
        failed = 2,
        )

    status_map = dict()
    for k,v in status_options.iteritems():
        status_map[v] = k
    #end for

    test_list = []
    test_dict = {}


    @staticmethod
    def setup():
        NexusTestBase.launch_path = os.path.getcwd()
    #end def setup


    def __init__(self,name,operation):
        if not isinstance(name,str):
            raise NexusTestMisconstructed
        #end if
        self.name         = name
        self.operation    = operation
        self.exception    = None
        self.status  = NexusTest.status_options['unused']
        NexusTest.test_list.append(self)
        NexusTest.test_dict[self.name] = self
    #end def __init__

    @property
    def unused(self):
        return self.status==NexusTest.status_options['unused']
    #end def unused

    @property
    def passed(self):
        return self.status==NexusTest.status_options['passed']
    #end def passed

    @property
    def failed(self):
        return self.status==NexusTest.status_options['failed']
    #end def failed

    def run(self):
        NexusTestBase.current_test  = self.name
        NexusTestBase.current_label = ''
        NexusTestBase.test_count   += 1
        NexusTestBase.label_count   = 0
        os.chdir(NexusTestBase.launch_path)
        try:
            self.operation()
            self.status=NexusTest.status_options['passed']
        except Exception,e:
            self.exception = e
            self.traceback = sys.exc_info()[2]
            self.status=NexusTest.status_options['failed']
        #end try
    #end def run

    def message(self):
        s = ''
        s+='test name    : {0}\n'.format(self.name)
        status = NexusTest.status_map[self.status]
        s+='test status  : {0}\n'.format(status)
        if self.failed and self.exception is not None:# and not isinstance(self.exception,NexusTestFail):
            if len(NexusTestBase.current_label)>0:
                s+='test sublabel: {0}\n'.format(NexusTestBase.current_label)
            #end if
            e = self.exception
            btrace = traceback.format_tb(self.traceback)
            if isinstance(e,NexusTestFail):
                btrace = btrace[:-1]
            elif isinstance(e,NexusTestMisconstructed):
                btrace = btrace[:-1]
                s+='exception    : Nexus test is misconstructed.  Please contact developers.\n'
            else:
                s+='exception    : "{0}"\n'.format(e.__class__.__name__+': '+str(e).replace('\n','\n              '))
            #end if
            s+='backtrace    :\n'
            for s2 in btrace: 
                s+=s2
            #end for
        #end if
        return s
    #end def message
#end class NexusTest
