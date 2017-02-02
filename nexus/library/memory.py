##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################

######################################################################
# The following is adapted from
# http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/286222
# Python Cookbook, by Jean Brouwers
######################################################################

#====================================================================#
#  memory.py                                                         #
#    Calculate memory used by the Nexus host process.                #
#                                                                    #
#  Content summary:                                                  #
#    memory                                                          #
#      Return memory usage of the current process.                   #
#                                                                    #
#    resident                                                        #
#      Return resident memory usage.                                 #
#                                                                    #
#    stacksize                                                       #
#      Return stack size.                                            #
#                                                                    #
#====================================================================#


import os

_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
          'KB': 1024.0, 'MB': 1024.0*1024.0}

def _VmB(VmKey, pid=None):
    '''Private.
    '''
    global _scale
     # get pseudo file  /proc/<pid>/status
 
    if not pid:
        pid = os.getpid()
    proc_status = '/proc/%d/status' % pid
    try:
        t = open(proc_status)
        v = t.read()
        t.close()
    except:
        return 0.0  # non-Linux?
     # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    try:
        i = v.index(VmKey)
        v = v[i:].split(None, 3)  # whitespace
    except ValueError:
        return 0.0

    if len(v) < 3:
        return 0.0  # invalid format?
     # convert Vm value to bytes
    return float(v[1]) * _scale[v[2]]

def get_children(pid):
    proc_children = '/proc/%d/task/%d/children'%(pid,pid)
    try:
        t = open(proc_children,'r')
        v = t.read()
        t.close()
    except:
        return []

    children = [int(c) for c in v.split()]
    return children


def all_children(pid=None):
    if not pid:
        pid = os.getpid()

    child_list = get_children(pid)
    all_list = child_list[:]
    for child_pid in child_list:
        child2 = all_children(child_pid)
        all_list.extend(child2)
    return all_list


def memory(since=0.0, children=False):
    '''Return memory usage in bytes.
    '''
    mem = 0.0
    if children:
        for child_pid in all_children():
            mem += _VmB('VmSize:',pid=child_pid)
    mem += _VmB('VmSize:') - since
    return mem


def resident(since=0.0, children=False):
    '''Return resident memory usage in bytes.
    '''
    mem = 0.0
    if children:
        for child_pid in all_children():
            mem += _VmB('VmRSS:',pid=child_pid)

    mem += _VmB('VmRSS:') - since
    return mem
    #return _VmB('VmRSS:') - since


def stacksize(since=0.0, children=False):
    '''Return stack size in bytes.
    '''
    mem = 0.0
    if children:
        for child_pid in all_children():
            mem += _VmB('VmStk:',pid=child_pid)
    mem += _VmB('VmStk:') - since
    return mem
