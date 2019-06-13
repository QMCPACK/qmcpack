##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  pwscf_data_reader.py                                              #
#    Supports data processing of PWSCF's xml output.                 #
#                                                                    #
#  Content summary:                                                  #
#    read_qexml                                                      #
#      Read a data-file.xml and convert to nested object structure.  #
#                                                                    #
#    QEXML                                                           #
#      Class to represent an xml element.                            #
#                                                                    #
#    readval                                                         #
#      Function converts an attribute value string to numeric form.  #
#                                                                    #
#====================================================================#


import os
from numpy import array
from generic import obj
from developer import DevBase
from debug import *


class QEXML(DevBase):
    def __init__(self):
        self._value = []
    #end def __init__

    array_keys = set('type size columns len _value'.split())
    def finalize(self):
        keys = list(self.keys())
        enums = obj()
        for k in keys:
            v = self[k]
            if isinstance(k,str) and '.' in k and k.split('.',1)[1].isdigit():
                name,index = k.split('.',1)
                index = int(index)
                if not name in enums:
                    enums[name] = QEXML()
                #end if
                enums[name][index] = v
                del self[k]
                continue
            #end if
            if isinstance(v,QEXML):
                if len(set(v.keys())-self.array_keys)==0:
                    a = array(v._value)
                    if len(a)==1:
                        a = a[0]
                    elif 'columns' in v and v.size%v.columns==0:
                        a.shape = v.size/v.columns,v.columns
                    #end if
                    self[k] = a
                else:
                    v.finalize()
                #end if
            #end if
        #end for
        for k,v in enums.iteritems():
            self[k] = v
            v.finalize()
        #end for
        if len(self._value)==0:
            del self._value
        elif not 'value' in self:
            if len(self._value)==1:
                self.value = self._value[0]
            else:
                self.value = array(self._value)
            #end if
            del self._value
        else:
            if len(self._value)==1:
                self._value = self._value[0]
            else:
                self._value = array(self._value)
            #end if
        #end if
    #end def finalize
#end class QEXML



bools = obj(F=False,T=True)
def readval(s):
    s = s.strip()
    v = None
    if s in bools:
        v = bools[s]
    else:
        try:
            v = int(s)
        except:
            try:
                v = float(s)
            except:
                if not ' ' in s:
                    v=s
                else:
                    s = s.split()
                    try:
                        v=list(array(s,dtype=int))
                    except:
                        try:
                            v=list(array(s,dtype=float))
                        except:
                            v=s
                        #end try
                    #end try
                #end if
            #end try
        #end try
    #end if
    return v
#end readval
                    
        


def read_qexml(inp):
    if isinstance(inp,list):
        rawlines = inp
    elif isinstance(inp,str):# and os.path.exists(inp):
        rawlines = open(inp,'r').read().splitlines()
    else:
        print 'read_qexml error: input can only be filename or list of lines'
        print '  input received: ',inp
        return None
    #end if

    lines = []
    incomment = False
    for line in rawlines:
        ls = line.strip()
        if ls.startswith('<!--'):
            incomment = True
        #end if
        if ls.endswith('-->'):
            incomment = False
            continue
        #end if
        if not ls.startswith('<?') and not incomment and ls!='':
            lines.append(ls)
        #end if
    #end for

    root = QEXML()
    stack = [root]
    name = ''
    level = -1
    if len(lines)>0:
        for line in lines:
            if line.startswith('</'):
                #print level*' '+'end',line.strip('</>').lower()
                level-=1
                stack.pop()
                cur = stack[-1]
            elif line.startswith('<'):
                level+=1
                base = stack[level]
                ls = line.strip('<>')
                ended = ls.endswith('/')
                if ended:
                    ls = ls[0:-1]
                #end if
                instr=False
                lsl = list(tuple(ls))
                ls = ''
                for i in range(len(lsl)):
                    c = lsl[i]
                    if c=='"':
                        instr = not instr
                    elif c==' ' and instr:
                        c='@'
                    #end if
                    ls+=c
                #end for
                tokens = ls.split()
                name = tokens[0].lower().replace('-','_')
                attrs = tokens[1:]                
                cur = QEXML()
                base[name] = cur
                for attr in attrs:
                    n,v = attr.replace('@',' ').split('=',1)
                    v=readval(v.strip('"'))
                    cur[n.lower()]=v
                #end for
                if not ended:
                    #print level*' '+name
                    stack.append(cur)
                else:
                    #print level*' '+name+'  end',name
                    level-=1
                    cur = base
                #end if
            else:
                v = readval(line)
                #print level*' ',v
                if isinstance(v,list):
                    cur._value.extend(v)
                else:
                    cur._value.append(v)
                #end if
            #end if
        #end for
    #end if

    root.finalize()

    return root
#end def read_qexml
