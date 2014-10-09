
import os
from numpy import array,zeros
from generic import obj
from simulation import Simulation,SimulationAnalyzer
from vasp_input import VaspInput
from developer import DevBase
from debug import *



class VXML(DevBase):
    basic_types = set('i v dimension field set time'.split())

    data_types = obj(int=int,string=str,float=float)
    
    def __init__(self,tag,attr=None):
        self._tag   = tag
        self._lines = []
        self._attr  = None
        self._value = None
        if attr!=None and len(attr)>0:
            self._attr = obj()
            tokens = attr.split('"')
            next=0
            for token in tokens:
                next+=1
                if len(token)>0 and token[-1]=='=':
                    name = token.strip()[:-1]
                    self._attr[name]=tokens[next].replace(' ','_').replace('-','_').lower()
                #end if
            #end for
        #end if
    #end def __init__

    def _is_empty(self):
        return len(self)-4==0
    #end def _is_empty

    def _add(self,new):
        tag = new._tag
        if not tag in self:
            self[tag] = new
        else:
            cur = self[tag]
            if 0 in cur:
                cur._append(new)
            else:
                coll = VXMLcoll(tag)
                coll._append(cur)
                coll._append(new)
                self[tag] = coll
            #end if
        #end if
    #end def _add

    def _parse(self):
        # rename sub objects if name is present
        for name in list(self.keys()):
            value = self[name]
            if isinstance(value,VXML):
                if value._attr!=None and 'name' in value._attr:
                    del self[name]
                    self[value._attr.name] = value
                    del value._attr.name
                elif isinstance(value,VXMLcoll):
                    for n in list(value.keys()):
                        v = value[n]
                        if isinstance(v,VXML):
                            if v._attr!=None and 'name' in v._attr:
                                del value[n]
                                self[v._attr.name] = v
                                del v._attr.name
                            #end if
                        #end if
                    #end for
                    if len(value)==0:
                        del self[name]
                    else:
                        value._reorder()
                    #end if
                #end if
            #end if
        #end for
        # have all sub objects parse (read fields, set _value)
        for v in self:
            if isinstance(v,VXML):
                v._parse()
            #end if
        #end for
        lines = self._lines
        if len(lines)>0:
            #fail = False
            #try:
            if self._tag=='array':
                self._parse_array(lines)
            else:
                self._parse_values(lines)
            #end if
            #except Exception,e:
            #    print '_parse failed!'
            #    print e
            #    print self
            #    print tag
            #    print attr
            #    print lines[0:min(10,len(lines))]
            #    print
            #    print
            #    fail = True
            ##end try
            #if fail:
            #    self.error('parse failed please see information above','read_vxml')
            ##end if
        #end if

        # if sub-objects resolve to a value, replace with that value
        for name in list(self.keys()):
            value = self[name]
            if isinstance(value,VXML) and value._value!=None: 
                self[name] = value._value
            #end if
        #end for

        # assign attributes
        if self._attr!=None:
            if 'type' in self._attr:
                del self._attr.type
            #end if
            self.transfer_from(self._attr)
        #end if

        return
    #end def _parse


    def _parse_values(self,lines):
        if len(lines)==1 and not '<' in lines[0]:
            self._value = readval(lines[0])
        else:
            arr = None
            for line in lines:
                start = line.startswith('<') and not line.startswith('</')
                end   = line.endswith('>')
                if start:
                    fline = line
                else:
                    fline += line
                #end if
                if end:
                    tokens = fline.replace('<','|').replace('>','|').split('|')
                    tn = tokens[1]
                    val = readval(tokens[2].strip())
                    if 'name=' in tn:
                        name = tn.split('name="',1)[1].split('"')[0].lower()
                        self[name] = val
                    elif arr is None:
                        arr = [val]
                    else:
                        arr.append(val)
                    #end if
                #end if
            #end for
            if arr!=None:
                self._value = array(arr)
            #end if
        #end if
    #end def parse_values


    def _parse_array(self,lines):
        #print 'parsing array'
        dims       = obj()
        fields     = obj()
        dim_counts = None
        field_list = []
        level      = -1
        set_dims   = False
        for line in lines:
            if line.startswith('<'):
                if line.startswith('<dimension'):
                    tokens = line.replace('<','|').replace('>','|').split('|')
                    tn = tokens[1]
                    dname = tokens[2].lower().replace(' ','_').replace('-','_')
                    if 'dim=' in tn:
                        d = int(tn.split('dim="',1)[1].split('"')[0])
                        dims[d] = dname
                    else:
                        dims.append(dname)
                    #end if
                elif line.startswith('<field'):
                    tokens = line.replace('<','|').replace('>','|').split('|')
                    tn = tokens[1]
                    fname = tokens[2].lower().replace(' ','_').replace('-','_')
                    if 'type=' in tn:
                        t = tn.split('type="',1)[1].split('"')[0]
                        if t in VXML.data_types:
                            dtype = VXML.data_types[t]
                        else:
                            self.error('field type {0} is unrecognized: {1}'.format(t,line))
                        #end if
                    else:
                        dtype = float
                    #end if
                    fields.append(obj(name=fname,dtype=dtype))
                elif line.startswith('<set'):
                    if not set_dims:
                        dims   = dims.list()
                        dims.reverse()
                        dims = tuple(dims)
                        dim_counts = zeros((len(dims),),dtype=int)
                        set_dims = True
                    #end if
                    level += 1
                    dim_counts[level]=0
                elif line.startswith('</set'):
                    level -= 1
                    if level!=-1:
                        dim_counts[level]+=1
                    #end if
                else:
                    self.error('array parsing failed\n unrecognized xml encountered: {0}'.format(line),'read_vxml')
                #end if
            else:
                dim_counts[level]+=1
                field_list.append(line.split())
            #end if
        #end for
        self.dims   = dims
        for findex,field in fields.iteritems():
            lst = []
            for field_vals in field_list:
                lst.append(field_vals[findex])
            #end for
            arr = array(lst,dtype=field.dtype).ravel()
            arr.shape = tuple(dim_counts)
            self[field.name] = arr
        #end for
        #print '  done'
    #end def _parse_array


    def _remove_empty(self):
        for n in list(self.keys()):
            v = self[n]
            if isinstance(v,VXML):
                v._remove_empty()
                if isinstance(v,VXMLcoll) and v._is_empty():
                    del self[n]
                #end if
            #end if
        #end for 
    #end def _remove_empty


    def _remove_hidden(self):
        del self._tag
        del self._attr
        del self._lines
        del self._value
        for v in self:
            if isinstance(v,VXML):
                v._remove_hidden()
            #end if
        #end for
    #end def _remove_hidden()
#end class VXML


class VXMLcoll(VXML):
    def _append(self,new):
        index = len(self)-4
        self[index]=new
    #end def _append

    def _reorder(self):
        n=0
        for key in sorted(self.keys()):
            value = self[key]
            if isinstance(value,VXML):
                del self[key]
                self[n]=value
                n+=1
            #end if
        #end for
    #end def _reorder
#end class VXMLcoll


booldict = dict(T=True,F=False)

def readval(val):
    fail = False
    split = False
    if isinstance(val,str):
        split = ' ' in val
    #end if
    if isinstance(val,list) or split:
        if split:
            val = val.split()
        #end if
        try:
            v = array(val,dtype=int)
        except:
            try:
                v = array(val,dtype=float)
            except:
                try:
                    v = array(val,dtype=str)
                except:
                    fail = True
                #end try
            #end try
        #end try
    elif val in booldict:
        v = booldict[val]
    else:
        try:
            v = int(val)
        except:
            try:
                v = float(val)
            except:
                v = val
            #end try
        #end try
    #end if
    if fail:
        VXML.class_error('failed to read value: "{0}"'.format(val),'read_vxml')
    #end if
    return v
#end def readval
            


def read_vxml(filepath):
    if not os.path.exists(filepath):
        VXML.class_error('file {0} does not exist'.format(filepath),'read_vxml')
    #end if
    #print 'read'
    contents = open(filepath,'r').read()
    #print 'replace'
    contents = contents.replace('<rc>',' ').replace('</rc>',' ')
    contents = contents.replace('<r>' ,' ').replace('</r>' ,' ')
    contents = contents.replace('<c>' ,' ').replace('</c>' ,' ')
    #print 'split lines'
    lines    = contents.splitlines()

    #print 'process lines'
    root = VXML('vasprun')
    stack = [root]
    cur = stack[0]
    for line in lines:
        ls = line.strip()
        if ls.startswith('</'):
            tag = ls[2:-1]
            if tag==cur._tag:
                stack.pop()
                cur = stack[-1]
                #print len(stack)*'  '+'end '+tag
            else:
                cur._lines.append(ls)
            #end if
        elif ls.startswith('<?'):
            None
        elif ls.startswith('<'):
            ta,rest = ls[1:].split('>',1)
            tokens = ta.split(' ',1)
            tag = tokens[0]
            if not tag in VXML.basic_types:
                if len(tokens)==1:
                    attr = None
                else:
                    attr = tokens[1].strip()
                #end if
                if ls.endswith('</{0}>'.format(tag)):
                    new = VXML(tag,attr)
                    new._lines.append(ls.replace('<','|').replace('>','|').split('|')[2])
                    cur._add(new)
                    #print len(stack)*'  '+'end '+tag
                else:
                    #print len(stack)*'  '+tag
                    new = VXML(tag,attr)
                    cur._add(new)
                    cur = new
                    stack.append(cur)
                #end if
            else:
                cur._lines.append(ls)
            #end if
        else:
            cur._lines.append(ls)
        #end if
    #end for

    if len(stack)!=1:
        VXML.class_error('read failed\nxml tree did not seem to close')
    #end if

    #print 'parse'
    root._parse()
    root._remove_empty()
    root._remove_hidden()

    #print 'done'
    return root
#end def read_vxml



class VaspAnalyzer(SimulationAnalyzer):
    def __init__(self,arg0=None,xml=False,analyze=False):
        self.info = obj(xml=xml)
        prefix = None
        if isinstance(arg0,Simulation):
            sim = arg0
            infile = sim.infile
            path   = sim.locdir
        elif arg0!=None:
            path,infile = os.path.split(arg0)
            if infile=='':
                infile = None
            #end if
            if infile!=None:
                if not infile.endswith('INCAR'):
                    self.error('please provide the path to an INCAR file')
                #end if
                prefix = infile.replace('INCAR','').strip()
                if prefix=='':
                    prefix=None
                #end if
            #end if
        else:
            self.info.xml = False
            return
        #end if
        self.info.set(
            path   = path,
            infile = infile,
            prefix = prefix
            )
        if analyze:
            self.analyze()
        #end if
    #end def __init__

    def analyze(self):
        if self.info.xml:
            xmlfile = 'vasprun.xml'
            if self.info.prefix!=None:
                xmlfile = self.info.prefix+xmlfile
            #end if
            self.xmldata = read_vxml(os.path.join(self.info.path,xmlfile))
        #end if
    #end def analyze
#end class VaspAnalyzer


