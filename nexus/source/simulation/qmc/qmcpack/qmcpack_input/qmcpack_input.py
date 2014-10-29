import os
import types as pytypes
import keyword
from numpy import fromstring,empty,array,float64,\
    loadtxt,ndarray,dtype,sqrt,pi,arange,exp,eye,\
    ceil,mod,dot,abs,identity
from StringIO import StringIO
from superstring import string2val
from generic import obj
from xmlreader import XMLreader,XMLelement
from project_base import Pobj
from structure import Structure
from physical_system import PhysicalSystem
from simulation import SimulationInput
from pwscf_input import array_to_string as pwscf_array_string
from debug import *


yesno_dict     = {True:'yes' ,False:'no'}
truefalse_dict = {True:'true',False:'false'}
onezero_dict   = {True:'1'   ,False:'0'}
boolmap={'yes':True,'no':False,'true':True,'false':False,'1':True,'0':False}

def is_int(var):
    try:
        int(var)
        return True
    except ValueError:
        return False
#end def is_int

def is_float(var):
    try:
        float(var)
        return True
    except ValueError:
        return False
#end def is_float

def is_array(var,type):
    try:
        if isinstance(var,str):
            array(var.split(),type)
        else:
            array(var,type)
        #end if
        return True
    except ValueError:
        return False
#end def is_float_array


def attribute_to_value(attr):
    if is_int(attr):
        val = int(attr)
    elif is_float(attr):
        val = float(attr)
    elif is_array(attr,int):
        val = array(attr.split(),int)
        if val.size==9:
            val.shape = 3,3
        #end if
    elif is_array(attr,float):
        val = array(attr.split(),float)
    else:
        val = attr
    #end if
    return val
#end def attribute_to_value



#local write types
def yesno(var):
    if var:
        return 'yes'
    else:
        return 'no'
    #end if
#end def yesno

def onezero(var):
    if var:
        return '1'
    else:
        return '0'
    #end if
#end def onezero

def truefalse(var):
    if var:
        return 'true'
    else:
        return 'false'
    #end if
#end def onezero

bool_write_types = set([yesno,onezero,truefalse])




class QIobj(Pobj):
    None
#end class QIobj


class meta(obj):
    None
#end class meta


class section(QIobj):
    def __init__(self,*args,**kwargs):
        self.args   = args
        self.kwargs = kwargs
    #end def __init__
#end class section


class collection(QIobj):
    def __init__(self,*elements,**kwargs):
        if len(elements)==1 and isinstance(elements[0],list):
            elements = elements[0]
        #end if
        for element in elements:
            identifier = element.identifier
            if isinstance(identifier,str):
                key = element[identifier]
            else:
                key = ''
                for ident in identifier:
                    if ident in element:
                        key+=element[ident]
                    #end if
                #end for
            #end if
            self[key] = element
        #end for
        for key,element in kwargs.iteritems():
            self[key] = element
        #end for
    #end def __init__

    def get_single(self,preference):
        if len(self)>0:
            if preference in self:
                return self[preference]
            else:
                return self.list()[0]
            #end if
        else:
            return self
        #end if
    #end def get_single
#end class collection
def make_collection(elements):
    if len(elements)>0 and 'identifier' in elements[0].__class__.__dict__.keys():
        return collection(*elements)
    else:
        coll = collection()
        for i in range(len(elements)):
            coll[i] = elements[i]
        #end for
        return coll
        #return elements
    #end if
#end def make_collection


class classcollection(QIobj):
    def __init__(self,*classes):
        if len(classes)==1 and isinstance(classes[0],list):
            classes = classes[0]
        #end if
        self.classes = classes
    #end def __init__
#end class classcollection


class Names(QIobj):
    condensed_names = obj()
    expanded_names = None

    escape_names = set(keyword.kwlist)
    escaped_names = list(escape_names)
    for i in range(len(escaped_names)):
        escaped_names[i]+='_'
    #end for
    escaped_names = set(escaped_names)

    @staticmethod
    def set_expanded_names(**kwargs):
        Names.expanded_names = obj(**kwargs)
    #end def set_expanded_names


    def expand_name(self,condensed):
        expanded = condensed
        cname = self.condense_name(condensed)
        if cname in self.escaped_names:
            cname = cname[:-1]
            expanded = cname
        #end if
        if cname in self.expanded_names:
            expanded = self.expanded_names[cname]
        #end if
        return expanded
    #end def expand_name

    def condense_name(self,expanded):
        condensed = expanded
        condensed = condensed.replace('___','_').replace('__','_')
        condensed = condensed.replace('-','_').replace(' ','_')
        condensed = condensed.lower()
        if condensed in self.escape_names:
            condensed += '_'
        #end if
        self.condensed_names[expanded]=condensed
        return condensed
    #end def condense_name

    def condense_names(self,*namelists):
        out = []
        for namelist in namelists:
            exp = obj()
            for expanded in namelist:
                condensed = self.condense_name(expanded)
                exp[condensed]=expanded
            #end for
            out.append(exp)
        #end for
        return out
    #end def condense_names

    def condensed_name_report(self):
        print
        print 'Condensed Name Report:'
        print '----------------------'
        keylist = array(self.condensed_names.keys())
        order = array(self.condensed_names.values()).argsort()
        keylist = keylist[order]
        for expanded in keylist: 
            condensed = self.condensed_names[expanded]
            if expanded!=condensed:
                print "    {0:15} = '{1}'".format(condensed,expanded)
            #end if
        #end for
        print
        print
    #end def condensed_name_report
#end class Names




class QIxml(Names):

    def init_from_args(self,args):
        print
        print args
        print
        self.not_implemented()
    #end def init_from_args



    @classmethod
    def init_class(cls):
        vars = cls.__dict__.keys()
        init_vars = dict(tag         = cls.__name__,
                         attributes  = [],
                         elements    = [],
                         text        = None,
                         parameters  = [],
                         attribs     = [],
                         costs       = [],
                         h5tags      = [],
                         types       = obj(),
                         write_types = obj(),
                         defaults    = obj()
                         )
        for var,val in init_vars.iteritems():
            if not var in vars:
                cls.__dict__[var] = val
            #end if
        #end for
        for v in ['attributes','elements','parameters','attribs','costs','h5tags']:
            names = cls.__dict__[v]
            for i in range(len(names)):
                if names[i] in cls.escape_names:
                    names[i]+='_'
                #end if
            #end for
        #end for
        cls.params = cls.parameters + cls.attribs + cls.costs + cls.h5tags
        cls.plurals_inv = obj()
        for e in cls.elements:
            if e in plurals_inv:
                cls.plurals_inv[e] = plurals_inv[e]
            #end if
        #end for
        cls.plurals = cls.plurals_inv.inverse()
    #end def init_class


    def write(self,indent_level=0,pad='   ',first=False):
        self.check_junk(exit=True)
        indent  = indent_level*pad
        ip = indent+pad
        ipp= ip+pad
        expanded_tag = self.expand_name(self.tag)
        c = indent+'<'+expanded_tag
        for a in self.attributes:
            if a in self:
                val = self[a]
                if isinstance(val,str):
                    val = self.expand_name(val)
                #end if
                c += ' '+self.expand_name(a)+'='
                if a in self.write_types:
                    c += '"'+self.write_types[a](val)+'"'
                else:
                    c += '"'+param.write(val)+'"'
                #end if
            #end if
        #end for
        #if first:
        #    c+=' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://www.mcc.uiuc.edu/qmc/schema/molecu.xsd"'
        ##end if
        #no_contents = len(set(self.keys())-set(self.elements)-set(self.plurals.keys()))==0
        no_contents = len(set(self.keys())-set(self.attributes))==0
        if no_contents:
            c += '/>\n'
        else:
            c += '>\n'
            for v in self.h5tags:
                if v in self:
                    if v in self.write_types:
                        write_type = self.write_types[v]
                    else:
                        write_type = None
                    #end if
                    c += param.write(self[v],name=self.expand_name(v),tag='h5tag',mode='elem',pad=ip,write_type=write_type)
                #end if
            #end for
            for v in self.costs:
                if v in self:
                    c += param.write(self[v],name=self.expand_name(v),tag='cost',mode='elem',pad=ip)
                #end if
            #end for
            for p in self.parameters:
                if p in self:
                    if p in self.write_types:
                        write_type = self.write_types[p]
                    else:
                        write_type = None
                    #end if
                    c += param.write(self[p],name=self.expand_name(p),mode='elem',pad=ip,write_type=write_type)
                #end if
            #end for
            for a in self.attribs:
                if a in self:
                    if a in self.write_types:
                        write_type = self.write_types[a]
                    else:
                        write_type = None
                    #end if
                    c += param.write(self[a],name=self.expand_name(a),tag='attrib',mode='elem',pad=ip,write_type=write_type)
                #end if
            #end for
            for e in self.elements:
                if e in self:
                    elem = self[e]
                    if isinstance(elem,QIxml):
                        c += self[e].write(indent_level+1)
                    else:
                        begin = '<'+e+'>'                        
                        contents = param.write(elem)
                        end = '</'+e+'>'
                        if contents.strip()=='':
                            c += ip+begin+end+'\n'
                        else:                            
                            c += ip+begin+'\n'
                            c += ipp+contents+'\n'
                            c += ip+end+'\n'
                        #end if
                    #end if
                elif e in plurals_inv and plurals_inv[e] in self:
                    coll = self[plurals_inv[e]]
                    coll_len = len(coll)
                    if 0 in coll and coll_len-1 in coll:
                        for i in range(coll_len):
                            instance = coll[i]
                            c += instance.write(indent_level+1)
                        #end for
                    else:
                        for instance in coll:
                            c += instance.write(indent_level+1)
                        #end for
                    #end if
                #end if
            #end for
            if self.text!=None:
                c = c.rstrip('\n')
                c+=param.write(self[self.text],mode='elem',pad=ip,tag=None,normal_elem=True)
            #end if
            c+=indent+'</'+expanded_tag+'>\n'
        #end if
        return c
    #end def write


    def __init__(self,*args,**kwargs):
        if Param.metadata==None:
            Param.metadata = meta()
        #end if
        if len(args)==1:
            a = args[0]
            if isinstance(a,XMLelement):
                self.init_from_xml(a)
            elif isinstance(a,section):
                self.init_from_inputs(a.args,a.kwargs)
            elif isinstance(a,self.__class__):
                self.transfer_from(a)
            else:
                self.init_from_inputs(args,kwargs)
            #end if
        else:
            self.init_from_inputs(args,kwargs)
        #end if
    #end def __init__


    def init_from_xml(self,xml):
        al,el = self.condense_names(xml._attributes.keys(),xml._elements.keys())
        xa,sa = set(al.keys()) , set(self.attributes)
        attr = xa & sa
        junk = xa-attr
        junk_elem = []
        for e,ecap in el.iteritems():
            value = xml._elements[ecap]
            if (isinstance(value,list) or isinstance(value,tuple)) and e in self.plurals_inv.keys():
                p = self.plurals_inv[e]
                plist = []
                for instance in value:
                    plist.append(types[e](instance))
                #end for
                self[p] = make_collection(plist)
            elif e in self.elements:
                self[e] = types[e](value)
            elif e in ['parameter','attrib','cost','h5tag']:
                if isinstance(value,XMLelement):
                    value = [value]
                #end if
                for p in value:
                    name = self.condense_name(p.name)
                    if name in self.params:
                        self[name] = param(p)
                    else:
                        junk_elem.append(name)
                    #end if
                #end for
            else:
                junk_elem.append(e)
            #end if
        #end for
        junk = junk | set(junk_elem)
        self.check_junk(junk)
        for a in attr:
            if a in self.write_types and self.write_types[a] in bool_write_types:
                aval = xml._attributes[al[a]]
                if aval in boolmap:
                    self[a] = boolmap[aval]
                else:
                    self.error('{0} is not a valid value for boolean attribute {1}\n  valid values are: {2}'.format(aval,a,boolmap.keys()))
                #end if
            else:
                self[a] = attribute_to_value(xml._attributes[al[a]])
            #end if
        #end for
        if self.text!=None:
            self[self.text] = param(xml)
        #end if
    #end def init_from_xml


    def init_from_inputs(self,args,kwargs):
        if len(args)>0:
            if len(args)==1 and isinstance(args[0],self.__class__):
                self.transfer_from(args[0])
            else:
                self.init_from_args(args)
            #end if
        #end if
        self.init_from_kwargs(kwargs)
    #end def init_from_inputs


    def init_from_kwargs(self,kwargs):
        ks=[]
        kmap = dict()
        for key,val in kwargs.iteritems():
            ckey = self.condense_name(key)
            ks.append(ckey)
            kmap[ckey] = val
        #end for
        ks = set(ks)
        kwargs = kmap
        h5tags     = ks & set(self.h5tags)
        costs      = ks & set(self.costs)
        parameters = ks & set(self.parameters)
        attribs    = ks & set(self.attribs)
        attr = ks & set(self.attributes)
        elem = ks & set(self.elements)
        plur = ks & set(self.plurals.keys())
        if self.text!=None:
            text = ks & set([self.text])
        else:
            text = set()
        #end if
        junk = ks -attr -elem -plur -h5tags -costs -parameters -attribs -text
        self.check_junk(junk)

        for v in h5tags:
            self[v] = param(kwargs[v])
        #end for
        for v in costs:
            self[v] = param(kwargs[v])
        #end for
        for v in parameters:
            self[v] = param(kwargs[v])
        #end for
        for v in attribs:
            self[v] = param(kwargs[v])
        #end for
        for a in attr:
            self[a] = kwargs[a]
        #end for
        for e in elem:
            self[e] = types[e](kwargs[e])
        #end for
        for p in plur:
            plist = []
            e = self.plurals[p]
            kwcoll = kwargs[p]
            if isinstance(kwcoll,collection):
                cobj = collection()
                for name,instance in kwcoll.iteritems():
                    iobj = types[e](instance)
                    if isinstance(iobj.identifier,str):
                        iobj[iobj.identifier]=name
                    #end if
                    cobj[name] = iobj
                #end for 
                self[p] = cobj
            else:
                for instance in kwargs[p]:
                    plist.append(types[e](instance))
                #end for
                self[p] = make_collection(plist)
            #end if
        #end for
        for t in text:
            self[t] = kwargs[t]
        #end for
    #end def init_from_kwargs


    def incorporate_defaults(self,elements=False,overwrite=False,propagate=True):
        for name,value in self.defaults.iteritems():
            valtype = type(value)
            defval=None
            if isinstance(value,classcollection):
                if elements:
                    coll=[]
                    for cls in value.classes:
                        ins = cls()
                        ins.incorporate_defaults()
                        coll.append(ins)
                    #end for
                    defval = make_collection(coll)
                #end if
            elif valtype==pytypes.ClassType:
                if elements:
                    defval = value()
                #end if
            elif valtype==pytypes.FunctionType:
                defval = value()
            else:
                defval = value
            #end if
            if defval!=None:
                if overwrite or not name in self:
                    self[name] = defval
                #end if
            #end if
        #end for
        if propagate:
            for name,value in self.iteritems():
                if isinstance(value,QIxml):
                    value.incorporate_defaults(elements,overwrite)
                elif isinstance(value,collection):
                    for v in value:
                        if isinstance(v,QIxml):
                            v.incorporate_defaults(elements,overwrite)
                        #end if
                    #end for
                #end if
            #end for
        #end if
    #end def incorporate_defaults                    


    def check_junk(self,junk=None,exit=False):
        if junk is None:
            ks = set(self.keys())
            h5tags     = ks & set(self.h5tags)
            costs      = ks & set(self.costs)
            parameters = ks & set(self.parameters)
            attribs    = ks & set(self.attribs)
            attr = ks & set(self.attributes)
            elem = ks & set(self.elements)
            plur = ks & set(self.plurals.keys())
            if self.text!=None:
                text = ks & set([self.text])
            else:
                text = set()
            #end if
            junk = ks -attr -elem -plur -h5tags -costs -parameters -attribs -text
        #end if
        if len(junk)>0:
            oname = ''
            if self.tag!=self.__class__.__name__:
                oname = ' ('+self.__class__.__name__+')'
            #end if
            msg = '{0}{1} does not have the following attributes/elements:\n'.format(self.tag,oname)
            for jname in junk:
                msg+='    '+jname+'\n'
            #end for
            self.error(msg,'QmcpackInput',exit=exit,trace=exit)
        #end if
    #end def check_junk


    def get_single(self,preference):
        return self
    #end def get_single


    def get(self,names,namedict=None,host=False,root=True):
        if namedict is None:
            namedict = {}
        #end if
        if isinstance(names,str):
            names = [names]
        #end if
        if root and not host:
            if 'identifier' in self.__class__.__dict__ and self.identifier in self:
                identity = self[self.identifier]
            else:
                identity = None
            #end if
            for name in names:
                if name==self.tag:
                    namedict[name]=self
                elif name==identity:
                    namedict[name]=self
                #end if
            #end for
        #end if
        for name in names:
            loc = None
            if name in self:
                loc = name
            elif name in plurals_inv and plurals_inv[name] in self:
                loc = plurals_inv[name]
            #end if
            name_absent = not name in namedict 
            not_element = False
            if not name_absent:
                not_xml  = not isinstance(namedict[name],QIxml)
                not_coll = not isinstance(namedict[name],collection)
                not_element = not_xml and not_coll
            #end if
            if loc!=None and (name_absent or not_element):
                if host:
                    namedict[name] = self
                else:
                    namedict[name] = self[loc]
                #end if
            #end if
        #end for
        for name,value in self.iteritems():
            if isinstance(value,QIxml):
                value.get(names,namedict,host,root=False)
            elif isinstance(value,collection):
                for n,v in value.iteritems():
                    name_absent = not n in namedict 
                    not_element = False
                    if not name_absent:
                        not_xml  = not isinstance(namedict[n],QIxml)
                        not_coll = not isinstance(namedict[n],collection)
                        not_element = not_xml and not_coll
                    #end if
                    if n in names and (name_absent or not_element):
                        if host:
                            namedict[n] = value
                        else:
                            namedict[n] = v
                        #end if
                    #end if
                    if isinstance(v,QIxml):
                        v.get(names,namedict,host,root=False)
                    #end if
                #end if
            #end if
        #end for
        if root:
            namelist = []
            for name in names:
                if name in namedict:
                    namelist.append(namedict[name])
                else:
                    namelist.append(None)
                #end if
            #end for
            if len(namelist)==1:
                return namelist[0]
            else:
                return namelist
            #end if
        #end if
    #end def get

    def remove(self,*names):
        if len(names)==1 and not isinstance(names[0],str):
            names = names[0]
        #end if
        remove = []
        for name in names:
            attempt = True
            if name in self:
                rname = name
            elif name in plurals_inv and plurals_inv[name] in self:
                rname = plurals_inv[name]
            else:
                attempt = False
            #end if
            if attempt:
                val = self[rname]
                if isinstance(val,QIxml) or isinstance(val,collection):
                    remove.append(rname)
                #end if
            #end if
        #end for
        for name in remove:
            del self[name]
        #end for
        for name,value in self.iteritems():
            if isinstance(value,QIxml):
                value.remove(*names)
            elif isinstance(value,collection):
                for element in value:
                    if isinstance(element,QIxml):
                        element.remove(*names)
                    #end if
                #end for
            #end if
        #end for
    #end def remove


    def assign(self,**kwargs):
        for var,vnew in kwargs.iteritems():
            if var in self:
                val = self[var]
                not_coll = not isinstance(val,collection)
                not_xml  = not isinstance(val,QIxml)
                not_arr  = not isinstance(val,ndarray)
                if not_coll and not_xml and not_arr:
                    self[var] = vnew
                #end if
            #end if
        #end for
        for vname,val in self.iteritems():
            if isinstance(val,QIxml):
                val.assign(**kwargs)
            elif isinstance(val,collection):
                for v in val:
                    if isinstance(v,QIxml):
                        v.assign(**kwargs)
                    #end if
                #end for
            #end if
        #end for
    #end def assign


    def replace(self,*args,**kwargs):
        if len(args)==2 and isinstance(args[0],str) and isinstance(args[1],str):
            vold,vnew = args
            args = [(vold,vnew)]
        #end for
        for valpair in args:
            vold,vnew = valpair
            for var,val in self.iteritems():
                not_coll = not isinstance(val,collection)
                not_xml  = not isinstance(val,QIxml)
                not_arr  = not isinstance(val,ndarray)
                if not_coll and not_xml and not_arr and val==vold:
                    self[var] = vnew
                #end if
            #end for
        #end for
        for var,valpair in kwargs.iteritems():
            vold,vnew = valpair
            if var in self:
                val = self[var]
                if vold==None:
                    self[var] = vnew
                else:
                    not_coll = not isinstance(val,collection)
                    not_xml  = not isinstance(val,QIxml)
                    not_arr  = not isinstance(val,ndarray)
                    if not_coll and not_xml and not_arr and val==vold:
                        self[var] = vnew
                    #end if
                #end if
            #end if
        #end for
        for vname,val in self.iteritems():
            if isinstance(val,QIxml):
                val.replace(*args,**kwargs)
            elif isinstance(val,collection):
                for v in val:
                    if isinstance(v,QIxml):
                        v.replace(*args,**kwargs)
                    #end if
                #end for
            #end if
        #end for
    #end def replace


    def combine(self,other):
        #elemental combine only
        for name,element in other.iteritems():
            plural = isinstance(element,collection)
            single = isinstance(element,QIxml)
            if single or plural:
                elem = []
                single_name = None
                plural_name = None
                if single:
                    elem.append(element)
                    single_name = name
                    if name in plurals_inv:
                        plural_name = plurals_inv[name]
                    #end if
                else:
                    elem.extend(element.values())
                    plural_name = name
                    single_name = plurals[name]
                #end if
                if single_name in self:
                    elem.append(self[single_name])
                    del self[single_name]
                elif plural_name!=None and plural_name in self:
                    elem.append(self[plural_name])
                    del self[plural_name]
                #end if
                if len(elem)==1:
                    self[single_name]=elem[0]
                elif plural_name==None:
                    self.error('attempting to combine non-plural elements: '+single_name)
                else:
                    self[plural_name] = make_collection(elem)
                #end if
            #end if
        #end for
    #end def combine

                    
    def move(self,**elemdests):        
        names = elemdests.keys()
        hosts = self.get_host(names)
        dests = self.get(elemdests.values())
        if len(names)==1:
            hosts = [hosts]
            dests = [dests]
        #end if
        for i in range(len(names)):
            name = names[i]
            host = hosts[i]
            dest = dests[i]
            if host!=None and dest!=None and id(host)!=id(dest):
                if not name in host:
                    name = plurals_inv[name]
                #end if
                dest[name] = host[name]
                del host[name]
            #end if
        #end for
    #end def move



    def pluralize(self):
        make_plural = []
        for name,value in self.iteritems():
            if name in plurals_inv:
                make_plural.append(name)
            #end if
            if isinstance(value,QIxml):
                value.pluralize()
            elif isinstance(value,collection):
                for v in value:
                    if isinstance(v,QIxml):
                        v.pluralize()
                    #end if
                #end for
            #end if
        #end for
        for name in make_plural:
            value = self[name]
            del self[name]
            plural_name = plurals_inv[name]
            self[plural_name] = make_collection([value])
        #end for
    #end def pluralize


    def difference(self,other,root=True):
        if root:
            q1 = self.copy()
            q2 = other.copy()
        else:
            q1 = self
            q2 = other
        #end if
        if q1.__class__!=q2.__class__:
            different = True
            diff = None
            d1 = q1
            d2 = q2
        else:
            cls = q1.__class__
            s1 = set(q1.keys())
            s2 = set(q2.keys())
            shared  = s1 & s2
            unique1 = s1 - s2
            unique2 = s2 - s1
            different = len(unique1)>0 or len(unique2)>0
            diff = cls()
            d1 = cls()
            d2 = cls()
            d1.transfer_from(q1,unique1)
            d2.transfer_from(q2,unique2)
            for k in shared:
                value1 = q1[k]
                value2 = q2[k]
                is_coll1 = isinstance(value1,collection)
                is_coll2 = isinstance(value2,collection)
                is_qxml1 = isinstance(value1,QIxml)
                is_qxml2 = isinstance(value2,QIxml)
                if is_coll1!=is_coll2 or is_qxml1!=is_qxml2:
                    self.error('values for '+k+' are of inconsistent types\n  difference could not be taken')
                #end if
                if is_qxml1 and is_qxml2:
                    kdifferent,kdiff,kd1,kd2 = value1.difference(value2,root=False)
                elif is_coll1 and is_coll2:
                    ks1 = set(value1.keys())
                    ks2 = set(value2.keys())
                    kshared  = ks1 & ks2
                    kunique1 = ks1 - ks2
                    kunique2 = ks2 - ks1
                    kdifferent = len(kunique1)>0 or len(kunique2)>0
                    kd1 = collection()
                    kd2 = collection()
                    kd1.transfer_from(value1,kunique1)
                    kd2.transfer_from(value2,kunique2)
                    kdiff = collection()
                    for kk in kshared:
                        v1 = value1[kk]
                        v2 = value2[kk]
                        if isinstance(v1,QIxml) and isinstance(v2,QIxml):
                            kkdifferent,kkdiff,kkd1,kkd2 = v1.difference(v2,root=False)
                            kdifferent = kdifferent or kkdifferent
                            if kkdiff!=None:
                                kdiff[kk]=kkdiff
                            #end if
                            if kkd1!=None:
                                kd1[kk]=kkd1
                            #end if
                            if kkd2!=None:
                                kd2[kk]=kkd2
                            #end if
                        #end if
                    #end for
                else:
                    if isinstance(value1,ndarray):
                        a1 = value1.ravel()
                    else:
                        a1 = array([value1])
                    #end if
                    if isinstance(value2,ndarray):
                        a2 = value2.ravel()
                    else:
                        a2 = array([value2])
                    #end if
                    if len(a1)!=len(a2):
                        kdifferent = True
                    elif len(a1)==0:
                        kdifferent = False
                    elif (isinstance(a1[0],float) or isinstance(a2[0],float)) and not  (isinstance(a1[0],str)   or isinstance(a2[0],str)):
                        kdifferent = abs(a1-a2).max()/max(1e-99,abs(a1).max(),abs(a2).max()) > 1e-6
                    else:
                        kdifferent = not (a1==a2).all()
                    #end if
                    if kdifferent:
                        kdiff = (value1,value2)
                        kd1   = value1
                        kd2   = value2
                    else:
                        kdiff = None
                        kd1   = None
                        kd2   = None
                    #end if
                #end if
                different = different or kdifferent
                if kdiff!=None:
                    diff[k] = kdiff
                #end if
                if kd1!=None:
                    d1[k] = kd1
                #end if
                if kd2!=None:
                    d2[k] = kd2  
                #end if
            #end for
        #end if
        if root:
            if diff!=None:
                diff.remove_empty()
            #end if
            d1.remove_empty()
            d2.remove_empty()
        #end if 
        return different,diff,d1,d2
    #end def difference

    def remove_empty(self):
        names = list(self.keys())
        for name in names:
            value = self[name]
            if isinstance(value,QIxml):
                value.remove_empty()
                if len(value)==0:
                    del self[name]
                #end if
            elif isinstance(value,collection):
                ns = list(value.keys())
                for n in ns:
                    v = value[n]
                    if isinstance(v,QIxml):
                        v.remove_empty()
                        if len(v)==0:
                            del value[n]
                        #end if
                    #end if
                #end for
                if len(value)==0:
                    del self[name]
                #end if
            #end if
        #end for
    #end def remove_empty

    def get_host(self,names):
        return self.get(names,host=True)
    #end def get_host
#end class QIxml



class QIxmlFactory(Names):
    def __init__(self,name,types,typekey='',typeindex=-1,typekey2='',default=None):
        self.name = name
        self.types = types
        self.typekey = typekey
        self.typekey2 = typekey2
        self.typeindex = typeindex
        self.default = default
    #end def __init__

    def __call__(self,*args,**kwargs):
        #emulate QIxml.__init__
        #get the value of the typekey
        a  = args
        kw = kwargs
        found_type = False
        if len(args)>0:
            v = args[0]
            if isinstance(v,XMLelement):
                kw = v._attributes
            elif isinstance(v,section):
                a  = v.args
                kw = v.kwargs
            elif isinstance(v,tuple(self.types.values())):
                found_type = True
                type = v.__class__.__name__
            #end if
        #end if
        if not found_type:
            if self.typekey in kw.keys():
                type = kw[self.typekey]
            elif self.typekey2 in kw.keys():
                type = kw[self.typekey2]
            elif self.default!=None:
                type = self.default
            else:
                type = a[self.typeindex]
            #end if
        #end if
        type = self.condense_name(type)
        if type in self.types:
            return self.types[type](*args,**kwargs)
        else:
            msg = self.name+' factory is not aware of the following subtype:\n'
            msg+= '    '+type+'\n'
            self.error(msg,exit=False,trace=False)
        #end if
    #end def __call__

    def init_class(self):
        None # this is for compatibility with QIxml only (do not overwrite)
    #end def init_class
#end class QIxmlFactory



class Param(Names):        
    metadata = None

    def __call__(self,*args,**kwargs):
        if len(args)==0:
            self.error('no arguments provided, should have recieved one XMLelement')
        elif not isinstance(args[0],XMLelement):
            return args[0]
            #self.error('first argument is not an XMLelement')
        #end if
        return self.read(args[0])
    #end def __call__

    def read(self,xml):
        val = ''
        attr = set(xml._attributes.keys())
        other_attr = attr-set(['name'])
        if 'name' in attr and len(other_attr)>0:
            oa = obj()
            for a in other_attr:
                oa[a] = xml._attributes[a]
            #end for
            self.metadata[xml.name] = oa
        #end if
        if 'text' in xml:
            token = xml.text.split('\n',1)[0].split(None,1)[0]
            if is_int(token):
                val = loadtxt(StringIO(xml.text),int)
            elif is_float(token):
                val = loadtxt(StringIO(xml.text),float)
            else:
                val = array(xml.text.split())
            #end if
            if val.size==1:
                val = val.ravel()[0]
            #end if
        #end if
        return val
    #end def read


    def write(self,value,mode='attr',tag='parameter',name=None,pad='   ',write_type=None,normal_elem=False):
        c = ''
        attr_mode = mode=='attr'
        elem_mode  = mode=='elem'
        if not attr_mode and not elem_mode:
            self.error(mode+' is not a valid mode.  Options are attr,elem.')
        #end if
        if isinstance(value,list) or isinstance(value,tuple):
            value = array(value)
        #end if
        if attr_mode:
            if isinstance(value,ndarray):
                arr = value.ravel()
                for v in arr:
                    c+=str(v)+' '
                #end for
                c=c[:-1]
            else:
                c = str(value)
            #end if
        elif elem_mode:
            c+=pad
            is_array = isinstance(value,ndarray)
            is_single = not (is_array and value.size>1)
            if tag!=None:
                if is_single:
                    max_len = 20
                    rem_len = max(0,max_len-len(name))
                else:
                    rem_len = 0
                #end if
                other=''
                if name in self.metadata:
                    for a,v in self.metadata[name].iteritems():
                        other +=' '+self.expand_name(a)+'="'+str(v)+'"'
                    #end for
                #end if
                c+='<'+tag+' name="'+name+'"'+other+rem_len*' '+'>'
                pp = pad+'   '
            else:
                pp = pad
            #end if
            if is_array:
                if normal_elem:
                    c+='\n'
                #end if
                if tag!=None:
                    c+='\n'
                #end if
                ndim = len(value.shape)
                if ndim==1:
                    if tag!=None:
                        c+=pp
                    #end if
                    for v in value:
                        c+=str(v)+' '
                    #end for
                    c=c[:-1]+'\n'
                elif ndim==2:
                    nrows,ncols = value.shape
                    fmt=pp
                    if value.dtype == dtype(float):
                        vfmt = ':16.8f'
                    else:
                        vfmt = ''
                    #end if
                    for nc in range(ncols):
                        fmt+='{'+str(nc)+vfmt+'}  '
                    #end for
                    fmt = fmt[:-2]+'\n'
                    for nr in range(nrows):
                        c+=fmt.format(*value[nr])
                    #end for
                else:
                    self.error('only 1 and 2 dimensional arrays are supported for xml formatting.\n  Received '+ndim+' dimensional array.')
                #end if
            else:
                if write_type!=None:
                    val = write_type(value)
                else:
                    val = value
                #end if
                c += '    '+str(val)
            #end if
            if tag!=None:
                c+=pad+'</'+tag+'>\n'
            #end if
        #end if
        return c
    #end def write
            

    def init_class(self):
        None
    #end def init_class
#end class Param
param = Param()




class simulation(QIxml):
    elements   = ['project','random','qmcsystem','particleset','wavefunction','hamiltonian','init','traces','qmc','loop']
    write_types = obj(random=yesno)
#end class simulation


class project(QIxml):
    attributes = ['id','series']
    elements   = ['application','host','date','user']
#end class project

class application(QIxml):
    attributes = ['name','role','class','version']
#end class application


class random(QIxml):
    attributes = ['seed','parallel']
    write_types= obj(parallel=truefalse)
#end class random


class qmcsystem(QIxml):
    elements = ['simulationcell','particleset','wavefunction','hamiltonian']
#end class qmcsystem



class simulationcell(QIxml):
    attributes = ['name']
    parameters = ['lattice','reciprocal','bconds','lr_dim_cutoff']
#end class simulationcell

class particleset(QIxml):
    attributes = ['name','size','random','random_source','randomsrc','charge']
    elements   = ['group']
    attribs    = ['ionid','position']
    write_types= obj(random=yesno)
    identifier = 'name'
#end class particleset

class group(QIxml):
    attributes = ['name','size']
    parameters = ['charge','valence','atomicnumber','mass']
    attribs    = ['position']
    identifier = 'name'
#end class group



class wavefunction(QIxml):
    attributes = ['name','target']
    elements   = ['determinantset','jastrow']
    identifier = 'name'
#end class wavefunction

class determinantset(QIxml):
    attributes = ['type','href','sort','tilematrix','twistnum','source','version','meshfactor','gpu','transform','precision','truncate','lr_dim_cutoff']
    elements   = ['basisset','slaterdeterminant']
    h5tags     = ['twistindex','twistangle']
    write_types = obj(gpu=yesno,sort=onezero,transform=yesno,truncate=yesno)
#end class determinantset

class basisset(QIxml):
    attributes = ['ecut']
    elements   = ['grid','atomicbasisset']
#end class basisset

class grid(QIxml):
    attributes = ['dir','npts','closed']
    identifier = 'dir'
#end class grid

class atomicbasisset(QIxml):
    attributes = ['type','elementtype','expandylm','href']
    elements   = ['basisgroup']
    write_types= obj(expandylm=yesno)
#end class atomicbasisset

class basisgroup(QIxml):
    attributes = ['rid','ds','n','l','m','zeta']
#end class basisgroup

class slaterdeterminant(QIxml):
    elements   = ['determinant']
#end class slaterdeterminant

class determinant(QIxml):
    attributes = ['id','size','ref','spin']
    elements   = ['occupation']
    identifier = 'id'
#end class determinant

class occupation(QIxml):
    attributes = ['mode','spindataset']
#end class occupation

class jastrow1(QIxml):
    tag = 'jastrow'
    attributes = ['type','name','function','source','print']
    elements   = ['correlation']
    identifier = 'name'
    write_types = obj(print_=yesno)
#end class jastrow1

class jastrow2(QIxml):
    tag = 'jastrow'
    attributes = ['type','name','function','print','spin']
    elements   = ['correlation']
    parameters = ['b']
    identifier = 'name'
    write_types = obj(print_=yesno)
#end class jastrow2

class jastrow3(QIxml):
    tag = 'jastrow'
    attributes = ['type','name','function','print','source']
    elements   = ['correlation']
    identifier = 'name'
    write_types = obj(print_=yesno)
#end class jastrow3

class correlation(QIxml):
    attributes = ['elementtype','speciesa','speciesb','size','ispecies','especies',
                  'especies1','especies2','isize','esize','rcut','cusp']
    elements   = ['coefficients','var']
    identifier = 'speciesa','speciesb','elementtype','especies1','especies2','ispecies'
#end class correlation

class var(QIxml):
    attributes = ['id','name']
    text       = 'value'
    identifier = 'id'
#end class var

class coefficients(QIxml):
    attributes = ['id','type','optimize']
    text       = 'coeff'
    write_types= obj(optimize=yesno)
#end class coefficients

jastrow = QIxmlFactory(
    name = 'jastrow',
    types   = dict(one_body=jastrow1,two_body=jastrow2,jastrow1=jastrow1,jastrow2=jastrow2,eei=jastrow3,jastrow3=jastrow3),
    typekey = 'type'
    )



class hamiltonian(QIxml):
    attributes = ['name','type','target'] 
    elements   = ['pairpot','constant','estimator']
    identifier = 'name'
#end class hamiltonian

class coulomb(QIxml):
    tag = 'pairpot'
    attributes  = ['type','name','source','target','physical']
    write_types = obj(physical=yesno)
    identifier  = 'name'
#end class coulomb

class constant(QIxml):
    attributes = ['type','name','source','target']
#end class constant

class pseudopotential(QIxml):
    tag = 'pairpot'
    attributes = ['type','name','source','wavefunction','format']
    elements   = ['pseudo']
    identifier = 'name'
#end class pseudopotential

class pseudo(QIxml):
    attributes = ['elementtype','href','format']
    identifier = 'elementtype'
#end class pseudo

class mpc(QIxml):
    tag='pairpot'
    attributes=['type','name','source','target','ecut','physical']
    write_types = obj(physical=yesno)    
    identifier='name'
#end class mpc

pairpot = QIxmlFactory(
    name  = 'pairpot',
    types = dict(coulomb=coulomb,pseudo=pseudopotential,
                 pseudopotential=pseudopotential,mpc=mpc),
    typekey = 'type'
    )



class localenergy(QIxml):
    tag = 'estimator'
    attributes = ['name','hdf5']
    write_types= obj(hdf5=yesno)
    identifier = 'name'
#end class localenergy

class energydensity(QIxml):
    tag = 'estimator'
    attributes = ['type','name','dynamic','static']
    elements   = ['spacegrid']
    identifier  = 'name'
#end class energydensity

class spacegrid(QIxml):
    attributes = ['coord','min_part','max_part']
    elements   = ['origin','axis']
#end class spacegrid

class origin(QIxml):
    attributes = ['p1','p2']
#end class origin

class axis(QIxml):
    attributes = ['p1','p2','scale','label','grid']
    identifier = 'label'
#end class axis

class chiesa(QIxml):
    tag = 'estimator'
    attributes = ['name','type','source','psi']
    identifier = 'name'
#end class chiesa

class density(QIxml):
    tag = 'estimator'
    attributes = ['name','type','delta']
    identifier = 'type'
#end class density

estimator = QIxmlFactory(
    name = 'estimator',
    types   = dict(localenergy=localenergy,energydensity=energydensity,chiesa=chiesa,density=density),
    typekey = 'type',
    typekey2= 'name'
    )


class init(QIxml):
    attributes = ['source','target']
#end class


class scalar_traces(QIxml):
    attributes  = ['defaults']
    text        = 'quantities'
    write_types = obj(defaults=yesno)
#end class scalar_traces

class particle_traces(QIxml):
    attributes  = ['defaults']
    text        = 'quantities'
    write_types = obj(defaults=yesno)
#end class particle_traces

class traces(QIxml):
    attributes = ['write','format','verbose','scalar','particle',
                  'scalar_defaults','particle_defaults']
    elements = ['scalar_traces','particle_traces']
    write_types = obj(write=yesno,verbose=yesno,scalar=yesno,particle=yesno,
                      scalar_defaults=yesno,particle_defaults=yesno)
#end class


class loop(QIxml):
    attributes = ['max']
    elements = ['qmc']
    def unroll(self):
        calculations=[]
        calcs = []
        if 'qmc' in self:
            calcs = [self.qmc]
        elif 'calculations' in self:
            calcs = self.calculations
        #end if
        for n in range(self.max):
            for i in range(len(calcs)):
                calculations.append(calcs[i].copy())
            #end for
        #end for
        return make_collection(calculations)
    #end def unroll
#end class loop

class linear(QIxml):
    tag = 'qmc'
    attributes = ['method','move','checkpoint','gpu','trace']
    elements   = ['estimator']
    parameters = ['blocks','warmupsteps','stepsbetweensamples','timestep',
                  'samples','minwalkers','maxweight','usedrift','minmethod',
                  'beta','exp0','bigchange','alloweddifference','stepsize',
                  'stabilizerscale','nstabilizers','max_its','cgsteps','eigcg',
                  'walkers','nonlocalpp','usebuffer','gevmethod','steps','substeps']
    costs      = ['energy','unreweightedvariance','reweightedvariance']
    write_types = obj(gpu=yesno,usedrift=yesno,nonlocalpp=yesno,usebuffer=yesno)
#end class linear

class cslinear(QIxml):
    tag = 'qmc'
    attributes = ['method','move','checkpoint','gpu','trace']
    elements   = ['estimator']
    parameters = ['blocks','warmupsteps','stepsbetweensamples','steps','samples','timestep','usedrift',
                  'minmethod','gevmethod','exp0','nstabilizers','stabilizerscale',
                  'stepsize','alloweddifference','beta','bigchange','minwalkers',
                  'usebuffer','maxweight','nonlocalpp','max_its','walkers','substeps']
    costs      = ['energy','unreweightedvariance','reweightedvariance']
    write_types = obj(gpu=yesno,usedrift=yesno,nonlocalpp=yesno,usebuffer=yesno)
#end class cslinear

class vmc(QIxml):
    tag = 'qmc'
    attributes = ['method','multiple','warp','move','gpu','checkpoint','trace']
    elements   = ['estimator']
    parameters = ['walkers','blocks','steps','substeps','timestep','usedrift','warmupsteps','samples','nonlocalpp','stepsbetweensamples','samplesperthread']
    write_types = obj(gpu=yesno,usedrift=yesno,nonlocalpp=yesno)
#end class vmc

class dmc(QIxml):
    tag = 'qmc'
    attributes = ['method','move','gpu','multiple','warp','checkpoint','trace']
    elements   = ['estimator']
    parameters = ['walkers','blocks','steps','timestep','nonlocalmove','nonlocalmoves','warmupsteps','pop_control']
    write_types = obj(gpu=yesno,nonlocalmoves=yesno)
#end class dmc

qmc = QIxmlFactory(
    name = 'qmc',
    types   = dict(linear=linear,cslinear=cslinear,vmc=vmc,dmc=dmc,loop=loop),
    typekey = 'method',
    default = 'loop'
    )


class gen(QIxml):
    attributes = []
    elements   = []
#end class gen


classes = [   #standard classes
    simulation,project,application,random,qmcsystem,simulationcell,particleset,
    group,hamiltonian,constant,pseudopotential,coulomb,pseudo,mpc,chiesa,density,
    localenergy,energydensity,spacegrid,origin,axis,wavefunction,
    determinantset,slaterdeterminant,basisset,grid,determinant,occupation,
    jastrow1,jastrow2,jastrow3,
    correlation,coefficients,loop,linear,cslinear,vmc,dmc,
    atomicbasisset,basisgroup,init,var,traces,scalar_traces,particle_traces
    ]
types = dict( #simple types and factories
    host      = param,
    date      = param,
    user      = param,
    pairpot   = pairpot,
    estimator = estimator,
    jastrow   = jastrow,
    qmc       = qmc
    )
plurals = obj(
    particlesets  = 'particleset',
    groups        = 'group',    
    hamiltonians  = 'hamiltonian',
    pairpots      = 'pairpot',
    pseudos       = 'pseudo',
    estimators    = 'estimator',
    spacegrids    = 'spacegrid',
    axes          = 'axis',
    wavefunctions = 'wavefunction',
    grids         = 'grid',
    determinants  = 'determinant',
    correlations  = 'correlation',
    jastrows      = 'jastrow',
    basisgroups   = 'basisgroup',
    calculations  = 'qmc',
    vars          = 'var'
    )
plurals_inv = plurals.inverse()
plural_names = set(plurals.keys())
single_names = set(plurals.values())
Names.set_expanded_names(
    elementtype     = 'elementType',
    energydensity   = 'EnergyDensity',
    gevmethod       = 'GEVMethod',
    localenergy     = 'LocalEnergy',
    lr_dim_cutoff   = 'LR_dim_cutoff',
    minmethod       = 'MinMethod',
    one_body        = 'One-Body',
    speciesa        = 'speciesA',
    speciesb        = 'speciesB',
    substeps        = 'subSteps',
    two_body        = 'Two-Body',
    usedrift        = 'useDrift',
    maxweight       = 'maxWeight',
    warmupsteps     = 'warmupSteps',
    twistindex      = 'twistIndex',
    twistangle      = 'twistAngle',
    usebuffer       = 'useBuffer',
    mpc             = 'MPC',
    kecorr          = 'KEcorr',
    ionion          = 'IonIon',
    elecelec        = 'ElecElec',
    pseudopot       = 'PseudoPot',
    posarray        = 'posArray',
    array           = 'Array',
    atomicbasisset  = 'atomicBasisSet',
    basisgroup      = 'basisGroup',
    expandylm       = 'expandYlm',
    mo              = 'MO',
    numerical       = 'Numerical'
    )
for c in classes:
    c.init_class()
    types[c.__name__] = c
#end for


#set default values
simulation.defaults.set(
    project      = project,
    qmcsystem    = qmcsystem,
    calculations = lambda:list()
    )
project.defaults.set(
    series=0,
    application = application
    )
application.defaults.set(
    name='qmcapp',role='molecu',class_='serial',version='1.0'
    )
simulationcell.defaults.set(
    bconds = 'p p p',lr_dim_cutoff=15
    )
wavefunction.defaults.set(
    name='psi0',target='e'
    )
determinantset.defaults.set(
    type='einspline',tilematrix=lambda:eye(3,dtype=int),twistnum=0,meshfactor=1.,gpu=False,precision='double'
    )
occupation.defaults.set(
    mode='ground',spindataset=0
    )
jastrow1.defaults.set(
    name='J1',type='one-body',function='bspline',print_=True,source='ion0',
    correlation=correlation
    )
jastrow2.defaults.set(
    name='J2',type='two-body',function='bspline',print_=True,
    correlation=correlation
    )
jastrow3.defaults.set(
    name='J3',type='eeI',function='polynomial',print_=True,source='ion0',
    correlation=correlation
    )
correlation.defaults.set(
    coefficients=coefficients
    )
coefficients.defaults.set(
    type='Array'
    )
hamiltonian.defaults.set(
    name='h0',type='generic',target='e',
    constant = constant,
    pairpots = classcollection(coulomb,pseudopotential,mpc),
    estimators = classcollection(chiesa),
    )
coulomb.defaults.set(
    name='ElecElec',type='coulomb',source='e',target='e'
    )
constant.defaults.set(
    name='IonIon',type='coulomb',source='ion0',target='ion0'
    )
pseudopotential.defaults.set(
    name='PseudoPot',type='pseudo',source='ion0',wavefunction='psi0',format='xml'
    )
mpc.defaults.set(
    name='MPC',type='MPC',ecut=60.0,source='e',target='e',physical=False
    )
localenergy.defaults.set(
    name='LocalEnergy',hdf5=True
    )
chiesa.defaults.set(
    name='KEcorr',type='chiesa',source='e',psi='psi0'
    )
energydensity.defaults.set(
    type='EnergyDensity',name='EDvoronoi',dynamic='e',static='ion0',
    spacegrid = spacegrid
    )
spacegrid.defaults.set(
    coord='voronoi'
    )




linear.defaults.set(
     method = 'linear',move='pbyp',checkpoint=-1,gpu=False,
     estimators = classcollection(localenergy)
#  #jtk
#    method='linear',move='pbyp',checkpoint=-1,gpu=True,
#    energy=0, reweightedvariance=0, unreweightedvariance=0,
#    warmupsteps       = 20,
#    usedrift          = True,
#    timestep          = .5,
#    minmethod         ='rescale',
#    stepsize          = .5,
#    beta              = 0.05,
#    alloweddifference = 1e-8,
#    bigchange         = 1.1,
#    cgsteps           = 3,
#    eigcg             = 1,
#    exp0              = -6,
#    maxweight         = 1e9,
#    minwalkers        = .5,
#    nstabilizers      = 10,
#    stabilizerscale   = .5,
#    usebuffer         = True,
    )
cslinear.defaults.set(
    method='cslinear', move='pbyp', checkpoint=-1, gpu=False,
    estimators = classcollection(localenergy)
  #jtk
    #method='cslinear',move='pbyp',checkpoint=-1,gpu=True,
    #energy=0,reweightedvariance=0,unreweightedvariance=1.,
    #warmupsteps=5,steps=2,usedrift=True,timestep=.5,
    #minmethod='quartic',gevmethod='mixed',exp0=-15,
    #nstabilizers=5,stabilizerscale=3,stepsize=.35,
    #alloweddifference=1e-5,beta=.05,bigchange=5.,
    #estimators=classcollection(localenergy)
  #lschulen
    #method='cslinear', move='pbyp', checkpoint=-1, gpu=True,
    #energy=0, reweightedvariance=0, unreweightedvariance=0,
    #warmupsteps       = 20,
    ##steps             = 5,
    #usedrift          = True,
    #timestep          = .8,
    #nonlocalpp        = False,
    #minmethod         = 'rescale',
    #stepsize          = .4,
    #beta              = .05,
    #gevmethod         = 'mixed',
    #alloweddifference = 1e-4,
    #bigchange         = 9.,
    #exp0              = -16,
    #max_its           = 1,
    #maxweight         = 1e9,
    #minwalkers        = .5,
    #nstabilizers      = 3,
    #stabilizerscale   = 1,
    #usebuffer         = False,
    #estimators = classcollection(localenergy)
  #jmm
    #method='cslinear', move='pbyp', checkpoint=-1, gpu=True,
    #energy=0, reweightedvariance=0, unreweightedvariance=0,
    #warmupsteps       = 20,
    #usedrift          = True,
    #timestep          = .5,
    #nonlocalpp        = True,
    #minmethod         = 'quartic',
    #stepsize          = .4,
    #beta              = 0.0,
    #gevmethod         = 'mixed',
    #alloweddifference = 1.0e-4,
    #bigchange         = 9.0,
    #exp0              = -16,
    #max_its           = 1,
    #maxweight         = 1e9,
    #minwalkers        = 0.5,
    #nstabilizers      = 3,
    #stabilizerscale   = 1.0,
    #usebuffer         = True,
    #estimators = classcollection(localenergy)
    )
vmc.defaults.set(
    method='vmc',move='pbyp',gpu=False,
    walkers     = 1,
    warmupsteps = 50,
    substeps    = 3,
    usedrift    = True,
    timestep    = .5,
    estimators = classcollection(localenergy)
    )
dmc.defaults.set(
    method='dmc',move='pbyp',gpu=False,
    warmupsteps   = 20,
    timestep      = .01,
    nonlocalmoves = True,
    estimators = classcollection(localenergy)
    )



opt_defaults = obj(
    linear = obj(
        jmm = linear(
            warmupsteps         = 100,
            timestep            = 0.5,
            stepsbetweensamples = 10,
            minwalkers          = 0.0,
            bigchange           = 15.0,
            alloweddifference   = 1.e-4
            )
        ),
    cslinear = obj(
        )
    )





class QmcpackInput(SimulationInput,Names):
    
    opt_methods = set(['opt','linear','cslinear'])

    simulation_type = simulation

    default_metadata = meta(
        lattice    = dict(units='bohr'),
        reciprocal = dict(units='2pi/bohr'),
        ionid      = dict(datatype='stringArray'),
        position   = dict(datatype='posArray', condition=0)
        )

    def __init__(self,arg0=None,arg1=None):
        Param.metadata = None
        filepath = None
        metadata = None
        element  = None
        if arg0==None and arg1==None:
            None
        elif isinstance(arg0,str) and arg1==None:
            filepath = arg0
        elif isinstance(arg0,QIxml) and arg1==None:
            element = arg0
        elif isinstance(arg0,meta) and isinstance(arg1,QIxml):
            metadata = arg0
            element  = arg1
        else:
            self.error('input arguments of types '+arg0.__class__.__name__+' and '+arg0.__class__.__name__+' cannot be used to initialize QmcpackInput')
        #end if
        if metadata!=None:
            self._metadata = metadata
        else:
            self._metadata = meta()
        #end if
        if filepath!=None:
            self.read(filepath)
        elif element!=None:
            #simulation = arg0
            #self.simulation = self.simulation_type(simulation)
            elem_class = element.__class__
            if 'identifier' in elem_class.__dict__:
                name = elem_class.identifier
            else:
                name = elem_class.__name__
            #end if
            self[name] = elem_class(element)
        #end if
        Param.metadata = None
    #end def __init__

    def get_base(self):
        elem_names = list(self.keys())
        elem_names.remove('_metadata')
        if len(elem_names)>1:
            self.error('qmcpack input cannot have more than one base element\n  You have provided '+str(len(elem_names))+': '+str(elem_names))
        #end if
        return self[elem_names[0]]
    #end def get_base

    def get_basename(self):
        elem_names = list(self.keys())
        elem_names.remove('_metadata')
        if len(elem_names)>1:
            self.error('qmcpack input cannot have more than one base element\n  You have provided '+str(len(elem_names))+': '+str(elem_names))
        #end if
        return elem_names[0]
    #end def get_basename

    def read(self,filepath=None,xml=None):
        if xml!=None or os.path.exists(filepath):
            element_joins=['qmcsystem']
            element_aliases=dict(loop='qmc')
            xml = XMLreader(filepath,element_joins,element_aliases,warn=False,xml=xml).obj
            xml.condense()
            self._metadata = meta() #store parameter/attrib attribute metadata
            Param.metadata = self._metadata
            if 'simulation' in xml:
                self.simulation = simulation(xml.simulation)
            else:
                #try to determine the type
                elements = []
                keys = []
                error = False
                for key,value in xml.iteritems():
                    if isinstance(key,str) and key[0]!='_':
                        if key in types:
                            elements.append(types[key](value))
                            keys.append(key)
                        else:
                            self.error('element '+key+' is not a recognized type',exit=False)
                            error = True
                        #end if
                    #end if
                #end for
                if error:
                    self.error('cannot read input xml file')
                #end if
                if len(elements)==0:
                    self.error('no valid elements were found for input xml file')
                #end if
                for i in range(len(elements)):
                    elem = elements[i]
                    key  = keys[i]
                    if isinstance(elem,QIxml):
                        if 'identifier' in elem.__class__.__dict__:
                            name = elem.identifier
                        else:
                            name = elem.tag
                        #end if
                    else:
                        name = key
                    #end if
                    self[name] = elem
                #end for
            #end if
            Param.metadata = None
        else:
            self.error('the filepath you provided does not exist.\n  Input filepath: '+filepath)
        #end if
        return self
    #end def read


    def write_contents(self):
        c = ''
        header = '''<?xml version="1.0"?>
'''
        c+= header
        if len(self._metadata)==0:
            Param.metadata = self.default_metadata
        else:
            Param.metadata = self._metadata
        #end if
        base = self.get_base()
        c+=base.write(first=True)
        Param.metadata = None
        return c
    #end def write_contents


    def unroll_calculations(self,modify=True):
        qmc = []
        sim = self.simulation
        if 'calculations' in sim:
            calcs = sim.calculations
        elif 'qmc' in sim:
            calcs = [sim.qmc]
        else:
            calcs = []
        #end if
        for i in range(len(calcs)):
            c = calcs[i]
            if isinstance(c,loop):
                qmc.extend(c.unroll())
            else:
                qmc.append(c)
            #end if
        #end for
        qmc = make_collection(qmc)
        if modify:
            self.simulation.calculations = qmc
        #end if
        return qmc
    #end def unroll_calculations

    def get(self,*names):
        base = self.get_base()
        return base.get(names)
    #end def get

    def remove(self,*names):
        base = self.get_base()
        base.remove(*names)
    #end def remove
    
    def assign(self,**kwargs):
        base = self.get_base()
        base.assign(**kwargs)
    #end def assign
    
    def replace(self,*args,**kwargs):# input is list of keyword=(oldval,newval)
        base = self.get_base()
        base.replace(*args,**kwargs)
    #end def replace

    def move(self,**elemdests):
        base = self.get_base()
        base.move(**elemdests)
    #end def move
            

    def get_host(self,names):
        base = self.get_base()
        return base.get_host(names)
    #end if

    def incorporate_defaults(self,elements=False,overwrite=False,propagate=False):
        base = self.get_base()
        base.incorporate_defaults(elements,overwrite,propagate)
    #end def incorporate_defaults

    def pluralize(self):
        base = self.get_base()
        base.pluralize()
    #end def pluralize

    def standard_placements(self):
        self.move(particleset='qmcsystem',wavefunction='qmcsystem',hamiltonian='qmcsystem')
    #end def standard_placements

    def difference(self,other):
        s1 = self.copy()
        s2 = other.copy()
        b1 = s1.get_basename()
        b2 = s2.get_basename()
        q1 = s1[b1]
        q2 = s2[b2]
        if b1!=b2:
            different = True
            d1 = q1
            d2 = q2
            diff = None
        else:
            s1.standard_placements()
            s2.standard_placements()
            s1.pluralize()
            s2.pluralize()
            different,diff,d1,d2 = q1.difference(q2,root=False)
        #end if
        if diff!=None:
            diff.remove_empty()
        #end if
        d1.remove_empty()
        d2.remove_empty()
        return different,diff,d1,d2
    #end def difference

    def remove_empty(self):
        base = self.get_base()
        base.remove_empty()
    #end def remove_empty

    def read_xml(self,filepath=None,xml=None):
        if os.path.exists(filepath):
            element_joins=['qmcsystem']
            element_aliases=dict(loop='qmc')
            if xml is None:
                xml = XMLreader(filepath,element_joins,element_aliases,warn=False).obj
            else:
                xml = XMLreader(None,element_joins,element_aliases,warn=False,xml=xml).obj
            #end if
            xml.condense()
        else:
            self.error('the filepath you provided does not exist.\n  Input filepath: '+filepath)
        #end if
        return xml
    #end def read_xml

    def include_xml(self,xmlfile,replace=True,exists=True):
        xml = self.read_xml(xmlfile)
        Param.metadata = self._metadata
        for name,exml in xml.iteritems():
            if not name.startswith('_'):
                qxml = types[name](exml)
                qname = qxml.tag
                host = self.get_host(qname)
                if host==None and exists:
                    self.error('host xml section for '+qname+' not found','QmcpackInput')
                #end if
                if qname in host:
                    section_name = qname
                elif qname in plurals_inv and plurals_inv[qname] in host:
                    section_name = plurals_inv[qname]
                else:
                    section_name = None
                #end if
                if replace:
                    if section_name!=None:
                        del host[section_name]
                    #end if
                    host[qname] = qxml
                else:
                    if section_name==None:
                        host[qname] = qxml
                    else:
                        section = host[section_name]
                        if isinstance(section,collection):
                            section[qxml.identifier] = qxml
                        elif section_name in plurals_inv:
                            coll = collection()
                            coll[section.identifier] = section
                            coll[qxml.identifier]    = qxml
                            del host[section_name]
                            host[plurals_inv[section_name]] = coll
                        else:
                            section.combine(qxml)
                        #end if
                    #end if
                #end if
            #end if
        #end for
        Param.metadata = None
    #end def include_xml

    def get_output_info(self,*requests):
        project = self.simulation.project
        prefix = project.id
        series = project.series
        qmc_ur = self.unroll_calculations(modify=False)

        qmc = []
        calctypes = set()
        outfiles = []
        n=0
        for qo in qmc_ur:
            q = obj()
            q.prefix = prefix
            q.series = series+n
            n+=1
            method = qo.method
            if method in self.opt_methods:
                q.type = 'opt'
            else:
                q.type = method
            #end if
            calctypes.add(q.type)
            q.method = method
            fprefix = prefix+'.s'+str(q.series).zfill(3)+'.'
            files = obj()
            files.scalar = fprefix+'scalar.dat'
            files.stat   = fprefix+'stat.h5'
            # apparently this one is no longer generated by default as of r5756
            #files.config = fprefix+'storeConfig.h5' 
            if q.type=='opt':
                files.opt = fprefix+'opt.xml'
            elif q.type=='dmc':
                files.dmc = fprefix+'dmc.dat'
            #end if
            outfiles.extend(files.values())
            q.files = files
            qmc.append(q)
        #end for
        res = dict(qmc=qmc,calctypes=calctypes,outfiles=outfiles)

        values = []
        for req in requests:
            if req in res:
                values.append(res[req])
            else:
                self.error(req+' is not a valid output info request')
            #end if
        #end for
        if len(values)==1:
            return values[0]
        else:
            return values
        #end if
    #end def get_output_info


    def generate_jastrows(self,size=None,j1func='bspline',j1size=8,j2func='bspline',j2size=8):
        if size!=None:
            j1size = size
            j2size = size
        #end if

        #self.remove('jastrow')
        lattice,particlesets,wavefunction = self.get('lattice','particleset','wavefunction')
        no_lattice = lattice==None
        no_particleset = particlesets==None
        no_wavefunction = wavefunction==None
        if no_lattice:
            self.error('a simulationcell lattice must be present to generate jastrows',exit=False)
        #end if
        if no_particleset:
            self.error('a particleset must be present to generate jastrows',exit=False)
        #end if
        if no_wavefunction:
            self.error('a wavefunction must be present to generate jastrows',exit=False)
        #end if
        if no_lattice or no_particleset or no_wavefunction:
            self.error('jastrows cannot be generated')
        #end if
        if isinstance(particlesets,QIxml):
            particlesets = make_collection([particlesets])
        #end if
        if not 'e' in particlesets:
            self.error('electron particleset (e) not found\n particlesets: '+str(particlesets.keys()))
        #end if


        jastrows = collection()

        cell = Structure(lattice)
        volume = cell.volume()
        rcut   = cell.rmin()

        #use the rpa jastrow for electrons (modeled after Luke's tool)
        size = j2size
        e = particlesets.e
        nelectrons = 0
        for g in e.groups:
            nelectrons += g.size
        #end for
        density = nelectrons/volume
        wp = sqrt(4*pi*density)
        dr = rcut/size
        r = .02 + dr*arange(size)
        uuc = .5/(wp*r)*(1.-exp(-r*sqrt(wp/2)))*exp(-(2*r/rcut)**2)
        udc = .5/(wp*r)*(1.-exp(-r*sqrt(wp)))*exp(-(2*r/rcut)**2)
        jastrows.J2 = jastrow2(
            name = 'J2',type='Two-Body',function=j2func,print_='yes',
            correlations = collection(
                uu = correlation(speciesA='u',speciesB='u',size=size,
                                 coefficients=section(id='uu',type='Array',coeff=uuc)),
                ud = correlation(speciesA='u',speciesB='d',size=size,
                                 coefficients=section(id='ud',type='Array',coeff=udc))
                )
            )

        #generate electron-ion jastrows, if ions present
        ions = []
        for name in particlesets.keys():
            if name=='i' or name.startswith('ion'):
                ions.append(name)
            #end if
        #end for
        if len(ions)>0:
            size = j1size
            j1 = []
            for ion in ions:
                i = particlesets[ion]
                if 'group' in i:
                    groups = [i.group]
                else:
                    groups = i.groups
                #end if
                corr = []
                for g in groups:
                    elem = g.name
                    c=correlation(
                        elementtype=elem,
                        cusp=0.,
                        size=size,
                        coefficients=section(
                            id='e'+elem,
                            type='Array',
                            coeff=size*[0]
                            )
                        )
                    corr.append(c)
                #end for
                j=jastrow1(
                    name='J1_'+ion,
                    type='One-Body',
                    function=j1func,
                    source=ion,
                    print_='yes',
                    correlations = corr
                    )
                j1.append(j)
            #end for
            if len(j1)==1:
                j1[0].name='J1'
            #end if
            for j in j1:
                jastrows[j.name]=j
            #end for
        #end if

        if 'J2' in wavefunction.jastrows:
            J2 = wavefunction.jastrows.J2
            if 'function' in J2 and J2.function.lower()=='bspline':
                c = wavefunction.jastrows.J2.correlations
                ctot = abs(array(c.uu.coefficients.coeff)).sum() + abs(array(c.ud.coefficients.coeff)).sum()
                if ctot < 1e-3:
                    wavefunction.jastrows.J2 = jastrows.J2
                #end if
            #end if
        #end if

        #only add the jastrows if ones of the same type 
        # (one-body,two-body,etc) are not already present
        for jastrow in jastrows:
            jtype = jastrow.type.lower().replace('-','_')
            has_jtype = False
            for wjastrow in wavefunction.jastrows:
                wjtype = wjastrow.type.lower().replace('-','_')
                has_jtype = has_jtype or wjtype==jtype
            #end for
            if not has_jtype:
                wavefunction.jastrows[jastrow.name] = jastrow
            #end if
        #end for
    #end def generate_jastrows


    def incorporate_system(self,system):
        system.check_folded_system()
        system.change_units('B')
        system.structure.group_atoms()
        particles  = system.particles
        structure  = system.structure
        net_charge = system.net_charge
        net_spin   = system.net_spin

        qs,sc,ham,ps = self.get('qmcsystem','simulationcell','hamiltonian','particleset')

        old_eps_name = None
        old_ips_name = None
        if ps!=None:
            if isinstance(ps,particleset):
                ps = make_collection([ps])
            #end if
            for pname,pset in ps.iteritems():
                g0name = pset.groups.keys()[0]
                g0 = pset.groups[g0name]
                if abs(-1-g0.charge)<1e-2:
                    old_eps_name = pname
                elif 'ionid' in pset:
                    old_ips_name = pname
                #end if
            #end for
        #end if
        del ps
        self.remove('particleset')
        if qs==None:
            qs = qmcsystem()
            qs.incorporate_defaults(elements=False,propagate=False)
            self.simulation.qmcsystem = qs
        #end if
        if sc==None:
            sc = simulationcell()
            sc.incorporate_defaults(elements=False,propagate=False)
            qs.simulationcell = sc
        #end if
        if ham==None:
            ham = hamiltonian()
            ham.incorporate_defaults(elements=False,propagate=False)
            qs.hamiltonian = ham
        elif isinstance(ham,collection):
            if 'h0' in ham:
                ham = ham.h0
            elif len(ham)==1:
                ham = ham.list()[0]
            else:
                self.error('cannot find hamiltonian for system incorporation')
            #end if
        #end if

        elem = structure.elem
        pos  = structure.pos

        if len(structure.axes)>0: #exclude systems with open boundaries
            #setting the 'lattice' (cell axes) requires some delicate care
            #  qmcpack will fail if this is even 1e-10 off of what is in 
            #  the wavefunction hdf5 file from pwscf
            if structure.folded_structure!=None:
                fs = structure.folded_structure
                axes = array(pwscf_array_string(fs.axes).split(),dtype=float)
                axes.shape = fs.axes.shape
                axes = dot(axes,structure.tmatrix)
                if abs(axes-structure.axes).sum()>1e-5:
                    self.error('supercell axes do not match tiled version of folded cell axes\n  you may have changed one set of axes (super/folded) and not the other\n  folded cell axes:\n'+str(fs.axes)+'\n  supercell axes:\n'+str(structure.axes)+'\n  folded axes tiled:\n'+str(axes))
                #end if
            else:
                axes = array(pwscf_array_string(structure.axes).split(),dtype=float)
                axes.shape = structure.axes.shape
            #end if
            structure.adjust_axes(axes)

            sc.lattice = axes
        #end if    

        elns = particles.get_electrons()
        ions = particles.get_ions()
        eup  = elns.up_electron
        edn  = elns.down_electron

        particlesets = []
        eps = particleset(
            name='e',random=True,
            groups = [
                group(name='u',charge=-1,mass=eup.mass,size=eup.count),
                group(name='d',charge=-1,mass=edn.mass,size=edn.count)
                ]
            )
        particlesets.append(eps)
        if len(ions)>0:
            eps.randomsrc = 'ion0'
            ips = particleset(
                name='ion0',
                )
            groups = []
            ham.pluralize()
            pseudos = ham.get('pseudo')
            if pseudos==None:
                pp = ham.get('PseudoPot')
                if pp!=None:
                    pseudos = collection()
                    pp.pseudos = pseudos
                #end if
            #end if
            for ion in ions:
                gpos = pos[elem==ion.name]
                g = group(
                    name         = ion.name,
                    charge       = ion.charge,
                    valence      = ion.charge,
                    atomicnumber = ion.protons,
                    mass         = ion.mass,
                    position     = gpos,
                    size         = len(gpos)
                    )
                groups.append(g)
                if pseudos!=None and not ion.name in pseudos:
                    pseudos[ion.name] = pseudo(elementtype=ion.name,href='MISSING.xml')
                #end if
            #end for
            ips.groups = make_collection(groups)
            particlesets.append(ips)
        #end if
        qs.particlesets = make_collection(particlesets)

        if old_eps_name!=None:
            self.replace(old_eps_name,'e')
        #end if
        if old_ips_name!=None and len(ions)>0:
            self.replace(old_ips_name,'ion0')
        #end if
            
        udet,ddet = self.get('updet','downdet')

        if udet!=None:
            udet.size = elns.up_electron.count
        #end if
        if ddet!=None:
            ddet.size = elns.down_electron.count
        #end if

        if abs(net_spin) > 1e-1:
            if ddet!=None:
                ddet.occupation.spindataset = 1 #jtk mark check
            #end if
        #end if
    #end def incorporate_system
        

    def return_system(self):
        input = self.copy()
        input.pluralize()
        axes,ps = input.get('lattice','particlesets')

        if ps is None:
            return None
        #end if

        ions = None
        elns = None
        ion_list = []
        for name,p in ps.iteritems():
            if 'ionid' in p:
                ion_list.append(p)
            elif name.startswith('e'):
                elns = p
            #end if
        #end for
        if len(ion_list)==0: #try to identify ions by positive charged groups
            for name,p in ps.iteritems():
                if 'groups' in p:
                    for g in p.groups:
                        if 'charge' in g and g.charge>0:
                            ion_list.append(p)
                            break
                        #end if
                    #end for
                #end if
            #end for
        #end if
        if len(ion_list)==1:
            ions = ion_list[0]
        elif len(ion_list)>1:
            self.error('ability to handle multiple ion particlesets has not been implemented')
        #end if

        if elns==None:
            self.error('could not find electron particleset')
        #end if
        if ions==None:
            self.error('could not find ion particleset')
        #end if

        elem = None
        if 'ionid' in ions:
            if isinstance(ions.ionid,str):
                elem = [ions.ionid]
            else:
                elem = list(ions.ionid)
            #end if
            pos  = ions.position
        elif 'size' in ions and ions.size==1:
            elem = [ions.groups.list()[0].name]
            pos  = [[0,0,0]]
        elif 'groups' in ions:
            elem = []
            pos  = []
            for group in ions.groups:
                if 'position' in group:
                    nions = len(group.position)
                    elem.extend(nions*[group.name])
                    pos.extend(list(group.position))
                #end if
            #end for
            if len(elem)==1:
                elem = None
                pos  = None
            else:
                elem = array(elem)
                pos  = array(pos)
                order = elem.argsort()
                elem = elem[order]
                pos  = pos[order]
            #end if
        #end if
        if elem is None:
            self.error('could not read ions from ion particleset')
        #end if

        if axes==None:
            center = (0,0,0)
        else:
            md = input._metadata
            if 'position' in md and 'condition' in md['position'] and md['position']['condition']==1:
                pos = dot(pos,axes)
            #end if
            center = axes.sum(0)/2
        #end if
        
        structure = Structure(axes=axes,elem=elem,pos=pos,center=center,units='B')
        
        valency = dict()
        ion_charge = 0
        for name,element in ions.groups.iteritems():
            if 'charge' in element:
                valence = element.charge
            elif 'valence' in element:
                valence = element.valence
            elif 'atomic_number' in element:
                valence = element.atomic_number
            else:
                self.error('could not identify valency of '+name)
            #end if
            valency[name] = valence
            count = list(elem).count(name)
            ion_charge += valence*count
        #end for

        net_spin = 0
        eln_charge = 0
        for spin,eln in elns.groups.iteritems():
            if spin[0]=='u':
                net_spin+=eln.size
            elif spin[0]=='d':
                net_spin-=eln.size
            #end if
            eln_charge += eln.charge*eln.size
        #end if

        net_charge = ion_charge + eln_charge

        system = PhysicalSystem(structure,net_charge,net_spin,**valency) 
        
        return system
    #end def return_system


    def get_ion_particlesets(self):
        ions = obj()
        ps = self.get('particlesets')
        #try to identify ions by positive charged groups
        for name,p in ps.iteritems():
            if name.startswith('ion') or name.startswith('atom'):
                ions[name] = p
            elif 'groups' in p:
                for g in p.groups:
                    if 'charge' in g and g.charge>0:
                        ions[name] = p
                        break
                    #end if
                #end for
            #end if
        #end for
        return ions
    #end def get_ion_particlesets


    def get_pp_files(self):
        pp_files = []
        h = self.get('hamiltonian')
        if h != None:
            pp = None
            if 'pairpots' in h:
                for pairpot in h.pairpots:
                    if 'type' in pairpot and pairpot.type=='pseudo':
                        pp = pairpot
                    #end if
                #end for
            elif 'pairpot' in h and 'type' in h.pairpot and h.pairpot.type=='pseudo':
                pp = h.pairpot
            #end if
            if pp!=None:
                if 'pseudo' in pp and 'href' in pp.pseudo:
                    pp_files.append(pp.pseudo.href)
                elif 'pseudos' in pp:
                    for pseudo in pp.pseudos:
                        if 'href' in pseudo:
                            pp_files.append(pseudo.href)
                        #end if
                    #end for
                #end if
            #end if
        #end if
        return pp_files
    #end def get_pp_files


    def remove_physical_system(self):
        qs = self.simulation.qmcsystem
        if 'simulationcell' in qs:
            del qs.simulationcell
        #end if
        if 'particlesets' in qs:
            del qs.particlesets
        #end if
        for name in qs.keys():
            if isinstance(qs[name],particleset):
                del qs[name]
            #end if
        #end for
        self.replace('ion0','i')
    #end def remove_physical_system


    def bundle(self,inputs,filenames):
        return BundledQmcpackInput(inputs,filenames)
    #end def bundle

    
    def trace(self,quantity,values):
        return TracedQmcpackInput(quantity,values,self)
    #end def trace

    
    def twist_average(self,twistnums):
        return self.trace('twistnum',twistnums)
    #end def twist_average
#end class QmcpackInput



class BundledQmcpackInput(SimulationInput):
    
    def __init__(self,inputs,filenames):
        self.inputs = obj()
        for input in inputs:
            self.inputs.append(input)
        #end for
        self.filenames = filenames
    #end def __init__


    def get_output_info(self,*requests):
        outfiles = []

        for index,input in self.inputs.iteritems():
            outfs = input.get_output_info('outfiles')
            infile = self.filenames[index]
            outfile= infile.rsplit('.',1)[0]+'.g'+str(index).zfill(3)+'.qmc'
            outfiles.append(infile)
            outfiles.append(outfile)
            for outf in outfs:
                prefix,rest = outf.split('.',1)
                outfiles.append(prefix+'.g'+str(index).zfill(3)+'.'+rest)
            #end for
        #end for

        values = []
        for req in requests:
            if req=='outfiles':
                values.append(outfiles)
            else:
                values.append(None)
            #end if
        #end for
        if len(values)==1:
            return values[0]
        else:
            return values
        #end if
    #end def get_output_info

        
    def generate_filenames(self):
        None
    #end def generate_filenames
        

    def write(self,filepath=None):
        if filepath!=None and not 'filenames' in self:
            infile = os.path.split(filepath)[1]
            if not infile.endswith('.xml'):
                infile+='.xml'
            #end if
            self.generate_filenames(infile)
        #end if
        if filepath==None:
            c = ''
            for i in range(len(self.inputs)):
                c += self.filenames[i]+'\n'
            #end for
            return c
        else:
            path,file  = os.path.split(filepath)
            #if file!=self.filenames[-1]:
            #    self.error('main filenames do not match\n  internal: '+self.filenames[-1]+'\n  inputted: '+file)
            ##end if
            c = ''
            for i in range(len(self.inputs)):
                input = self.inputs[i]
                bfile = self.filenames[i]
                c += bfile+'\n'
                bfilepath = os.path.join(path,bfile)
                input.write(bfilepath)
            #end for
            fobj = open(filepath,'w')
            fobj.write(c)
            fobj.close()
        #end if
    #end def write
#end class BundledQmcpackInput



class TracedQmcpackInput(BundledQmcpackInput):
    def __init__(self,quantity=None,values=None,input=None):
        self.quantities = obj()
        self.variables = obj()
        self.inputs = obj()
        if quantity!=None and values!=None and input!=None:
            self.bundle_inputs(quantity,values,input)
        #end if
    #end def __init__

    def bundle_inputs(self,quantity,values,input):
        range = len(self.inputs),len(self.inputs)+len(values)
        self.quantities.append(obj(quantity=quantity,range=range))
        for value in values:
            inp = input.copy()
            qhost = inp.get_host(quantity)                               
            if qhost!=None:
                qhost[quantity] = value
            else:
                self.error('quantity '+quantity+' was not found in '+input.__class__.__name__)
            #end if
            self.variables.append(obj(quantity=quantity,value=value))
            self.inputs.append(inp)
        #end for
    #end def bundle_inputs


    def generate_filenames(self,infile):
        prefix,ext = infile.split('.',1)
        if not ext.endswith('xml'):
            ext+='.xml'
        #end if
        self.filenames = []
        for i in range(len(self.variables)):
            var = self.variables[i]
            q = var.quantity
            v = var.value
            bfile = prefix+'.g'+str(i).zfill(3)+'.'+q+'_'+str(v)+'.'+ext
            self.filenames.append(bfile)
        #end if
        self.filenames.append(prefix+'.in')
    #end def generate_filenames
#end class TracedQmcpackInput



#def generate_jastrows(descriptors,system=None):
#    js = []
#    for d in descriptors:
#        js.append(generate_jastrow(d,system))
#    #end for
#    jastrows = make_collection(js)
#    return jastrows
##end def generate_jastrows


def generate_jastrow(descriptor,*args,**kwargs):
    keywords = set(['function','size','rcut','elements','coeff','cusp','ename',
                    'iname','spins','density','Buu','Bud'])
    if not 'system' in kwargs:
        kwargs['system'] = None
    #end if
    system = kwargs['system']
    del kwargs['system']
    if isinstance(descriptor,str):
        descriptor = [descriptor]
    #end if
    ikw=0
    for i in range(len(descriptor)):
        if descriptor[i] in keywords:
            break
        #end if
        ikw += 1
    #end for
    dargs = descriptor[1:ikw]
    if len(dargs)>0:
        args = dargs
    #end if
    for i in range(ikw,len(descriptor)):
        d = descriptor[i]
        if isinstance(d,str):
            if d in keywords:
                kwargs[d] = descriptor[i+1]
            else:
                QmcpackInput.class_warn('keyword {0} is unrecognized and ignored\n  valid options are: {1}'.format(d,str(keywords)),'generate_jastrow')
            #end if
        #end if
    #end for
    kwargs['system'] = system
    jtype = descriptor[0]
    if jtype=='J1':
        jastrow = generate_jastrow1(*args,**kwargs)
    elif jtype=='J2':
        jastrow = generate_jastrow2(*args,**kwargs)
    elif jtype=='J3':
        jastrow = generate_jastrow3(*args,**kwargs)
    else:
        QmcpackInput.class_error('jastrow type unrecognized: '+jtype)
    #end if
    return jastrow
#end def generate_jastrow



def generate_jastrows(jastrows,system=None,return_list=False):
    jin = []
    if isinstance(jastrows,str):
        jorders = set(jastrows.replace('generate',''))
        if '1' in jorders:
            jin.append(
                generate_jastrow('J1','bspline',8,system=system)
                )
        #end if
        if '2' in jorders:
            jin.append(
                generate_jastrow('J2','bspline',8,system=system)
                )
        #end if
        if '3' in jorders:
            jin.append(
                generate_jastrow('J3','polynomial',4,4,5.0,system=system)
                )
        #end if
        if len(jin)==0:
            QmcpackInput.class_error('jastrow generation requested but no orders specified (1,2,and/or 3)')
        #end if
    else:
        jset = set(['J1','J2','J3'])
        for jastrow in jastrows:
            if isinstance(jastrow,QIxml):
                jin.append(jastrow)
            elif isinstance(jastrow,dict) or isinstance(jastrow,obj):
                jdict = dict(**jastrow)
                if not 'type' in jastrow:
                    QmcpackInput.class_error("could not determine jastrow type from input\n  field 'type' must be 'J1', 'J2', or 'J3'\n  object you provided: "+str(jastrow))
                #end if
                jtype = jdict['type']
                if not jtype in jset:
                    QmcpackInput.class_error("invalid jastrow type provided\n  field 'type' must be 'J1', 'J2', or 'J3'\n  object you provided: "+str(jdict))
                #end if
                del jdict['type']
                if 'system' in jdict:
                    jsys = jdict['system']
                    del jdict['system']
                else:
                    jsys = system
                #end if
                jin.append(generate_jastrow(jtype,system=jsys,**jdict))
                del jtype
                del jsys
            elif jastrow[0] in jset:
                jin.append(generate_jastrow(jastrow,system=system))
            else:
                QmcpackInput.class_error('starting jastrow unrecognized:\n  '+str(jastrow))
            #end if
        #end for
    #end if
    if return_list:
        return jin
    else:
        wf = wavefunction(jastrows=jin)
        wf.pluralize()
        return wf.jastrows
    #end if
#end def generate_jastrows




def generate_jastrow1(function='bspline',size=8,rcut=None,coeff=None,cusp=0.,ename='e',iname='ion0',elements=None,system=None,**elemargs):
    noelements = elements is None
    nosystem   = system is None
    noelemargs = len(elemargs)==0
    isopen = False
    if noelements and nosystem and noelemargs:
        QmcpackInput.class_error('must specify elements or system','generate_jastrow1')
    #end if
    if noelements:
        elements = []
    #end if
    if not nosystem:
        elements.extend(list(set(system.structure.elem)))
        isopen = system.structure.volume() is None
    #end if
    if not noelemargs:
        elements.extend(elemargs.keys())
    #end if
    # remove duplicate elements
    eset = set()
    elements = [ e for e in elements if e not in eset and not eset.add(e) ]     
    corrs = []
    for i in range(len(elements)):
        element = elements[i]
        if cusp is 'Z':
            QmcpackInput.class_error('need to implement Z cusp','generate_jastrow1')
        else:
            lcusp  = cusp
        #end if
        lrcut  = rcut
        lcoeff = size*[0]
        if coeff!=None:
            if element in coeff:
                lcoeff = coeff[element]
            else:
                lcoeff = coeff[i]
            #end if
        #end if
        if element in elemargs:
            v = elemargs[element]
            if 'cusp' in v:
                lcusp = v['cusp']
            #end if
            if 'rcut' in v:
                lrcut = v['rcut']
            #end if
            if 'size' in v and not 'coeff' in v:
                lcoeff = v['size']*[0]
            #end if
            if 'coeff' in v:
                lcoeff = v['coeff']
            #end if
        #end if
        corr = correlation(
            elementtype = element,
            size        = len(lcoeff),
            cusp        = cusp,
            coefficients=section(
                id    = ename+element,
                type  = 'Array',
                coeff = lcoeff
                )
            )            
        if lrcut!=None:
            corr.rcut = lrcut
        elif isopen:
            QmcpackInput.class_error('rcut must be provided for an open system','generate_jastrow1')
        #end if
        corrs.append(corr)
    #end for
    j1 = jastrow1(
        name         = 'J1',
        type         = 'One-Body',
        function     = function,
        source       = iname,
        print_       = True,
        correlations = corrs
        )
    return j1
#end def generate_jastrow1



def generate_bspline_jastrow2(size=8,rcut=None,coeff=None,spins=('u','d'),density=None,system=None,init='rpa'):
    if coeff is None and system is None and (init=='rpa' and density is None or rcut is None):
        QmcpackInput.class_error('rcut and density or system must be specified','generate_bspline_jastrow2')
    #end if
    set_rcut = rcut!=None
    if system!=None:
        cell = system.structure
        volume = cell.volume()
        if volume is None: #assume it is an open system
            if rcut is None:
                QmcpackInput.class_error('rcut must be provided for an open system','generate_bspline_jastrow2')
            #end if
            if init=='rpa':
                init = 'zero'
            #end if
        else:
            if rcut is None:
                rcut = cell.rmin()
            #end if
            nelectrons = system.particles.count_electrons()
            density = nelectrons/volume
        #end if
    elif init=='rpa':
        init = 'zero'
    #end if
    if coeff is None:
        if init=='rpa':
            wp = sqrt(4*pi*density)
            dr = rcut/size
            r = .02 + dr*arange(size)
            uuc = .5/(wp*r)*(1.-exp(-r*sqrt(wp/2)))*exp(-(2*r/rcut)**2)
            udc = .5/(wp*r)*(1.-exp(-r*sqrt(wp))  )*exp(-(2*r/rcut)**2)
            coeff = [uuc,udc]
        elif init=='zero' or init==0:
            coeff = [size*[0],size*[0]]
        else:
            QmcpackInput.class_error(str(init)+' is not a valid value for parameter init\n  valid options are: rpa, zero','generate_bspline_jastrow2')
        #end if
    elif len(coeff)!=2:
        QmcpackInput.class_error('must provide 2 sets of coefficients (uu,ud)','generate_bspline_jastrow2')
    #end if
    size = len(coeff[0])
    uname,dname = spins
    uuname = uname+uname
    udname = uname+dname
    corrs = [
        correlation(speciesA=uname,speciesB=uname,size=size,
                    coefficients=section(id=uuname,type='Array',coeff=coeff[0])),
        correlation(speciesA=uname,speciesB=dname,size=size,
                    coefficients=section(id=udname,type='Array',coeff=coeff[1]))
        ]
    if set_rcut:
        for corr in corrs:
            corr.rcut=rcut
        #end for
    #end if
    j2 = jastrow2(
        name = 'J2',type='Two-Body',function='bspline',print_=True,
        correlations = corrs
        )
    return j2
#end def generate_bspline_jastrow2


def generate_pade_jastrow2(Buu=None,Bud=None,spins=('u','d'),system=None):
    if Buu is None:
        Buu = 2.0
    #end if
    if Bud is None:
        Bud = float(Buu)
    #end if
    uname,dname = spins
    uuname = uname+uname
    udname = uname+dname
    cuu = var(id=uuname+'_b',name='B',value=Buu)
    cud = var(id=udname+'_b',name='B',value=Bud)
    corrs = [
        correlation(speciesA=uname,speciesB=uname,
                    vars=[cuu]),
        correlation(speciesA=uname,speciesB=dname,
                    vars=[cud])
        ]
    j2 = jastrow2(
        name = 'J2',type='Two-Body',function='pade',
        correlations = corrs
        )
    return j2
#end def generate_pade_jastrow2



def generate_jastrow2(function='bspline',*args,**kwargs):
    if not 'spins' in kwargs:
        kwargs['spins'] = ('u','d')
    #end if
    spins = kwargs['spins']
    if not isinstance(spins,tuple) and not isinstance(spins,list):
        QmcpackInput.class_error('spins must be a list or tuple of u/d spin names\n  you provided: '+str(spins))
    #end if
    if len(spins)!=2:
        QmcpackInput.class_error('name for up and down spins must be specified\n  you provided: '+str(spins))
    #end if
    if not isinstance(function,str):
        QmcpackInput.class_error('function must be a string\n  you provided: '+str(function),'generate_jastrow2')
    #end if
    if function=='bspline':
        j2 = generate_bspline_jastrow2(*args,**kwargs)
    elif function=='pade':
        j2 = generate_pade_jastrow2(*args,**kwargs)
    else:
        QmcpackInput.class_error('function is invalid\n  you provided: {0}\n  valid options are: bspline or pade'.format(function),'generate_jastrow2')
    #end if
    return j2
#end def generate_jastrow2



def generate_jastrow3(function='polynomial',esize=4,isize=4,rcut=5.,coeff=None,iname='ion0',spins=('u','d'),elements=None,system=None):
    if elements is None and system is None:
        QmcpackInput.class_error('must specify elements or sytem to generate jastrow3')
    elif elements is None:
        elements = list(system.structure.elem)
    #end if
    if coeff!=None:
        QmcpackInput.class_error('handling coeff is not yet implemented for generate jastrow3')
    #end if
    if len(spins)!=2:
        QmcpackInput.class_error('must specify name for up and down spins\n  provided: '+str(spins))
    #end if
    uname,dname = spins
    uuname = uname+uname
    udname = uname+dname
    corrs=[]
    for element in elements:
        corrs.append(
            correlation(
                especies1=uname,especies2=uname,ispecies=element,esize=esize,
                isize=isize,rcut=rcut,
                coefficients=section(id=uuname+element,type='Array',optimize=True))
            )
        corrs.append(
            correlation(
                especies1=uname,especies2=dname,ispecies=element,esize=esize,
                isize=isize,rcut=rcut,
                coefficients=section(id=udname+element,type='Array',optimize=True))
            )
    #end for
    jastrow = jastrow3(
        name = 'J3',type='eeI',function=function,print_=True,source=iname,
        correlations = corrs
        )
    return jastrow
#end def generate_jastrow3



def count_jastrow_params(jastrows):
    if isinstance(jastrows,QIxml):
        jastrows = [jastrows]
    #end if
    params = 0
    for jastrow in jastrows:
        name = jastrow.name
        if 'type' in jastrow:
            type = jastrow.type.lower()
        else:
            type = ''
        #end if
        jastrow.pluralize()
        if name=='J1' or type=='one-body':
            for correlation in jastrow.correlations:
                params += correlation.size
            #end for
        elif name=='J2' or type=='two-body':
            for correlation in jastrow.correlations:
                params += correlation.size
            #end for
        elif name=='J3' or type=='eeI':
            for correlation in jastrow.correlations:
                params += correlation.esize
                params += correlation.isize
            #end for
        #end if
    #end for
    return params
#end def count_jastrow_params


opt_map = dict(linear=linear,cslinear=cslinear)
def generate_opt(method,
                 repeat           = 1,
                 energy           = None,
                 rw_variance      = None,
                 urw_variance     = None,
                 params           = None,
                 jastrows         = None,
                 processes        = None,
                 walkers_per_proc = None,
                 threads          = None,
                 blocks           = 2000,
                 #steps            = 5,
                 decorr           = 10,
                 min_walkers      = None, #use e.g. 128 for gpu's
                 timestep         = .5,
                 nonlocalpp       = False,
                 sample_factor    = 1.0):
    if not method in opt_map:
        QmcpackInput.class_error('section cannot be generated for optimization method '+method)
    #end if
    if energy is None and rw_variance is None and urw_variance is None:
        QmcpackInput.class_error('at least one cost parameter must be specified\n options are: energy, rw_variance, urw_variance')
    #end if
    if params is None and jastrows is None:
        QmcpackInput.class_error('must provide either number of opt parameters (params) or a list of jastrow objects (jastrows)')
    #end if
    if processes is None:
        QmcpackInput.class_error('must specify total number of processes')
    elif walkers_per_proc is None and threads is None:
        QmcpackInput.class_error('must specify walkers_per_proc or threads')
    #end if

    if params is None:
        params = count_jastrow_params(jastrows)
    #end if
    samples = max(100000,100*params**2)
    samples = int(round(sample_factor*samples))
    samples_per_proc = int(round(float(samples)/processes))

    if walkers_per_proc is None:
        walkers = 1
        walkers_per_proc = threads
    else:
        walkers = walkers_per_proc
    #end if
    tot_walkers = processes*walkers_per_proc
    if min_walkers!=None:
        tot_walkers = max(min_walkers,tot_walkers)
        walkers = int(ceil(float(tot_walkers)/processes-.001))
        if threads!=None and mod(walkers,threads)!=0:
            walkers = threads*int(ceil(float(walkers)/threads-.001))
        #end if
    #end if
    #blocks = int(ceil(float(decorr*samples)/(steps*tot_walkers)))
    blocks = min(blocks,samples_per_proc*decorr)

    opt = opt_map[method]()

    opt.set(
        walkers    = walkers,
        blocks     = blocks,
        #steps      = steps,
        samples    = samples,
        substeps   = decorr,
        timestep   = timestep,
        nonlocalpp = nonlocalpp,
        stepsbetweensamples = 1
        )
    if energy!=None:
        opt.energy = energy
    #end if
    if rw_variance!=None:
        opt.reweightedvariance = rw_variance
    #end if
    if urw_variance!=None:
        opt.unreweightedvariance = urw_variance
    #end if
    
    opt.incorporate_defaults(elements=True)

    if repeat>1:
        opt = loop(max=repeat,qmc=opt)
    #end if

    return opt
#end def generate_opt


def generate_opts(opt_reqs,**kwargs):
    opts = []
    for opt_req in opt_reqs:
        opts.append(generate_opt(*opt_req,**kwargs))
    #end for
    return opts
#end def generate_opts



def generate_qmcpack_input(selector,*args,**kwargs):
    if selector=='basic':
        return generate_basic_input(**kwargs)
    elif selector=='opt_jastrow':
        return generate_opt_jastrow_input(*args,**kwargs)
    else:
        QmcpackInput.class_error('selection '+str(selector)+' has not been implemented for qmcpack input generation')
    #end if
#end def generate_qmcpack_input




def generate_basic_input(id           = 'qmc',
                         series       = 0,
                         purpose      = '',
                         seed         = None,
                         bconds       = None,
                         remove_cell  = False,
                         meshfactor   = 1.0,
                         twistnum     = 0, 
                         orbitals_h5  = 'MISSING.h5',
                         system       = None,
                         pseudos      = None,
                         jastrows     = None,
                         corrections  = None,
                         observables  = None,
                         estimators   = None,
                         calculations = None):
    if bconds is None:
        bconds = array(['p','p','p'])
    else:
        bconds = array(tuple(bconds))
    #end if
    if pseudos is None:
        pseudopotentials = []
    else:
        pseudopotentials = pseudos
    #end if
    if jastrows is None:
        jastrows = [('J1','bspline','size',8),
                    ('J2','bspline','size',8)]
    #end if
    if corrections is None and tuple(bconds)==tuple('ppp'):
        corrections = ['mpc','chiesa']
    else:
        corrections = []
    #end if
    if observables is None:
        #observables = ['localenergy']
        observables = []
    #end if
    if estimators is None:
        estimators = []
    #end if
    if calculations is None:
        calculations = []
    #end if

    metadata = QmcpackInput.default_metadata.copy()
    qi = QmcpackInput(
        metadata,
        simulation(
            project = section(
                id=id,
                series = series,
                application = application()
                ),
            qmcsystem = section(
                simulationcell = section(bconds=array(bconds)),
                wavefunction   = section(
                    determinantset = section(
                        twistnum = twistnum,
                        meshfactor = meshfactor,
                        href       = orbitals_h5,
                        slaterdeterminant = section(
                            determinants = collection(
                                updet   = determinant(
                                    id='updet'  ,#ref='updet',
                                    occupation=section(mode='ground',spindataset=0)
                                    ),
                                downdet = determinant(
                                    id='downdet',#ref='downdet',
                                    occupation=section(mode='ground',spindataset=0)
                                    ),
                                )
                            )
                        )
                    ),
                hamiltonian = section(
                    constant = section(),
                    pairpots = collection(
                        coulomb(name='ElecElec'),
                        pseudopotential(name='PseudoPot')
                        )
                    )
                )
            )
        )


    if seed!=None:
        qi.simulation.random = random(seed=seed)
    #end if
    if system!=None:
        qi.incorporate_system(system)
    #end if
    ham = qi.simulation.qmcsystem.hamiltonian
    estimators = estimators + corrections + observables
    ests=[]
    for estimator in estimators:
        estname = estimator.lower().replace(' ','_').replace('-','_').replace('__','_')
        est=None
        if estname=='mpc':
            ham.pairpots.mpc = mpc()
        elif estname=='chiesa':
            est = chiesa(name='KEcorr')
        elif estname=='localenergy':
            est = localenergy(name='LocalEnergy')
        elif estname=='energydensity':
            est = energydensity(
                type='EnergyDensity',name='EDvoronoi',dynamic='e',static='ion0',
                spacegrid = spacegrid(coord='voronoi')
                )
        else:
            QmcpackInput.class_error('estimator '+estimator+' has not yet been enabled in generate_basic_input')
        #end if
        if est!=None:
            ests.append(est)
        #end if
    #end for
    if len(ests)>0:
        ham.estimators = make_collection(ests)
    #end if
    ham.pluralize()
    pseudos = ham.get('pseudo')
    for ppfile in pseudopotentials:
        #element = ppfile[0:2]
        element = ppfile.split('.')[0]
        if len(element)>2:
            element = element[0:2]
        #end if
        if element in pseudos:
            pseudos[element].href = ppfile
        else:
            QmcpackInput.class_error('hamiltonian is malformed\n  '+element+' is not in pseudopotential list')
        #end if
    #end for

    wf = qi.get('wavefunction')
    wf.jastrows = generate_jastrows(jastrows,system)

    for calculation in calculations:
        if isinstance(calculation,loop):
            calc = calculation.qmc
        else:
            calc = calculation
        #end if
        has_localenergy = False
        has_estimators = 'estimators' in calc
        if has_estimators:
            estimators = calc.estimators
            if not isinstance(estimators,collection):
                estimators = make_collection(estimators)
            #end if
            has_localenergy = 'localenergy' in estimators or 'LocalEnergy' in estimators
        else:
            estimators = collection()
        #end if
        if not has_localenergy:
            estimators.localenergy = localenergy(name='LocalEnergy')
            calc.estimators = estimators
        #end if
    #end for
    calc_coll = make_collection(calculations)
    qi.simulation.calculations = calc_coll.copy()

    qi.incorporate_defaults(elements=False,overwrite=False,propagate=True)

    if remove_cell:
        qi.remove_physical_system()
    #end if

    return qi
#end def generate_basic_input




def generate_opt_jastrow_input(id  = 'qmc',
                               series           = 0,
                               purpose          = '',
                               seed             = None,
                               bconds           = None,
                               remove_cell      = False,
                               meshfactor       = 1.0,
                               twistnum         = 0, 
                               orbitals_h5      = 'MISSING.h5',
                               system           = None,
                               pseudos          = None,
                               jastrows         = 'generateJ12',
                               corrections      = None,
                               observables      = None,
                               processes        = None,
                               walkers_per_proc = None,
                               threads          = None,
                               decorr           = 10,
                               min_walkers      = None, #use e.g. 128 for gpu's
                               timestep         = 0.5,
                               nonlocalpp       = False,
                               sample_factor    = 1.0,
                               opt_calcs        = None):
    jastrows = generate_jastrows(jastrows,system)

    if opt_calcs is None:
        opt_calcs = [
            ('linear', 4,  0,  0, 1.0),
            ('linear', 4, .8, .2,   0)
            ]
    #end if
    opts = []
    for opt_calc in opt_calcs:
        if isinstance(opt_calc,QIxml):
            opts.append(opt_calc)
        elif len(opt_calc)==5:
            if opt_calc[0] in opt_map:
                opts.append(
                    generate_opt(
                        *opt_calc,
                         jastrows         = jastrows,
                         processes        = processes,
                         walkers_per_proc = walkers_per_proc,
                         threads          = threads,
                         decorr           = decorr,
                         min_walkers      = min_walkers,
                         timestep         = timestep,
                         nonlocalpp       = nonlocalpp,
                         sample_factor    = sample_factor
                         )
                    )
            else:
                QmcpackInput.class_error('optimization method '+opt_calc[0]+' has not yet been implemented')
            #end if
        else:
            QmcpackInput.class_error('optimization calculation is ill formatted\n  opt calc provided: \n'+str(opt_calc))
        #end if
    #end if

    input = generate_basic_input(
        id           = id,
        series       = series       ,
        purpose      = purpose      ,
        seed         = seed         ,
        bconds       = bconds       ,
        remove_cell  = remove_cell  ,
        meshfactor   = meshfactor   ,
        twistnum     = twistnum     ,
        orbitals_h5  = orbitals_h5  ,
        system       = system       ,
        pseudos      = pseudos      ,
        jastrows     = jastrows     ,
        corrections  = corrections  ,
        observables  = observables  ,
        calculations = opts
        )

    return input
#end def generate_opt_jastrow_input







if __name__=='__main__':

    filepath = './example_input_files/c_boron/qmcpack.in.xml'

    element_joins=['qmcsystem']
    element_aliases=dict(loop='qmc')
    xml = XMLreader(filepath,element_joins,element_aliases,warn=False).obj
    xml.condense()

    qi = QmcpackInput()
    qi.read(filepath)


    s = qi.simulation
    q = s.qmcsystem
    c = s.calculations
    h = q.hamiltonian
    p = q.particlesets
    w = q.wavefunction
    j = w.jastrows
    co= j.J1.correlations.B.coefficients


    qi.write('./output/qmcpack.in.xml')

    #qi.condensed_name_report()
    #exit()


    test_ret_system    = 1
    test_gen_input     = 0
    test_difference    = 0
    test_moves         = 0
    test_defaults      = 0
    test_substitution  = 0
    test_generation    = 0


    if test_ret_system:
        from structure import generate_structure
        from physical_system import PhysicalSystem
        
        system = PhysicalSystem(
            structure = generate_structure('diamond','fcc','Ge',(2,2,2),scale=5.639,units='A'),
            net_charge = 1,
            net_spin   = 1,
            Ge = 4
        )

        gi = generate_qmcpack_input('basic',system=system)
        
        rsys = gi.return_system()

        print rsys

    #end if


    if test_gen_input:
        from structure import generate_structure
        from physical_system import PhysicalSystem
        
        system = PhysicalSystem(
            structure = generate_structure('diamond','fcc','Ge',(2,2,2),scale=5.639,units='A'),
            net_charge = 1,
            net_spin   = 1,
            Ge = 4
        )

        gi = generate_qmcpack_input('basic',system=system)
        
        print gi

        print gi.write()
    #end if



    if test_difference:
        tstep = QmcpackInput('./example_input_files/luke_tutorial/diamond-dmcTsteps.xml')
        opt   = QmcpackInput('./example_input_files/luke_tutorial/opt-diamond.xml')

        different,diff,d1,d2 = tstep.difference(tstep)
        different,diff,d1,d2 = tstep.difference(opt)
        
    #end if



    if test_moves:
        print 50*'='
        sim = qi.simulation
        print repr(sim)
        print repr(sim.qmcsystem)
        print 50*'='
        qi.move(particleset='simulation')
        print repr(sim)
        print repr(sim.qmcsystem)
        print 50*'='
        qi.standard_placements()
        print repr(sim)
        print repr(sim.qmcsystem)

        qi.pluralize()
    #end if



    if test_defaults:
        q=QmcpackInput(
            simulation(
                qmcsystem=section(
                    simulationcell = section(),
                    wavefunction = section(),
                    hamiltonian = section()             
                    ),
                calculations = [
                    cslinear(),
                    vmc(),
                    dmc()
                    ]
                )            
            )

        #q.simulation = simulation()

        q.incorporate_defaults(elements=True)

        print q
    #end if


    if test_substitution:
        q = qi.copy()

        q.remove('simulationcell','particleset','wavefunction')
        q.write('./output/qmcpack.remove.xml')
        q.include_xml('./example_input_files/energy_density/Si.ptcl.xml',replace=False)
        q.include_xml('./example_input_files/energy_density/Si.wfs.xml',replace=False)
        q.write('./output/qmcpack.replace.xml')

        qnj = QmcpackInput()
        qnj.read('./example_input_files/jastrowless/opt_jastrow.in.xml')

        qnj.generate_jastrows(size=6)
        qnj.write('./output/jastrow_gen.in.xml')

    #end if
    


    if test_generation:

        q=QmcpackInput(
            meta(
                lattice    = {'units':'bohr'},
                reciprocal = {'units':'2pi/bohr'},
                ionid      = {'datatype':'stringArray'},
                position   = {'datatype':'posArray', 'condition':0}
                ),
            simulation(
                project = section(
                    id='C16B',
                    series = 0,
                    application = section(
                        name = 'qmcapp',
                        role = 'molecu',
                        class_ = 'serial',
                        version = .2
                        ),
                    host = 'kraken',
                    date = '3 May 2012',
                    user = 'jtkrogel'
                    ),
                random = section(seed=13),
                qmcsystem = section(
                    simulationcell = section(
                        name = 'global',
                        lattice = array([[1,1,0],[1,0,1],[0,1,1]]),
                        reciprocal = array([[1,1,-1],[1,-1,1],[-1,1,1]]),
                        bconds = 'p p p',
                        LR_dim_cutoff = 15            
                        ),
                    particlesets = [
                        particleset(
                            name = 'ion0',
                            size = 32,
                            groups=[
                                group(
                                    name='C',
                                    charge=4.
                                    ),
                                group(
                                    name='B',
                                    charge = 3.
                                    )
                                ],
                            ionid = ['B','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C',
                                     'B','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C'],
                            position = array([
                                    [ 0.00, 0.00, 0.00],[ 1.68, 1.68, 1.68],[ 3.37, 3.37, 0.00],
                                    [ 5.05, 5.05, 1.68],[ 3.37, 0.00, 3.37],[ 5.05, 1.68, 5.05],
                                    [ 6.74, 3.37, 3.37],[ 8.42, 5.05, 5.05],[ 0.00, 3.37, 3.37],
                                    [ 1.68, 5.05, 5.05],[ 3.37, 6.74, 3.37],[ 5.05, 8.42, 5.05],
                                    [ 3.37, 3.37, 6.74],[ 5.05, 5.05, 8.42],[ 6.74, 6.74, 6.74],
                                    [ 8.42, 8.42, 8.42],[ 6.74, 6.74, 0.00],[ 8.42, 8.42, 1.68],
                                    [10.11,10.11, 0.00],[11.79,11.79, 1.68],[10.11, 6.74, 3.37],
                                    [11.79, 8.42, 5.05],[13.48,10.11, 3.37],[15.16,11.79, 5.05],
                                    [ 6.74,10.11, 3.37],[ 8.42,11.79, 5.05],[10.11,13.48, 3.37],
                                    [11.79,15.16, 5.05],[10.11,10.11, 6.74],[11.79,11.79, 8.42],
                                    [13.48,13.48, 6.74],[15.16,15.16, 8.42]])
                            ),
                        particleset(
                            name='e',
                            random = 'yes',
                            random_source = 'ion0',
                            groups=[
                                group(
                                    name='u',
                                    size=64,
                                    charge=-1
                                    ),
                                group(
                                    name='d',
                                    size=63,
                                    charge=-1
                                    )                    
                                ]
                            ),
                        ],
                    hamiltonians = [
                        hamiltonian(
                            name='h0',
                            type='generic',
                            target='e',
                            pairpots=[
                                pairpot(
                                    type = 'coulomb',
                                    name = 'ElecElec',
                                    source = 'e',
                                    target = 'e'
                                    ),
                                pairpot(
                                    type = 'pseudo',
                                    name = 'PseudoPot',
                                    source = 'ion0',
                                    wavefunction='psi0',
                                    format='xml',
                                    pseudos = [
                                        pseudo(
                                            elementtype='B',
                                            href='B.pp.xml'
                                            ),
                                        pseudo(
                                            elementtype='C',
                                            href='C.pp.xml'
                                            )
                                        ]
                                    )
                                ],
                            constant = section(
                                type='coulomb',
                                name='IonIon',
                                source='ion0',
                                target='ion0'
                                ),
                            estimators = [
                                estimator(
                                    type='energydensity',
                                    name='edvoronoi',
                                    dynamic='e',
                                    static='ion0',
                                    spacegrid = section(
                                        coord = 'voronoi'
                                        )
                                    ),
                                energydensity(
                                    name='edchempot',
                                    dynamic='e',
                                    static='ion0',
                                    spacegrid=spacegrid(
                                        coord='voronoi',
                                        min_part=-4,
                                        max_part=5
                                        )
                                    ),
                                estimator(
                                    type='energydensity',
                                    name='edcell',
                                    dynamic='e',
                                    static='ion0',
                                    spacegrid = section(
                                        coord = 'cartesian',
                                        origin = section(p1='zero'),
                                        axes   = (
                                            axis(label='x',p1='a1',scale=.5,grid='-1 (192) 1'),
                                            axis(label='y',p1='a2',scale=.5,grid='-1 (1) 1'),
                                            axis(label='z',p1='a3',scale=.5,grid='-1 (1) 1')
                                            )
        #                                axes   = collection(
        #                                    x = section(p1='a1',scale=.5,grid='-1 (192) 1'),
        #                                    y = section(p1='a2',scale=.5,grid='-1 (1) 1'),
        #                                    z = section(p1='a3',scale=.5,grid='-1 (1) 1')
        #                                    )
                                        )
                                    )
                                ]
                            )
                        ],
                    wavefunction = section(
                        name = 'psi0',
                        target = 'e',
                        determinantset = section(
                            type='bspline',
                            href='Si.pwscf.h5',
                            sort = 1,
                            tilematrix = array([[1,0,0],[0,1,0],[0,0,1]]),
                            twistnum = 0,
                            source = 'ion0',
                            slaterdeterminant = section(
                                determinants=[
                                    determinant(
                                        id='updet',
                                        size=64,
                                        occupation = section(
                                            mode='ground',
                                            spindataset=0
                                            )
                                        ),
                                    determinant(
                                        id='downdet',
                                        size=63,
                                        occupation = section(
                                            mode='ground',
                                            spindataset=1
                                            )
                                        )
                                    ]
                                ),
                            ),
                        jastrows = [
                            jastrow(
                                type='two-body',
                                name='J2',
                                function='bspline',
                                print_='yes',
                                correlations = [
                                    correlation(
                                        speciesA='u',
                                        speciesB='u',
                                        size=6,
                                        rcut=3.9,
                                        coefficients = section(
                                            id='uu',
                                            type='Array',
                                            coeff=[0,0,0,0,0,0]
                                            )
                                        ),
                                    correlation(
                                        speciesA='u',
                                        speciesB='d',
                                        size=6,
                                        rcut=3.9,
                                        coefficients = section(
                                            id='ud',
                                            type='Array',
                                            coeff=[0,0,0,0,0,0]
                                            )
                                        )
                                    ]
                                ),
                            jastrow(
                                type='one-body',
                                name='J1',
                                function='bspline',
                                source='ion0',
                                print_='yes',
                                correlations = [
                                    correlation(
                                        elementtype='C',
                                        size=6,
                                        rcut=3.9,
                                        coefficients = section(
                                            id='eC',
                                            type='Array',
                                            coeff=[0,0,0,0,0,0]
                                            )
                                        ),
                                    correlation(
                                        elementtype='B',
                                        size=6,
                                        rcut=3.9,
                                        coefficients = section(
                                            id='eB',
                                            type='Array',
                                            coeff=[0,0,0,0,0,0]
                                            )
                                        )
                                    ]
                                )
                            ]
                        ),

                    ),
                calculations=[
                    loop(max=4,
                         qmc=qmc(
                            method='cslinear',
                            move='pbyp',
                            checkpoint=-1,
                            gpu='no',
                            blocks = 3125,
                            warmupsteps = 5,
                            steps = 2,
                            samples = 80000,
                            timestep = .5,
                            usedrift = 'yes',
                            minmethod = 'rescale',
                            gevmethod = 'mixed',
                            exp0=-15,
                            nstabilizers = 5,
                            stabilizerscale = 3,
                            stepsize=.35,
                            alloweddifference=1e-5,
                            beta = .05,
                            bigchange = 5.,
                            energy = 0.,
                            unreweightedvariance = 0.,
                            reweightedvariance = 0.,
                            estimators=[
                                estimator(
                                    name='LocalEnergy',
                                    hdf5='no'
                                    )
                                ]
                            )
                        ),
                    qmc(
                        method = 'vmc',
                        multiple = 'no',
                        warp = 'no',
                        move = 'pbyp',
                        walkers = 1,
                        blocks = 2,
                        steps = 500,
                        substeps = 3,
                        timestep = .5,
                        usedrift = 'yes',
                        estimators=[
                            estimator(
                                name='LocalEnergy',
                                hdf5='yes'
                                )
                            ]
                        ),
                    qmc(
                        method='dmc',
                        move='pbyp',
                        walkers = 72,
                        blocks = 2,
                        steps = 50,
                        timestep = .01,
                        nonlocalmove = 'yes',
                        estimators=[
                            estimator(
                                name='LocalEnergy',
                                hdf5='no'
                                )
                            ]            
                        )
                    ]
                )
            )

        q.write('./output/gen.in.xml')





        #something broke this, check later
        exit()
        qs=QmcpackInput(
            simulation = section(
                project = section(
                    id='C16B',series = 0,
                    application = section(name='qmcapp',role='molecu',class_='serial',version=.2),
                    host='kraken',date='3 May 2012',user='jtkrogel'
                    ),
                random = section(seed=13),
                qmcsystem = section(
                    simulationcell = section(
                        name='global',bconds='p p p',lr_dim_cutoff=15,
                        lattice    = [[1,1,0] ,[1,0,1] ,[0,1,1]],
                        reciprocal = [[1,1,-1],[1,-1,1],[-1,1,1]],
                        ),
                    particlesets = collection(
                        ion0=particleset(
                            size=32,
                            groups=collection(
                                C = group(charge=4.),
                                B = group(charge=3.)),
                            ionid = ('B','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C',
                                     'B','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C'),
                            position = [[ 0.00, 0.00, 0.00],[ 1.68, 1.68, 1.68],[ 3.37, 3.37, 0.00],
                                        [ 5.05, 5.05, 1.68],[ 3.37, 0.00, 3.37],[ 5.05, 1.68, 5.05],
                                        [ 6.74, 3.37, 3.37],[ 8.42, 5.05, 5.05],[ 0.00, 3.37, 3.37],
                                        [ 1.68, 5.05, 5.05],[ 3.37, 6.74, 3.37],[ 5.05, 8.42, 5.05],
                                        [ 3.37, 3.37, 6.74],[ 5.05, 5.05, 8.42],[ 6.74, 6.74, 6.74],
                                        [ 8.42, 8.42, 8.42],[ 6.74, 6.74, 0.00],[ 8.42, 8.42, 1.68],
                                        [10.11,10.11, 0.00],[11.79,11.79, 1.68],[10.11, 6.74, 3.37],
                                        [11.79, 8.42, 5.05],[13.48,10.11, 3.37],[15.16,11.79, 5.05],
                                        [ 6.74,10.11, 3.37],[ 8.42,11.79, 5.05],[10.11,13.48, 3.37],
                                        [11.79,15.16, 5.05],[10.11,10.11, 6.74],[11.79,11.79, 8.42],
                                        [13.48,13.48, 6.74],[15.16,15.16, 8.42]]
                            ),
                        e=particleset(
                            random='yes',random_source='ion0',
                            groups = collection(
                                u=group(size=64,charge=-1),
                                d=group(size=63,charge=-1))
                            ),
                        ),
                    hamiltonian = section(
                        name='h0',type='generic',target='e',
                        pairpots=collection(
                            ElecElec = coulomb(name='ElecElec',source='e',target='e'),
                            PseudoPot = pseudopotential(
                                source='ion0',wavefunction='psi0',format='xml',
                                pseudos = collection(
                                    B = pseudo(href='B.pp.xml'),
                                    C = pseudo(href='C.pp.xml'))
                                )
                            ),
                        constant = section(type='coulomb',name='IonIon',source='ion0',target='ion0'),
                        estimators = collection(
                            edvoronoi = energydensity(
                                dynamic='e',static='ion0',spacegrid=section(coord ='voronoi')
                                ),
                            edchempot = energydensity(
                                dynamic='e',static='ion0',
                                spacegrid=section(coord='voronoi',min_part=-4,max_part=5)
                                ),
                            edcell = energydensity(
                                dynamic='e',static='ion0',
                                spacegrid = section(
                                    coord = 'cartesian',
                                    origin = section(p1='zero'),
                                    axes = collection(
                                        x = axis(p1='a1',scale=.5,grid='-1 (192) 1'),
                                        y = axis(p1='a2',scale=.5,grid='-1 (1) 1'),
                                        z = axis(p1='a3',scale=.5,grid='-1 (1) 1'))
                                    )
                                )
                            )
                        ),
                    wavefunction = section(
                        name = 'psi0',target = 'e',
                        determinantset = section(
                            type='bspline',href='Si.pwscf.h5',sort=1,twistnum=0,source='ion0',
                            tilematrix=(1,0,0,0,1,0,0,0,1),
                            slaterdeterminant = section(
                                determinants=collection(
                                    updet = determinant(
                                        size=64,
                                        occupation=section(mode='ground',spindataset=0)
                                        ),
                                    downdet = determinant(
                                        size=63,
                                        occupation = section(mode='ground',spindataset=1))
                                    )
                                ),
                            ),
                        jastrows = collection(
                            J2=jastrow2(
                                function='bspline',print_='yes',
                                correlations = collection(
                                    uu=correlation(
                                        speciesA='u',speciesB='u',size=6,rcut=3.9,
                                        coefficients = section(id='uu',type='Array',coeff=[0,0,0,0,0,0])
                                        ),
                                    ud=correlation(
                                        speciesA='u',speciesB='d',size=6,rcut=3.9,
                                        coefficients = section(id='ud',type='Array',coeff=[0,0,0,0,0,0])
                                        )
                                    )
                                ),
                            J1=jastrow1(
                                function='bspline',source='ion0',print_='yes',
                                correlations = collection(
                                    C=correlation(
                                        size=6,rcut=3.9,
                                        coefficients = section(
                                            id='eC',type='Array',coeff=[0,0,0,0,0,0])
                                        ),
                                    B=correlation(
                                        size=6,rcut=3.9,
                                        coefficients = section(id='eB',type='Array',coeff=[0,0,0,0,0,0])
                                        )
                                    )
                                )
                            )
                        ),
                    ),
                calculations=(
                    loop(max=4,
                         qmc=cslinear(
                            move='pbyp',checkpoint=-1,gpu='no',
                            blocks      = 3125,
                            warmupsteps = 5,
                            steps       = 2,
                            samples     = 80000,
                            timestep    = .5,
                            usedrift    = 'yes',
                            minmethod   = 'rescale',
                            gevmethod   = 'mixed',
                            exp0              = -15,
                            nstabilizers      =  5,
                            stabilizerscale   =  3,
                            stepsize          =  .35,
                            alloweddifference = 1e-5,
                            beta              = .05,
                            bigchange         = 5.,
                            energy               = 0.,
                            unreweightedvariance = 0.,
                            reweightedvariance   = 0.,
                            estimator = localenergy(hdf5='no')
                            )
                        ),
                    vmc(multiple='no',warp='no',move='pbyp',
                        walkers  =  1,
                        blocks   =  2,
                        steps    = 500,
                        substeps =  3,
                        timestep = .5,
                        usedrift = 'yes',
                        estimator = localenergy(hdf5='no')
                        ),
                    dmc(move='pbyp',
                        walkers  =  72,
                        blocks   =   2,
                        steps    =  50,
                        timestep = .01,
                        nonlocalmove = 'yes',
                        estimator = localenergy(hdf5='yes')
                        )
                    )
                )
            )

        qs.write('./output/simple.in.xml')


        est = qs.simulation.qmcsystem.hamiltonian.estimators
        sg = est.edcell.spacegrid
        print repr(est)

        exit()



















        q=QmcpackInput()
        q.simulation = section(
            project = section('C16B',0,
                application = section('qmcapp','molecu','serial',.2),
                host = 'kraken',
                date = '3 May 2012',
                user = 'jtkrogel'
                ),
            random = (13),
            qmcsystem = section(
                simulationcell = section(
                    units = 'bohr',
                    lattice = array([[1,1,0],[1,0,1],[0,1,1]]),
                    bconds = 'p p p',
                    LR_dim_cutoff = 15            
                    ),
                particlesets = [
                    particleset('ion0', ('C',4), ('B',3),
                        ionid = ['B','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C',
                                 'B','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C'],
                        position = array([
                                [ 0.00, 0.00, 0.00],[ 1.68, 1.68, 1.68],[ 3.37, 3.37, 0.00],
                                [ 5.05, 5.05, 1.68],[ 3.37, 0.00, 3.37],[ 5.05, 1.68, 5.05],
                                [ 6.74, 3.37, 3.37],[ 8.42, 5.05, 5.05],[ 0.00, 3.37, 3.37],
                                [ 1.68, 5.05, 5.05],[ 3.37, 6.74, 3.37],[ 5.05, 8.42, 5.05],
                                [ 3.37, 3.37, 6.74],[ 5.05, 5.05, 8.42],[ 6.74, 6.74, 6.74],
                                [ 8.42, 8.42, 8.42],[ 6.74, 6.74, 0.00],[ 8.42, 8.42, 1.68],
                                [10.11,10.11, 0.00],[11.79,11.79, 1.68],[10.11, 6.74, 3.37],
                                [11.79, 8.42, 5.05],[13.48,10.11, 3.37],[15.16,11.79, 5.05],
                                [ 6.74,10.11, 3.37],[ 8.42,11.79, 5.05],[10.11,13.48, 3.37],
                                [11.79,15.16, 5.05],[10.11,10.11, 6.74],[11.79,11.79, 8.42],
                                [13.48,13.48, 6.74],[15.16,15.16, 8.42]])
                        ),
                    particleset('e', ('u',-1,64), ('d',-1,63), random_source = 'ion0'),
                    ],
                hamiltonian = section('h0','e',
                    pairpots=[
                        coulomb('ElecElec','e','e'),
                        pseudopotential('PseudoPot','ion0','psi0',('B','B.pp.xml'),('C','C.pp.xml')),
                        coulomb('IonIon','ion0','ion0'),
                        ],
                    estimators = [
                        energydensity('edvoronoi','e','ion0','voronoi',-4,5),
                        energydensity('edcell','e','ion0',
                            spacegrid('cartesian',
                                origin = 'zero',
                                x = ('a1',.5,'-1 (192) 1'),
                                y = ('a2',.5,'-1 (1) 1'),
                                z = ('a3',.5,'-1 (1) 1')
                                )
                            )
                        ]
                    ),
                wavefunction = section('psi0','e',
                    determinantset = section('bspline','Si.pwscf.h5','ion0',
                        sort = 1,
                        tilematrix = array([[1,0,0],[0,1,0],[0,0,1]]),
                        twistnum = 0,
                        slaterdeterminant = [
                            determinant('updet',64,'ground',0),
                            determinant('downdet',63,'ground',1)
                            ],
                        jastrows = [
                            twobody('J2','bspline',
                                    ('u','u',3.9,[0,0,0,0,0,0]),
                                    ('u','d',3.9,[0,0,0,0,0,0])),
                            onebody('J1','bspline','ion0',
                                    ('C',3.9,[0,0,0,0,0,0]),                            
                                    ('B',3.9,[0,0,0,0,0,0]))
                            ]
                        )
                    )
                ),
            calculations=[
                loop(4,
                    cslinear(
                        blocks = 3125,
                        warmupsteps = 5,
                        steps = 2,
                        samples = 80000,
                        timestep = .5,
                        minmethod = 'rescale',
                        gevmethod = 'mixed',
                        exp0=-15,
                        nstabilizers = 5,
                        stabilizerscale = 3,
                        stepsize=.35,
                        alloweddifference=1e-5,
                        beta = .05,
                        bigchange = 5.,
                        energy = 0.,
                        unreweightedvariance = 0.,
                        reweightedvariance = 0.,
                        estimator = localenergy(hdf5='no') 
                        )
                    ),
                vmc(
                    blocks = 2,
                    steps = 500,
                    substeps = 3,
                    timestep = .5,
                    estimator = localenergy(hdf5='yes') 
                    ),
                dmc(
                    walkers = 72,
                    blocks = 2,
                    steps = 50,
                    timestep = .01,
                    nonlocalmove = 'yes',
                    estimator = localenergy(hdf5='no') 
                    )
                ]
            )
    #end if
#end if
