##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  sqd_input.py                                                      #
#    Supports input file I/O and manipulation for the SQD code.      #
#    Essentially an early fork of qmcpack_input.py.                  #
#                                                                    #
#  Content summary:                                                  #
#    SqdInput                                                        #
#      SimulationInput class for the SQD code.                       #
#                                                                    #
#    generate_sqd_input                                              #
#      User-facing function to generate SQD input files.             #
#                                                                    #
#    For descriptions of remaining contents, see the header of       #
#      qmcpack_input.py.                                             #
#                                                                    #
#                                                                    #
#====================================================================#


import os
import inspect
import keyword
from numpy import fromstring,empty,array,float64,\
    loadtxt,ndarray,dtype,sqrt,pi,arange,exp,eye,ceil,mod
from StringIO import StringIO
from superstring import string2val
from generic import obj
from xmlreader import XMLreader,XMLelement
from periodic_table import pt as periodic_table
from developer import DevBase
from structure import Structure
from physical_system import PhysicalSystem
from simulation import SimulationInput



yesno     = {True:'yes' ,False:'no'}
truefalse = {True:'true',False:'false'}
onezero   = {True:'1'   ,False:'0'}
boolmap={'yes':True,'no':False,'true':True,'false':False,'1':True,'0':False}

def SQDbool(b):
    return boolmap[b]
#end def SQDbool

def is_bool(var):
    return var in boolmap
#end def is_bool

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
    #if is_bool(attr):
    #    val = SQDbool(attr)
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
    



def is_list(s):
    try:
        t = eval(s)
        if isinstance(t,list):
            return True
        #end if
    except (ValueError,NameError):
        return False
#end def is_list

def string_to_val(s):
    if is_int(s):
        val = int(s)
    elif is_float(s):
        val = float(s)
    elif is_list(s):
        val = eval(s)
    else:
        val = s
    #end if
    return val
#end def string_to_val


class SQDobj(DevBase):
    None
#end class SQDobj


class meta(obj):
    None
#end class meta


class section(SQDobj):
    def __init__(self,*args,**kwargs):
        self.args   = args
        self.kwargs = kwargs
    #end def __init__
#end class section


class collection(SQDobj):
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
                if 'identifier_type' in element.__class__.__dict__:
                    identifier_type = element.identifier_type
                else:
                    identifier_type = str
                #end if
                if identifier_type==str:
                    for ident in identifier:
                        if ident in element:
                            key+=element[ident]
                        #end if
                    #end for
                elif identifier_type==tuple:
                    key = element.tuple(*identifier)
                else:
                    self.error('identifier_type '+str(identifier_type)+' has not yet been implemented')
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
    #end if
#end def make_collection


class classcollection(SQDobj):
    def __init__(self,*classes):
        if len(classes)==1 and isinstance(classes[0],list):
            classes = classes[0]
        #end if
        self.classes = classes
    #end def __init__
#end class classcollection


class Names(SQDobj):
    names = set([
            'atom','c','condition','eigensolve','grid','hamiltonian','id',
            'l','m','n','name','npts','num_closed_shells','orbital',
            'orbitalset','parameter','project','rf','ri','s','scale',
            'series','simulation','type'])

    bools = dict()

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
        if expanded in self.escaped_names:
            expanded = expanded[:-1]
        #end if
        if expanded in self.expanded_names:
            expanded = self.expanded_names[expanded]
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




class SQDxml(Names):

    def init_from_args(self,args):
        print
        print args
        print
        self.not_implemented()
    #end def init_from_args



    @classmethod
    def init_class(cls):
        cls.class_set_optional(
            tag        = cls.__name__,
            attributes = [],
            elements   = [],
            text       = None,
            parameters = [],
            attribs    = [],
            costs      = [],
            h5tags     = [],
            defaults   = obj()
            )
        for v in ['attributes','elements','parameters','attribs','costs','h5tags']:
            names = cls.class_get(v)
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
        indent  = indent_level*pad
        ip = indent+pad
        ipp= ip+pad
        c = indent+'<'+self.tag
        for a in self.attributes:
            if a in self:
                val = self[a]
                if isinstance(val,str):
                    val = self.expand_name(val)
                #end if
                c += ' '+self.expand_name(a)+'='
                if a in self.bools and (val==True or val==False):
                    c += '"'+self.bools[a][val]+'"'
                else:
                    c += '"'+param.write(val)+'"'
                #end if
            #end if
        #end for
        if first:
            None
            #c+=' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://www.mcc.uiuc.edu/qmc/schema/molecu.xsd"'
        #end if
        #no_contents = len(set(self.keys())-set(self.elements)-set(self.plurals.keys()))==0
        no_contents = len(set(self.keys())-set(self.attributes))==0
        if no_contents:
            c += '/>\n'
        else:
            c += '>\n'
            for v in self.h5tags:
                if v in self:
                    c += param.write(self[v],name=self.expand_name(v),tag='h5tag',mode='elem',pad=ip)
                #end if
            #end for
            for v in self.costs:
                if v in self:
                    c += param.write(self[v],name=self.expand_name(v),tag='cost',mode='elem',pad=ip)
                #end if
            #end for
            for p in self.parameters:
                if p in self:
                    c += param.write(self[p],name=self.expand_name(p),mode='elem',pad=ip)
                #end if
            #end for
            for a in self.attribs:
                if a in self:
                    c += param.write(self[a],name=self.expand_name(a),tag='attrib',mode='elem',pad=ip)
                #end if
            #end for
            for e in self.elements:
                if e in self:
                    elem = self[e]
                    if isinstance(elem,SQDxml):
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
                        keys = coll.keys()
                        keys.sort()
                        for key in keys:
                            instance = coll[key]
                            c += instance.write(indent_level+1)
                        #end for
                    #end if
                #end if
            #end for
            #for p in self.plurals.keys():
            #    if p in self:
            #        for instance in self[p]:
            #            c += instance.write(indent_level+1)
            #        #end for
            #    #end if
            ##end for
            if self.text!=None:
                #c+=ip+param.write(self[self.text])+'\n'
                c+=param.write(self[self.text],mode='elem',pad=ip,tag=None)
            #end if
            c+=indent+'</'+self.tag+'>\n'
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
            if a in self.bools:
                self[a] = boolmap[xml._attributes[al[a]]]
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
            elif inspect.isclass(value):
                if elements:
                    defval = value()
                #end if
            elif inspect.isfunction(value):
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
                if isinstance(value,SQDxml):
                    value.incorporate_defaults(elements,overwrite)
                elif isinstance(value,collection):
                    for v in value:
                        if isinstance(v,SQDxml):
                            v.incorporate_defaults(elements,overwrite)
                        #end if
                    #end for
                #end if
            #end for
        #end if
    #end def incorporate_defaults                    


    def check_junk(self,junk):
        if len(junk)>0:
            msg = self.tag+' does not have the following attributes/elements:\n'
            for jname in junk:
                msg+='    '+jname+'\n'
            #end for
            self.error(msg,'SqdInput',exit=False,trace=False)
            #self.error(msg,'SqdInput')
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
                not_xml  = not isinstance(namedict[name],SQDxml)
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
            if isinstance(value,SQDxml):
                value.get(names,namedict,host,root=False)
            elif isinstance(value,collection):
                for n,v in value.iteritems():
                    name_absent = not n in namedict 
                    not_element = False
                    if not name_absent:
                        not_xml  = not isinstance(namedict[n],SQDxml)
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
                    if isinstance(v,SQDxml):
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
                if isinstance(val,SQDxml) or isinstance(val,collection):
                    remove.append(rname)
                #end if
            #end if
        #end for
        for name in remove:
            del self[name]
        #end for
        for name,value in self.iteritems():
            if isinstance(value,SQDxml):
                value.remove(*names)
            elif isinstance(value,collection):
                for element in value:
                    if isinstance(element,SQDxml):
                        element.remove(*names)
                    #end if
                #end for
            #end if
        #end for
    #end def remove


    def replace(self,*args,**kwargs):
        if len(args)==2 and isinstance(args[0],str) and isinstance(args[1],str):
            vold,vnew = args
            args = [(vold,vnew)]
        #end for
        for valpair in args:
            vold,vnew = valpair
            for var,val in self.iteritems():
                not_coll = not isinstance(val,collection)
                not_xml  = not isinstance(val,SQDxml)
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
                    not_xml  = not isinstance(val,SQDxml)
                    not_arr  = not isinstance(val,ndarray)
                    if not_coll and not_xml and not_arr and val==vold:
                        self[var] = vnew
                    #end if
                #end if
            #end if
        #end for
        for vname,val in self.iteritems():
            if isinstance(val,SQDxml):
                val.replace(*args,**kwargs)
            elif isinstance(val,collection):
                for v in val:
                    if isinstance(v,SQDxml):
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
            single = isinstance(element,SQDxml)
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
            if isinstance(value,SQDxml):
                value.pluralize()
            elif isinstance(value,collection):
                for v in value:
                    if isinstance(v,SQDxml):
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
                is_qxml1 = isinstance(value1,SQDxml)
                is_qxml2 = isinstance(value2,SQDxml)
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
                        if isinstance(v1,SQDxml) and isinstance(v2,SQDxml):
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
            if isinstance(value,SQDxml):
                value.remove_empty()
                if len(value)==0:
                    del self[name]
                #end if
            elif isinstance(value,collection):
                ns = list(value.keys())
                for n in ns:
                    v = value[n]
                    if isinstance(v,SQDxml):
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
#end class SQDxml



class SQDxmlFactory(Names):
    def __init__(self,name,types,typekey='',typeindex=-1,typekey2='',default=None):
        self.name = name
        self.types = types
        self.typekey = typekey
        self.typekey2 = typekey2
        self.typeindex = typeindex
        self.default = default
    #end def __init__

    def __call__(self,*args,**kwargs):
        #emulate SQDxml.__init__
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
        None # this is for compatibility with SQDxml only (do not overwrite)
    #end def init_class
#end class SQDxmlFactory



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


    def write(self,value,mode='attr',tag='parameter',name=None,pad='   '):
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
                        vfmt = ':16.8e'
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
                cname = self.condense_name(name)
                if cname in self.bools and (value==0 or value==1):
                    val = self.bools[cname][value]
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




class simulation(SQDxml):
    elements   = ['project','atom','eigensolve']
#end class simulation


class project(SQDxml):
    attributes = ['id','series']
#end class project

class atom(SQDxml):
    attributes = ['name','num_closed_shells']
    elements   = ['grid','orbitalset','hamiltonian']
    #identifier = 'name'
#end class atom

class grid(SQDxml):
    attributes = ['type','ri','rf','npts','scale']
#end class grid

class orbitalset(SQDxml):
    attributes = ['condition']
    elements   = ['orbital']
#end class orbitalset

class orbital(SQDxml):
    attributes = ['n','l','m','s','c']
    identifier = 'n','l','m','s'
    identifier_type = tuple
#end class orbital

class hamiltonian(SQDxml):
    attributes = ['type']
    parameters = ['z']
#end class hamiltonian

class eigensolve(SQDxml):
    parameters = ['max_iter','etot_tol','eig_tol','mix_ratio']
#end class eigensolve




classes = [   #standard classes
    simulation,project,atom,grid,orbitalset,orbital,hamiltonian,eigensolve
    ]
types = dict( #simple types and factories
    )
plurals = obj(
    orbitals      = 'orbital'
    )
plurals_inv = plurals.inverse()
plural_names = set(plurals.keys())
single_names = set(plurals.values())
Names.set_expanded_names(
    )
for c in classes:
    c.init_class()
    types[c.__name__] = c
#end for


#set default values
simulation.defaults.set(
    project      = project,
    atom         = atom,
    eigensolve   = eigensolve
    )
project.defaults.set(
    series = 0
    )
atom.defaults.set(
    grid        = grid,
    orbitalset  = orbitalset,
    hamiltonian = hamiltonian
    )
grid.defaults.set(
    type='log',ri=1e-6,rf=400,npts=2001
    )
orbitalset.defaults.set(
    condition = 'spin_space'
    )
hamiltonian.defaults.set(
    type = 'nuclear'
    )
eigensolve.defaults.set(
    max_iter  = 1000,
    etot_tol  = 1e-8,
    eig_tol   = 1e-12,
    mix_ratio = .3
    )






class SqdInput(SimulationInput,Names):
    
    opt_methods = set([])

    simulation_type = simulation

    default_metadata = meta(
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
        elif isinstance(arg0,SQDxml) and arg1==None:
            element = arg0
        elif isinstance(arg0,meta) and isinstance(arg1,SQDxml):
            metadata = arg0
            element  = arg1
        else:
            self.error('input arguments of types '+arg0.__class__.__name__+' and '+arg0.__class__.__name__+' cannot be used to initialize SqdInput')
        #end if
        if metadata!=None:
            self._metadata = metadata
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
        if '_metadata' in self:
            elem_names.remove('_metadata')
        #end if
        if len(elem_names)>1:
            self.error('sqd input cannot have more than one base element\n  You have provided '+str(len(elem_names))+': '+str(elem_names))
        #end if
        return self[elem_names[0]]
    #end def get_base

    def get_basename(self):
        elem_names = list(self.keys())
        if '_metadata' in self:
            elem_names.remove('_metadata')
        #end if
        if len(elem_names)>1:
            self.error('sqd input cannot have more than one base element\n  You have provided '+str(len(elem_names))+': '+str(elem_names))
        #end if
        return elem_names[0]
    #end def get_basename

    def read(self,filepath):
        if os.path.exists(filepath):
            element_joins=[]
            element_aliases=dict()
            xml = XMLreader(filepath,element_joins,element_aliases,warn=False).obj
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
                    if isinstance(elem,SQDxml):
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
    #end def read


    def write_text(self,filepath=None):
        c = ''
        header = '''<?xml version="1.0"?>
'''
        c+= header
        if '_metadata' in self:
            Param.metadata = self._metadata
        elif Param.metadata == None:
            Param.metadata = self.default_metadata
        #end if
        base = self.get_base()
        c+=base.write(first=True)
        Param.metadata = None
        return c
    #end def write_text

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

    def read_xml(self,filepath):
        if os.path.exists(filepath):
            element_joins=['qmcsystem']
            element_aliases=dict(loop='qmc')
            xml = XMLreader(filepath,element_joins,element_aliases,warn=False).obj
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
                    self.error('host xml section for '+qname+' not found','SqdInput')
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


    def get_output_info(self,list=True):
        project = self.simulation.project
        prefix = project.id+'.s'+str(project.series).zfill(3)+'.'

        outfiles = obj(
            exchange = prefix+'exchange',
            h5       = prefix+'h5',
            hartree  = prefix+'hartree',
            log      = prefix+'log',
            orb      = prefix+'orb.dat',
            qmc      = prefix+'qmc.xml',
            Vext     = prefix+'Vext.xml'
            )
        if list:
            return outfiles.list()
        else:
            return outfiles
        #end if
    #end def get_output_info



    def incorporate_system(self,system):
        element = system.structure.elem[0]
        Z = periodic_table[element].atomic_number

        atom = self.simulation.atom
        atom.name = element
        atom.grid.scale = Z
        atom.hamiltonian.z = Z
    #end def incorporate_system
        

    def return_system(self):
        system = PhysicalSystem(
            structure = Structure(
                elem = [self.simulation.atom.name],
                pos  = [[0,0,0]]
                )
            )
        return system
    #end def return_system
#end class SqdInput







class Shells(DevBase):
    channel_names = tuple('spdfghik')
    channels = obj()
    for l in range(len(channel_names)):
        channels[channel_names[l]] = range(l,-l-1,-1)
    max_shells = 7
    all_shells = obj()
    for n in range(1,max_shells+1):
        shell = obj()
        for l in range(n):
            shell[channel_names[l]] = range(l,-l-1,-1)
        #end for
        all_shells[n] = shell
    #end for

    channel_indices = obj()
    for i in range(len(channel_names)):
        channel_indices[channel_names[i]] = i
    #end for

    core_names = ['He','Ne','Ar','Kr','Xe','Rn']
    core_orbitals = obj(
        He = '1s',
        Ne = '1s 2s 2p',
        Ar = '1s 2s 2p 3s 3p',
        Kr = '1s 2s 2p 3s 3p 4s 3d 4p',
        Xe = '1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p',
        Rn = '1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p'
        )
    ncore_shells = obj(
        He = 1,
        Ne = 2,
        Ar = 3,
        Kr = 4,
        Xe = 5,
        Rn = 6
        )
    for cname in core_orbitals.keys():
        co = core_orbitals[cname].split()
        olist = []
        for o in co:
            n  = int(o[0])
            ls = o[1]
            olist.append((n,ls))
        #end for
        core_orbitals[cname] = olist
    #end for
    cores = obj()
    for cname,orblist in core_orbitals.iteritems():
        core = obj()
        cores[cname] = core
        for (n,ls) in orblist:
            if not n in core:
                core[n] = obj()
            #end if
            shell = core[n]
            if not ls in shell:
                shell[ls] = list(all_shells[n][ls])
            #end if
        #end for
    #end for

    orbital_fill_order_list = '1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p'.split()
    orbital_fill_order = []
    for orbital in orbital_fill_order_list:
        n  = int(orbital[0])
        ls = orbital[1]
        l  = channel_indices[ls]
        mult = len(channels[ls])
        mlist = list(channels[ls])
        orbital_fill_order.append(obj(n=n,l=l,ls=ls,mult=mult,mlist=mlist,total_mult=2*mult))
    #end for
    #m_l fill order is max m_l to min m_l, i.e. l,l-1,l-2,...,-l+2,-l+1,-l


    @classmethod
    def hunds_rule_filling(cls,atom,net_charge=0,net_spin='ground',location='hunds_rule_filling'):
        if isinstance(atom,str) and atom in periodic_table:
            Z = periodic_table[atom].atomic_number
        elif isinstance(atom,int) and atom>0 and atom<110:
            Z = atom
        else:
            cls.class_error('expected atomic symbol or atomic number for atom\n  you provided '+str(atom),location)
        #end if

        nelectrons = Z-net_charge

        if isinstance(net_spin,int):
            nup   = float(nelectrons + net_spin)/2
            ndown = float(nelectrons - net_spin)/2        
            if abs(nup-int(nup))>1e-3:
                cls.class_error('requested spin state {0} incompatible with {1} electrons'.format(net_spin,nelectrons),location)
            #end if
            nup   = int(nup)
            ndown = int(ndown)
        elif net_spin=='ground' or net_spin==None:
            net_spin = None
            nup      = None
            ndown    = None
        else:
            cls.class_error("net_spin must be 'ground'/None or integer\n  you provided "+str(net_spin))
        #end if

        if net_spin!=None and nup+ndown!=nelectrons:
            cls.class_error('number of up and down electrons does add up to the total\n  this may reflect an error in your input or the code\n  please check: nel={0} nup={1} ndown={2} Z={3} net_charge={4} net_spin={5}'.format(nelectrons,nup,ndown,Z,net_charge,net_spin),location)
        #end if

        closed_orbitals = []
        open_orbital   = None
        nel = 0
        nud_closed = 0
        for orbital in cls.orbital_fill_order:
            nel_orb = orbital.total_mult
            if nel<nelectrons:
                if nel+nel_orb<nelectrons:
                    closed_orbitals.append(orbital)
                    nud_closed += orbital.mult
                else:
                    open_orbital = orbital
                    break
                #end if
            #end if
            nel+=nel_orb
        #end if


        up = ''
        down = ''
        updown = ''

        o = open_orbital
        if net_spin!=None:
            nup_open   = nup   - nud_closed
            ndown_open = ndown - nud_closed

            if nup_open>open_orbital.mult or ndown_open>open_orbital.mult:
                cls.class_error('more up or down electrons in open shell than will fit\n  open_shell={0}{1}, up/down_size={2}, nup={3}, ndown={4}'.format(o.n,o.ls,o.mult,nup_open,ndown_open),'Developer')
            #end if
        else:
            nopen = nelectrons - 2*nud_closed
            nup_open = min(nopen,o.mult)
            ndown_open = nopen - nup_open
        #end if

        for orbital in closed_orbitals:
            updown += '{0}{1}'.format(orbital.n,orbital.ls)
        #end for
        if nup_open>0:
            up   = '{0}{1}{2}'.format(o.n,o.ls,str(tuple(o.mlist[:nup_open]))).replace(',)',')')
        #end if
        if ndown_open>0:
            down = '{0}{1}{2}'.format(o.n,o.ls,str(tuple(o.mlist[:ndown_open]))).replace(',)',')')
        #end if
        
        if up=='':
            up=None
        #end if
        if down=='':
            down=None
        #end if
        if updown=='':
            updown=None
        #end if

        net_spin_found = nup_open-ndown_open
        if net_spin!=None and net_spin_found!=net_spin:
            cls.error('spin state determined incorrectly\n  net_spin requested: {0}\n  net_spin found: {1}\n  nup,ndown: {2},{3}'.format(net_spin,net_spin_found,nup_open,ndown_open))
        #end if
            
        return up,down,updown,net_spin_found
    #end def hunds_rule_filling



    

    def __init__(self,shells=None,location='Shells'):
        self.location = location
        self.shells   = obj()
        self.core     = None
        self.nclosed  = 0
        if shells is None:
            return
        #end if
        if isinstance(shells,str):
            self.read_shell_string(shells)
        elif isinstance(shells,obj) or isinstance(shells,dict):
            if set(shells.keys()) <= set(self.all_shells.keys()):
                self.shells.transfer_from(shells,copy=True)
            else:
                self.error('unexpected key values for shells\n  expected values: '+str(shells.keys())+'\n  you provided '+str(self.all_shells.keys()))
            #end if
        else:
            self.error('must provide a string, dict, or obj describing atomic shells\n  you provided '+str(shells))
        #end if
        self.check_shells()
    #end def __init__


    def read_shell_string(self,ss):
        self.shell_string = str(ss)
        ss = ss.replace('[',' ').replace(']',' ').replace('(',' [').replace(')','] ')
        ss = ss.replace(',   ',',').replace(',  ',',').replace(', ',',')
        ss = ss.replace('   ,',',').replace('  ,',',').replace(' ,',',')
        ss = ss.replace('[   ','[').replace('[  ','[').replace('[ ','[')
        ss = ss.replace('   ]',']').replace('  ]',']').replace(' ]',']')
        for core_name in self.core_names:
            ss = ss.replace(core_name,' '+core_name+' ')
        #end for
        for channel in self.channel_names:
            ss = ss.replace(channel,' '+channel+' ')
        #end for

        shells = self.shells
        sl = ss.split()
        if len(sl)>0:
            i = 0
            if sl[i] in self.cores:
                shells.transfer_from(self.cores[sl[i]],copy=True)
                i+=1
            #end if
            n=None
            l=None
            m=None
            shell = None
            channel = None
            while i<len(sl):
                v = string_to_val(sl[i])
                if isinstance(v,int):
                    if n!=None:
                        if n in shells:
                            shells[n].transfer_from(shell)
                        else:
                            shells[n] = shell
                        #end if
                        n=None
                    #end if
                    n = v
                    if n>self.max_shells:
                        self.error('maximum shell number is {0}\n you requested {1}'.format(self.max_shells,n))
                    #end if
                    shell = obj()
                elif isinstance(v,str):
                    l = v.lower()
                    if not l in self.channel_names:
                        self.error('you requested an invalid channel: '+str(l)+'\n allowed channels are '+str(self.channel_names))
                    #end if
                    if i+1<len(sl):
                        v1 = string_to_val(sl[i+1])
                        if isinstance(v1,list):
                            m=v1
                            ma = abs(array(m))
                            if ma.max()>self.channel_indices[l]:
                                self.error('maximum |m| for {0} channel is {1}\n  you requested {2}: {3}'.format(l,self.channel_indices[l],ma.max(),m))
                            #end if
                            channel = m
                            m=None
                            i+=1
                        else:
                            ln = self.channel_indices[l]
                            channel = range(-ln,ln+1)
                        #end if
                    else:
                        ln = self.channel_indices[l]
                        channel = range(-ln,ln+1)
                    #end if
                    if n!=None:
                        shell[l] = channel
                        l=None
                    #end if
                #end if
                i+=1
            #end while
            if n!=None:
                if n in shells:
                    shells[n].transfer_from(shell)
                else:
                    shells[n] = shell
                #end if
                n=None
            #end if
        #end if
    #end def read_shell_string


    def check_shells(self):
        ref = self.all_shells
        shells = self.shells
        rkn = set(ref.keys())
        skn = set(shells.keys())
        errors = False
        if not skn <= rkn:
            self.error('shell indices (n) are invalid\n  options for valid shell indices: '+str(list(rkn))+'\n  shell indices of self: '+str(list(skn)),exit=False)
            errors = True
        #end if
        for n,shell in shells.iteritems():
            rshell = ref[n]
            rkl = set(rshell.keys())
            skl = set(shell.keys())
            if not skl<=rkl:
                self.error('channel keys (l) are invalid\n options for valid channel keys: '+str(list(rkl))+'\n  channel keys of self: '+str(list(skl)),exit=False)
                errors = True
            #end if
            for l,channel in shell.iteritems():
                rchannel = rshell[l]
                rkm = set(rchannel)
                skm = set(channel)
                if not skm<=rkm:
                    self.error('azimuthal indices (m) are invalid\n  options for valid azimuthal indices: '+str(list(rkm))+'\n  azimuthal indices of self: '+str(list(skm)),exit=False)
                    errors = True
                #end if
            #end for
        #end for
        if errors:
            self.log('\nreference shells:\n'+str(ref))
            self.log('\nself shells:\n'+str(shells))
            self.error('encountered errors')
        #end if
    #end def check_shells


    def partition(self):
        shells = self.shells
        #find all closed subshells
        closed = []
        for n,shell in shells.iteritems():
            for l,channel in shell.iteritems():
                if set(channel)==set(self.all_shells[n][l]):
                    closed.append((n,l))
                #end if
            #end for
        #end for
        closed = set(closed)
        #find what the core is, He,Ne,Ar, etc
        core_orb = set()
        for i in range(len(self.core_names)-1,-1,-1):
            core_name = self.core_names[i]
            core_orbitals = set(self.core_orbitals[core_name])
            if core_orbitals<=closed:
                self.core = core_name
                self.nclosed = self.ncore_shells[core_name]
                core_orb = core_orbitals
                break
            #end if
        #end for
        #remove orbitals belonging to the core
        for (n,l) in core_orb:
            del shells[n][l]
        #end for
        for n in list(shells.keys()):
            if len(shells[n])==0:
                del shells[n]
            #end if
        #end for
    #end def partition
                

    def orbitals(self,spin):
        if spin!='up' and spin!='down' and spin!=1 and spin!=-1:
            self.error('spin must be up/1 or down/-1\n  you provided: '+str(spin))
        elif spin=='up':
            s = 1
        elif spin=='down':
            s = -1
        #end if
        
        orbitals = []
        for n,shell in self.shells.iteritems():
            for lname,channel in shell.iteritems():
                l = self.channel_indices[lname]
                for m in channel:
                    orbitals.append(orbital(n=n,l=l,m=m,s=s,c=1.0))
                #end for
            #end for
        #end for
        return orbitals
    #end def orbitals
#end class Shells
hunds_rule_filling = Shells.hunds_rule_filling



def generate_orbitalset(up=None,down=None,updown=None,location='generate_orbitalset'):
    nclosed  = 0
    uorb   = []
    dorb   = []
    if isinstance(updown,str):
        bshells = Shells(updown,location)
        bshells.partition()
        nclosed = bshells.nclosed
        uorb = bshells.orbitals('up')
        dorb = bshells.orbitals('down')
    #end if
    ushells = Shells(up,location)
    dshells = Shells(down,location)

    uorb += ushells.orbitals('up')
    dorb += dshells.orbitals('down')

    orbset = orbitalset(
        orbitals = make_collection(uorb+dorb)
        )

    return orbset,nclosed
#end def generate_orbitalset




def generate_sqd_input(id         = None,
                       series     = 0,
                       system     = None,
                       filling    = None,
                       net_spin   = 'none',
                       up         = None,
                       down       = None,
                       updown     = None,
                       grid_type  = 'log',
                       ri         = 1e-6,
                       rf         = 400,
                       npts       = 10001,
                       max_iter   = 1000,
                       etot_tol   = 1e-8,
                       eig_tol    = 1e-12,
                       mix_ratio  = 0.7  ):

    location          = 'generate_sqd_input'


    metadata = SqdInput.default_metadata.copy()
    si = SqdInput(
        metadata,
        simulation(
            project = section(
                series = series
                ),
            atom = section(
                grid = section(type=grid_type,ri=ri,rf=rf,npts=npts),
                orbitalset = section(),
                hamiltonian = section()
                ),
            eigensolve = section(
                max_iter  = max_iter,
                etot_tol  = etot_tol,
                eig_tol   = eig_tol,
                mix_ratio = mix_ratio
                )
            )
        )

    #set the atomic system
    if system==None:
        SqdInput.class_error('system (atom) must be provided',location)
    elif isinstance(system,str) and system in periodic_table:
        Z = periodic_table[system].atomic_number
        net_charge = Z
        net_spin   = 0
        for orbital in orbset.orbitals:
            net_charge -= 1
            net_spin   += orbital.s
        #end for
        system = PhysicalSystem(
            structure  = Structure(elem=[system],pos=[[0,0,0]]),
            net_charge = net_charge,
            net_spin   = net_spin
            )
        si.incorporate_system(system)
    elif isinstance(system,PhysicalSystem):
        si.incorporate_system(system)
    else:
        SqdInput.class_error('system must be an atomic symbol or a PhysicalSystem object','generate_sqd_input')
    #end if

    #generate spin state via Hund's rule if requested
    if filling is None and up is None and down is None and updown is None:
        SqdInput.class_error('filling or up/down/updown must be provided',location)
    #end if
    if filling !=None and not isinstance(filling,str):
        SqdInput.class_error('filling must be a string\n  you provided '+str(filling),location)
    elif isinstance(filling,str) and filling.lower()=='hund':
        atom = system.structure.elem[0]
        net_charge = system.net_charge
        if net_spin is 'none':
            net_spin = system.net_spin
        #end if
        up,down,updown,net_spin = hunds_rule_filling(
            atom       = atom,
            net_charge = net_charge,
            net_spin   = net_spin,
            location   = location
            ) 
    elif filling!=None:
        SqdInput.class_error("{0} is not a valid choice for filling\n valid options are: 'hund'".format(filling),location)
    #end if

    #generate the orbitalset
    orbset,nclosed = generate_orbitalset(
        up     = up,
        down   = down,
        updown = updown
        )


    if id is None:
        id = system.structure.elem[0]
    #end if
    sim = si.simulation
    sim.project.id = id
    sim.atom.set(
        num_closed_shells = nclosed,
        orbitalset = orbset
        )

    si.incorporate_defaults(elements=False,overwrite=False,propagate=True)

    return si
#end def generate_sqd_input








