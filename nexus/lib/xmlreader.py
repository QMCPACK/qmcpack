##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  xmlreader.py                                                      #
#    Provides support for reading general XML files and converting   #
#    contents to numeric form, where appropriate.  Based on the      #
#    expat XML parser.                                               #
#                                                                    #
#  Content summary:                                                  #
#    XMLreader                                                       #
#      Class reads an XML file, converting it into a structured      #
#      object format.                                                #
#                                                                    #
#    readxml                                                         #
#      Function interface to XMLreader                               #
#                                                                    #
#    XMLelement                                                      #
#      Class represents a single XML element                         #
#                                                                    #
#====================================================================#


from xml.parsers import expat
from numpy import array
import sys
import keyword
import re
import os
from inspect import getmembers

from superstring import \
    find_matching_pair, \
    remove_pair_sections, \
    remove_empty_lines, \
    valid_variable_name,\
    string2val

from generic import obj
from developer import DevBase


class XMLelement(DevBase):
    def _escape_name(self,name):
        if name in self._escape_names:
            name=name+'_'
        #end if
        return name
    #end def escape_name

    def _set_parent(self,parent):
        self._parent=parent
        return
    #end def set_parent

    def _add_xmlattribute(self,name,attribute):
        self._attributes[name]=attribute
        return 
    #end def add_attribute

    def _add_element(self,name,element):
        element._name=name
        self._elements[name]=element
        return 
    #end def add_element

    def _add_text(self,name,text):
        self._texts[name]=text
        return 
    #end def add_text

    def _to_string(self):
        s=''
        if len(self._attributes)>0:
            s+='  attributes:\n'
            for k,v in self._attributes.iteritems():
                s+= '    '+k+' = '+str(v)+'\n'
            #end for
        #end if
        if len(self._elements)>0:
            s+= '  elements:\n'
            for k,v in self._elements.iteritems():
                s+= '    '+k+'\n'
            #end for
        #end if
        if len(self._texts)>0:
            s+= '  texts:\n'
            for k,v in self._texts.iteritems():
                s+= '    '+k+'\n'
            #end for
        #end if
        return s
    #end def list

#    def __str__(self):
#        return self._to_string()
#    #end def __str__
#
#    def __repr__(self):
#        return self._to_string()
#    #end def __repr__

    def __init__(self):
        self._name=''
        self._parent=None
        self._elements=obj()
        self._texts=obj()
        self._attributes=obj()
        self._element_counts=obj()
        self._ntexts=0

        self._escape_names=None
        #self._escape_names=set(dict(getmembers(self)).keys()) | set(keyword.kwlist)
        self._escape_names=set(keyword.kwlist)
        return
    #end def __init__


    def condense(self):
        for name,elem in self._elements.iteritems():
            if isinstance(elem,XMLelement):
                elem.condense()
            #end if
        #end if
        cnames = []
        for name in self._elements.keys():
            if name[-1]=='1' and not name[-2].isdigit():
                cnames.append(name[:-1])
            #end if
        #end if
        for cname in cnames:
            cmax = 1
            for name,elem in self._elements.iteritems():
                ns = name.split(cname)
                if len(ns)==2 and ns[1].isdigit():
                    cmax = max(cmax,int(ns[1]))
                #end if
            #end if
            names = set()
            for n in range(1,cmax+1):
                name = cname+str(n)
                names.add(name)
            #end for
            not_present = names-set(self._elements.keys())
            if len(not_present)==0:
                collection = []
                for n in range(1,cmax+1):
                    name = cname+str(n)            
                    collection.append(self._elements[name])
                    del self._elements[name]
                    del self[name]
                #end for
                self._elements[cname] = collection
                self[cname] = collection
            #end if
        #end for
    #end def condense


    def convert_numeric(self):
        for name,attr in self._attributes.iteritems():
            self[name] = string2val(attr)
        #end for
        if 'text' in self:
            self.value = string2val(self.text)
            del self.text
        #end if
        texts = []
        for name,elem in self._elements.iteritems():
            if isinstance(elem,XMLelement):
                if 'text' in elem and len(elem._attributes)==0 and len(elem._elements)==0:
                    self[name] = string2val(elem.text)
                    texts.append(name)
                else:
                    elem.convert_numeric()
                #end if
            #end if
        #end if
        for name in texts:
            self._elements[name] = self[name]
        #end for
    #end def convert_numeric

                    
    def remove_hidden(self):
        for name,elem in self._elements.iteritems():
            if isinstance(elem,XMLelement):
                elem.remove_hidden()
            elif isinstance(elem,list):
                for e in elem:
                    if isinstance(e,XMLelement):
                        e.remove_hidden()
                    #end if
                #end for
            #end if
        #end for
        remove = []
        for name,value in self.iteritems():
            if str(name)[0]=='_':
                remove.append(name)
            #end if
        #end for
        for name in remove:
            del self[name]
        #end for
    #end def remove_hidden
#end class XMLelement



'''
  class XMLReader
    reads an xml file and creates a dynamic object out of its contents
'''
class XMLreader(DevBase):
    def __init__(self,fpath=None,element_joins=None,element_aliases=None,contract_names=False,strip_prefix=None,warn=True,xml=None):
        if element_joins is None:
            element_joins = []
        if element_aliases is None:
            element_aliases = {}

        #assign values
        self.fpath=fpath
        if fpath is None:
            self.base_path = None
        else:
            self.base_path = os.path.split(fpath)[0]
        #end if
        self.element_joins = set(element_joins)
        self.element_aliases = element_aliases
        self.contract_names = contract_names
        self.strip_prefix = strip_prefix
        self.warn = warn
        
        #create the parser
        self.parser = expat.ParserCreate()
        self.parser.StartElementHandler  = self.found_element_start
        self.parser.EndElementHandler    = self.found_element_end
        self.parser.CharacterDataHandler = self.found_text
        self.parser.AttlistDeclHandler   = self.found_attribute
        self.parser.returns_unicode = 0

        #read in xml file
        if xml is None:
            fobj = open(fpath,'r')
            self.xml = fobj.read()
        else:
            self.xml = xml
        #end if
        #remove all comments
        pair='<!--','-->'
        self.xml = remove_pair_sections(self.xml,pair)
        #process includes
        while self.xml.find('<include')!=-1:
            self.include_files()
            self.xml = remove_pair_sections(self.xml,pair)
        #end while
        #remove empty lines
        self.xml = remove_empty_lines(self.xml)
        #print self.xml

        #parse the xml and build the dynamic object
        self.nlevels=1
        self.ilevel=0
        self.pad=''
        #  Set the current xml element
        self.obj = XMLelement()
        self.cur=[self.obj]

        self.parser.Parse(self.xml,True)

        #the expat parser is troublesome in that it
        # -does not have typical class members
        # -is unpickleable
        # therefore it is removed after the dynamic object is built
        del self.parser

        return
    #end def __init__

    def include_files(self):
        pair = '<include','/>'
        qpair = '<?','?>'
        ir=0
        while ir!=-1:
            il,ir = find_matching_pair(self.xml,pair,ir)
            if ir!=-1:
                cont = self.xml[il:ir].strip(pair[0]).rstrip(pair[1])
                fname = cont.split('=',1)[1].strip().strip('"')
                fobj = open(os.path.join(self.base_path,fname),'r')
                fcont = fobj.read()
                fcont = remove_pair_sections(fcont,qpair)
                fobj.close()
                self.xml = self.xml.replace(self.xml[il:ir],fcont)
            #end if
        #end while
        return
    #end def include_files

    def increment_level(self):
        self.ilevel+=1
        self.nlevels = max(self.ilevel+1,self.nlevels)
        if self.ilevel+1==self.nlevels:
            self.cur.append(None)
        #end if
        self.pad = self.ilevel*'  '
        return
    #end def increment_level

    def decrement_level(self):
        self.ilevel-=1
        self.pad = self.ilevel*'  '
        return
    #end def decrement_level

    def found_element_start(self,ename,attributes):
        #print self.pad,name,attributes
        cur = self.cur[self.ilevel]
        if ename in self.element_aliases.keys():
            if self.element_aliases[ename].find('attributes')!=-1:
                exec 'name = '+self.element_aliases[ename]
            else:
                name = self.element_aliases[ename]
            #end if
        else:
            name=ename
        #end if

        #alter the name if it is a python keyword
        name = cur._escape_name(name)

        if self.contract_names:
            name = name.lower().replace('-','_')
        #end if
        if self.strip_prefix!=None:
            if name.startswith(self.strip_prefix):
                name = name.split(self.strip_prefix)[1]
            #end if
        #end if

        # joinable = in joins and no attributes
        # if in elements and joinable: don't add 
        # else if not in elements and joinable: add unnumbered
        # else if not in elements: add unnumbered
        # else: add numbered, if number==1: rename first element 
        joinable = name in self.element_joins and len(attributes.keys())==0
        epattern = re.compile(name+'\d+')
        in_elements=False
        for k in cur._elements.keys():
            if epattern.match(k) or k==name:
                in_elements=True
            #end if
        #end for
        #in_elements = name in cur._elements.keys()
        if in_elements and joinable:
            #check if there is a previous unjoinable element w/ same name
            if len(cur._elements[name]._attributes)!=0:
                #rename the prior element as a numbered one
                nelements=cur._element_counts[name]
                if nelements==1:
                    #it should be the first one
                    newname = name+str(1)
                    cur[newname]=cur[name]
                    cur._add_element(newname,cur[newname])
                    del  cur._elements[name]
                    del cur[name]
                else:
                    print 'prior unjoinable element is not the first'
                    print '  this should be impossible, stopping'
                    sys.exit()
                #end if
                #add the joinable element as unnumbered
                # later joinable elements will be joined to this one
                cur[name] = XMLelement()
                cur._add_element(name,cur[name])
            #end if
        elif not in_elements:
            #add unnumbered
            cur[name] = XMLelement()
            cur._add_element(name,cur[name])
            cur._element_counts[name]=1
        else:
            #add in a numbered way
            nelements=cur._element_counts[name]
            if nelements==1:
                #rename the first element
                newname = name+str(1)
                cur[newname]=cur[name]
                cur._add_element(newname,cur[newname])
                del  cur._elements[name]
                del cur[name]
            #end if
            nelements+=1
            newname = name + str(nelements)
            cur[newname] = XMLelement()
            cur._add_element(newname,cur[newname])
            cur._element_counts[name]=nelements
            name = newname
        #end if
        cur._elements[name]._parent = cur  #mark change
        self.increment_level()
        self.cur[self.ilevel] = cur._elements[name]
        cur = self.cur[self.ilevel]
        for kraw,v in attributes.iteritems():
            if self.contract_names:
                k = kraw.lower().replace('-','_')
            else:
                k = kraw
            #end if
            if valid_variable_name(k):
                kname = cur._escape_name(k)
                cur[kname] = v
                cur._add_xmlattribute(kname,cur[kname])
            else:
                if self.warn:
                    print 'xmlreader warning: attribute '+k+' is not a valid variable name and has been ignored'
                #end if
            #end if
        #end for
        return
    #end def found_element_start

    def found_element_end(self,name):
        self.cur[self.ilevel]=None
        self.decrement_level()
        #print self.pad,'end',name
        return
    #end def found_element_end

    def found_text(self,rawtext):
        text = rawtext.strip()
        if text!='':
            #print self.pad,text
            cur = self.cur[self.ilevel]
            if cur._ntexts>0:
                cur.text+='\n'+text
            else:
                cur.text = text
                cur._add_text('text',cur.text)
                cur._ntexts+=1
            #end if
        #end if
        return
    #end def found_text

    def found_attribute(self,ename,aname,atype,default,required):
        return
    #end def found_attribute
#end class XMLreader



def readxml(fpath=None,element_joins=None,element_aliases=None,contract_names=False,strip_prefix=None,warn=True,xml=None):
    xr = XMLreader(fpath,element_joins,element_aliases,contract_names,strip_prefix,warn,xml=xml)
    return xr.obj
#end def readxml
