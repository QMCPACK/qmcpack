##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  pseudopotential.py                                                #
#    Classes for reading pseudopotential data and converting         #
#    pseudopotentials between file formats.                          #
#                                                                    #
#  Content summary:                                                  #
#    Pseudopotentials                                                #
#      Class contains a list of all pseudopotentials available if    #
#      user provides pseudo_dir in settings.                         #
#                                                                    #
#    PseudoFile                                                      #
#      Represents a single pseudopotential file.  No functionality.  #
#                                                                    #
#    gamessPPFile                                                    #
#      Represents a pseudopotential file for GAMESS.                 #
#      Used during Nexus execution to separate basis and channel     #
#      information to fill in the input file.                        #
#                                                                    #
#    Classes below are user-facing, but are not used by Nexus itself.#
#                                                                    #
#    Pseudopotential                                                 #
#      Represents a generic pseudopotential.                         #
#                                                                    #
#    SemilocalPP                                                     #
#      Represents a semi-local pseudopotential.                      #
#      Contains data for each non-local channel.                     #
#      Supports plotting, gaussian fitting, and writing to QMCPACK's #
#        grid-based file format.                                     #
#      GaussianPP and QmcpackPP are derived classes to read/write    #
#        GAMESS/Gaussian and QMCPACK pseudopotential files.          #
#                                                                    #
#====================================================================#

import os
from subprocess import Popen
from execute import execute
import numpy as np
from numpy import linspace,array,zeros,append,mgrid,empty,exp,minimum,maximum,sqrt,arange

from fileio import TextFile
from xmlreader import readxml
from superstring import string2val,split_delims
from periodic_table import pt,is_element
from unit_converter import convert
from generic import obj
from developer import DevBase,unavailable,error
from basisset import process_gaussian_text,GaussianBasisSet
from physical_system import PhysicalSystem
from plotting import *
from debug import *
from testing import object_eq

try:
    import matplotlib.pyplot as plt
except:
    plt = unavailable('matplotlib','pyplot')
#end try


def pp_elem_label(filename,guard=False):
    el = ''
    for c in filename:
        if c=='.' or c=='_' or c=='-':
            break
        #end if
        el+=c
    #end for
    elem_label = el
    is_elem,symbol = is_element(el,symbol=True)
    if guard: 
        if not is_elem:
            error('cannot determine element for pseudopotential file: {0}\npseudopotential file names must be prefixed by an atomic symbol or label\n(e.g. Si, Si1, etc)'.format(filename))
        #end if
        return elem_label,symbol
    else:
        return elem_label,symbol,is_elem
    #end if
#end def pp_elem_label


# basic interface for nexus, only gamess really needs this for now
class PseudoFile(DevBase):
    def __init__(self,filepath=None):
        self.element       = None
        self.element_label = None
        self.filename      = None
        self.location      = None
        if filepath!=None:
            self.filename = os.path.basename(filepath)
            self.location = os.path.abspath(filepath)
            elem_label,symbol,is_elem = pp_elem_label(self.filename)
            if not is_elem:
                self.error('cannot determine element for pseudopotential file: {0}\npseudopotential file names must be prefixed by an atomic symbol or label\n(e.g. Si, Si1, etc)'.format(filepath))
            #end if
            self.element = symbol
            self.element_label = elem_label
            self.read(filepath)
        #end if
    #end def __init__

    def read(self,filepath):
        None
    #end def read
#end class PseudoFile



class gamessPPFile(PseudoFile):
    def __init__(self,filepath=None):
        self.pp_text    = None
        self.pp_name    = None
        self.basis_text = None
        PseudoFile.__init__(self,filepath)
    #end def __init__

    def read(self,filepath):
        lines = open(filepath,'r').read().splitlines()
        new_block  = True
        tokens     = []
        block      = ''
        nline = 0
        for line in lines:
            nline+=1
            ls = line.strip()
            if len(ls)>0 and ls[0]!='!' and ls[0]!='#':
                if new_block:
                    tokens = ls.split()
                    new_block = False
                    if len(tokens)!=5:
                        block+=line+'\n'
                    #end if
                else:
                    block+=line+'\n'
                #end if
            #end if
            if (len(ls)==0 or nline==len(lines)) and len(block)>0:
                block = block.rstrip()
                if len(tokens)==4:
                    self.pp_text = block
                    self.pp_name = tokens[0]
                elif len(tokens)==5:
                    self.basis_text = block
                else:
                    self.error('could not identify text block in {0} as pseudopotential or basis text\btext block:\n{1}'.format(self.filename,block))
                #end if
                new_block = True
                tokens    = []
                block     = ''
            #end if
        #end for
        if self.pp_text is None:
            self.error('could not find pseudopotential text in '+self.filename)
        #end if
        if self.basis_text is None:
            self.error('could not find basis text in '+self.filename)
        #end if
    #end def read
#end class gamessPPFile



class Pseudopotentials(DevBase):
    def __init__(self,*pseudopotentials):
        if len(pseudopotentials)==1 and isinstance(pseudopotentials[0],list):
            pseudopotentials = pseudopotentials[0]
        #end if
        ppfiles = []
        pps     = []
        errors = False
        for pp in pseudopotentials:
            if isinstance(pp,PseudoFile):
                pps.append(pp)
            elif isinstance(pp,str):
                ppfiles.append(pp)
            else:
                self.error('expected PseudoFile type or filepath, got '+str(type(pp)),exit=False)
                errors = True
            #end if
        #end for
        if errors:
            self.error('cannot create Pseudopotentials object')
        #end if
        if len(pps)>0:
            self.addpp(pps)
        #end if
        if len(ppfiles)>0:
            self.readpp(ppfiles)
        #end if
    #end def __init__


    def addpp(self,*pseudopotentials):
        if len(pseudopotentials)==1 and isinstance(pseudopotentials[0],list):
            pseudopotentials = pseudopotentials[0]
        #end if
        for pp in pseudopotentials:
            self[pp.filename] = pp
        #end for
    #end def addpp

        
    def readpp(self,*ppfiles):
        if len(ppfiles)==1 and isinstance(ppfiles[0],list):
            ppfiles = ppfiles[0]
        #end if
        pps = []
        self.log('\n  Pseudopotentials')
        for filepath in ppfiles:
            filename = os.path.basename(filepath)
            elem_label,symbol,is_elem = pp_elem_label(filename)
            is_file = os.path.isfile(filepath)
            if is_elem and is_file:
                self.log('    reading pp: ',filepath)
                ext = filepath.split('.')[-1].lower()
                if ext=='gms':
                    pp = gamessPPFile(filepath)
                else:
                    pp = PseudoFile(filepath)
                #end if
                pps.append(pp)
            elif not is_file:
                self.log('    ignoring directory: ',filepath)
            elif not is_elem:
                self.log('    ignoring file w/o atomic symbol: ',filepath)
            #end if
        #end for
        self.log(' ')
        self.addpp(pps)
    #end def readpp


    def pseudos_by_atom(self,*ppfiles):
        pps = obj()
        for ppfile in ppfiles:
            if ppfile in self:
                pp = self[ppfile]
                pps[pp.element_label] = pp
            else:
                self.error('pseudopotential file not found\nmissing file: {0}'.format(ppfile))
            #end if
        #end for
        return pps
    #end def pseudos_by_atom
#end class Pseudopotentials



# user interface to group sets of pseudopotentials together and refer to them by labels
#   labeling should eliminate the need to provide lists of pseudopotential files to each 
#   simulation object (e.g. via a generate_* call) separately
class PPset(DevBase):
    instance_counter = 0

    known_codes = set('pwscf gamess vasp qmcpack'.split())

    default_extensions = obj(
        pwscf   = ['ncpp','upf'],
        gamess  = ['gms'],
        vasp    = ['potcar'],
        qmcpack = ['xml'],
        )

    def __init__(self):
        if PPset.instance_counter!=0:
            self.error('cannot instantiate more than one PPset object\nintended use follows a singleton pattern')
        #end if
        PPset.instance_counter+=1
        self.pseudos = obj()
    #end def __init__

    def supports_code(self,code):
        return code in PPset.known_codes
    #end def supports_code

    def __call__(self,label,**code_pps):
        if not isinstance(label,str):
            self.error('incorrect use of ppset\nlabel provided must be a string\nreceived type instead: {0}\nwith value: {1}'.format(label.__class__.__name__,label))
        #end if
        if label in self.pseudos:
            self.error('incorrect use of ppset\npseudopotentials with label "{0}" have already been added to ppset'.format(label))
        #end if
        pseudos = obj()
        self.pseudos[label]=pseudos
        for code,pps in code_pps.items():
            clow = code.lower()
            if clow not in self.known_codes:
                self.error('incorrect use of ppset\ninvalid simulation code "{0}" provided with set labeled "{1}"\nknown simulation codes are: {2}'.format(code,label,sorted(self.known_codes)))
            #end if
            if not isinstance(pps,(list,tuple)):
                self.error('incorrect use of ppset\nmust provide a list of pseudopotentials for code "{0}" in set labeled "{1}"\ntype provided instead of list: {2}'.format(code,label,pps.__class__.__name__))
            #end if
            ppcoll = obj()
            for pp in pps:
                if not isinstance(pp,str):
                    self.error('incorrect use of ppset\nnon-filename provided with set labeled "{0}" for simulation code "{1}"\neach pseudopential file name must be a string\nreceived type: {2}\nwith value: {3}'.format(label,code,pp.__class__.__name__,pp))
                #end if
                elem_label,symbol,is_elem = pp_elem_label(pp)
                if not is_elem:
                    self.error('invalid filename provided to ppset\ncannot determine element for pseudopotential file: {0}\npseudopotential file names must be prefixed by an atomic symbol or label\n(e.g. Si, Si1, etc)'.format(pp))
                elif symbol in ppcoll:
                    self.error('incorrect use of ppset\nmore than one pseudopotential file provided for element "{0}" for code "{1}" in set labeled "{2}"\nfirst file: {3}\nsecond file: {4}'.format(symbol,code,label,ppcoll[symbol],pp))
                #end if
                ppcoll[symbol] = pp
            #end for
            pseudos[clow] = ppcoll
        #end for
    #end def __call__

    def has_set(self,label):
        return label in self.pseudos
    #end def has_set


    # test needed
    def get(self,label,code,system):
        if system is None or not system.pseudized:
            return []
        #end if
        if not isinstance(system,PhysicalSystem):
            self.error('system object must be of type PhysicalSystem')
        #end if
        species_labels,species = system.structure.species(symbol=True)
        if not isinstance(label,str):
            self.error('incorrect use of ppset\nlabel provided must be a string\nreceived type instead: {0}\nwith value: {1}'.format(label.__class__.__name__,label))
        #end if
        if not self.has_set(label):
            self.error('incorrect use of ppset\npseudopotential set labeled "{0}" is not present in ppset\nset labels present: {1}\nplease either provide pseudopotentials with label "{0}" or correct the provided label'.format(label,sorted(self.pseudos.keys())))
        #end if
        pseudos = self.pseudos[label]
        clow = code.lower()
        if clow not in self.known_codes:
            self.error('simulation code "{0}" is not known to ppset\nknown codes are: {1}'.format(code,sorted(self.known_codes)))
        elif clow not in pseudos:
            self.error('incorrect use of ppset\npseudopotentials were not provided for simulation code "{0}" in set labeled "{1}"\npseudopotentials are required for physical system with pseudo-elements: {2}\nplease add these pseudopotentials for code "{0}" in set "{1}"'.format(code,label,sorted(species)))
        #end if
        ppcoll = pseudos[clow]
        pps = []
        for symbol in species:
            if symbol not in ppcoll:
                self.error('incorrect use of ppset\npseudopotentials were not provided for element "{0}" code "{1}" in set labeled "{2}"\nphysical system encountered with pseudo-elements: {3}\nplease ensure that pseudopotentials are provided for these elements in set "{2}" for code "{1}"'.format(symbol,code,label,sorted(species)))
            #end if
            pps.append(ppcoll[symbol])
        #end for
        return pps
    #end def get
#end class PPset
ppset = PPset()











# real pseudopotentials
from plotting import *
show_plots = show
set_title  = title

class Pseudopotential(DevBase):

    requires_format = False
    formats = None

    def __init__(self,filepath=None,format=None):
        self.element = None
        self.core    = None
        self.Zval    = None
        self.Zcore   = None
        if filepath!=None:
            self.read(filepath,format)
        #end if
    #end def __init__

    
    def transfer_core_from(self,other):
        self.element = other.element
        self.core    = other.core   
        self.Zval    = other.Zval   
        self.Zcore   = other.Zcore  
    #end def transfer_core_from


    def read(self,filepath,format=None):
        if self.requires_format:
            if format is None:
                self.error('format keyword must be specified to read file {0}\nvalid options are: {1}'.format(filepath,self.formats))
            elif not format in self.formats:
                self.error('incorrect format requested: {0}\nvalid options are: {1}'.format(format,self.formats))
            #end if
        #end if
        if not os.path.exists(filepath):
            self.error('cannot read {0}, file does not exist'.format(filepath))
        #end if
        self.element = split_delims(os.path.split(filepath)[1])[0]
        text = open(filepath,'r').read()
        self.read_text(text,format,filepath=filepath)
    #end def read


    def write(self,filepath=None,format=None):
        if self.requires_format:
            if format is None:
                self.error('format keyword must be specified to write file {0}\nvalid options are: {1}'.format(filepath,self.formats))
            elif not format in self.formats:
                self.error('incorrect format requested: {0}\nvalid options are: {1}'.format(format,self.formats))
            #end if
        #end if
        text = self.write_text(format)
        if filepath!=None:
            open(filepath,'w').write(text)
        #end if
        return text
    #end def write


    def read_text(self,text,format=None,filepath=None):
        self.not_implemented()
    #end def read_text

    def write_text(self,format=None):
        self.not_implemented()
    #end def write_text

    def convert(self,format):
        self.not_implemented()
    #end def convert

    def plot(self,r=None,show=True):
        self.not_implemented()
    #end def plot
#end class Pseudopotential



class SemilocalPP(Pseudopotential):
    l_channels   = tuple('spdfghiklmnoqrtuvwxyz')
    channel_colors = obj(s='g',p='r',d='b',f='m',g='c',h='k',i='g',L2='k')

    numeric        = False
    interpolatable = True

    formats = ['qmcpack','casino']

    channel_indices = obj()
    for i,c in enumerate(l_channels):
        channel_indices[c] = i
    #end for
    channel_indices.L2 = -1

    def __init__(self,filepath=None,format=None,name=None,src=None):
        self.name = name
        self.rcut = None
        self.lmax = None
        self.local = None
        if self.numeric:
            self.r = None
        #end if
        self.rcut_L2 = None
        self.components = obj()
        Pseudopotential.__init__(self,filepath,format)
        if src!=None:
            self.transfer_core_from(src)
        #end if
    #end def __init__

    
    # test needed
    def transfer_core_from(self,other):
        self.name  = other.name
        self.rcut  = other.rcut
        self.lmax  = other.lmax
        self.local = other.local
        if self.numeric and other.numeric:
            self.r = other.r
        #end if
        Pseudopotential.transfer_core_from(self,other)
    #end def transfer_core_from


    def read(self,filepath,format=None):
        Pseudopotential.read(self,filepath,format)
        if self.rcut is None:
            self.update_rcut()
        #end if
    #end def read


    def has_component(self,l):
        return l in self.components
    #end def has_component


    # test needed
    def set_component(self,l,v,guard=False):
        if guard and l in self.components:
            self.error('cannot set requested component potential\nrequested potential is already present\nrequested potential: {0}'.format(l))
        #end if
        self.components[l] = v
    #end def set_component


    def get_component(self,l,guard=False):
        v = None
        if l in self.components:
            v = self.components[l]
        elif guard:
            self.error('cannot get requested component potential\nrequested potential is not present\nrequested potential: {0}\npotentials present: {1}'.format(l,list(self.components.keys())))
        #end if
        return v
    #end def get_component


    # test needed
    def remove_component(self,l,guard=False):
        if l in self.components:
            del self.components[l]
        elif guard:
            self.error('cannot remove requested component potential\nrequested potential is not present\nrequested potential: {0}\npotentials present: {1}'.format(l,list(self.components.keys())))
        #end if
    #end def remove_component


    def has_local(self):
        return self.local in self.components
    #end def has_local


    def has_nonlocal(self,l=None):
        vnl = self.get_nonlocal()
        if l is None:
            return len(vnl)>0
        else:
            return l in vnl
        #end if
    #end def has_nonlocal


    def has_L2(self):
        return 'L2' in self.components
    #end def has_L2


    def get_local(self):
        return self.get_component(self.local,guard=True)
    #end def get_local


    def get_nonlocal(self,l=None):
        vnl = obj()
        for lc,vc in self.components.items():
            if lc!=self.local and lc!='L2':
                vnl[lc] = vc
            #end if
        #end for
        if l is None:
            return vnl
        elif l in vnl:
            return vnl[l]
        else:
            self.error('cannot get nonlocal potential\nrequested potential is not nonlocal\nrequested potential: {0}\nnonlocal potentials present: {1}'.format(l,list(vnl.keys())))
        #end if
    #end def get_nonlocal


    # test needed
    def get_L2(self):
        return self.get_component('L2',guard=True)
    #end def get_L2


    # test needed
    def add_local(self,l,v):
        self.set_component(l,v,guard=True)
        self.local = l
    #end def add_local


    # test needed
    def add_nonlocal(self,l,v):
        if l==self.local:
            self.promote_local()
        #end if
        self.set_component(l,v,guard=True)
    #end def add_nonlocal


    # test needed
    def add_L2(self,v):
        self.set_component('L2',v,guard=True)
    #end def add_L2


    # test needed
    def remove_local(self):
        self.remove_component(self.local,guard=True)
        self.local = None
    #end def remove_local


    # test needed
    def remove_nonlocal(self,l=None):
        vnl = self.get_nonlocal()
        if l is None:
            for l in vnl.keys():
                self.remove_component(l,guard=True)
            #end for
        elif l in vnl:
            self.remove_component(l,guard=True)
        else:
            self.error('cannot remove nonlocal potential\nrequested potential is not present\nrequested potential: {0}\nnonlocal potentials present: {1}'.format(l,list(vnl.keys())))
        #end if
    #end def remove_nonlocal


    # test needed
    def remove_L2(self):
        self.remove_component('L2',guard=True)
    #end def remove_L2


    def assert_numeric(self,loc):
        if not self.numeric:
            self.error('failing at {0}\n{0} is only supported for numerical pseudopotential formats'.format(loc))
        #end if
    #end def assert_numeric


    # test needed
    def change_local(self,local):
        self.assert_numeric('change_local')
        if local==self.local:
            return
        #end if
        if self.has_component('L2'):
            self.error('cannot change local channel\nL2 potential is present')
        elif not self.has_component(self.local):
            self.error('cannot change local potential\ncurrent local potential is not present\ncurrent local potential: {0}\npotentials present: {1}'.format(self.local,list(self.components.keys())))
        elif not self.has_component(local):
            self.error('cannot change local potential\nrequested local potential is not present\nrequested local potential: {0}\npotentials present: {1}'.format(local,list(self.components.keys())))
        #end if
        vcs = self.components
        vloc = vcs[self.local]
        # get the l channels
        vls = obj()
        for l,vc in self.components.items():
            if l==self.local:
                vls[l] = vloc
            else:
                vls[l] = vc+vloc
            #end if
        #end for
        # from l channels, reconstruct nonlocal components
        vcs.clear()
        vloc = vls[local]
        for l,vl in vls.items():
            if l==local:
                vcs[l] = vloc
            else:
                vcs[l] = vls[l]-vloc
            #end if
        #end for
        self.local = local
    #end def change_local


    # test needed
    def promote_local(self):
        found = False
        for l in self.l_channels:
            if l not in self.components:
                found = True
                break
            #end if
        #end for
        if not found:
            self.error('could not promote local channel')
        #end if
        vloc = self.components[self.local]
        del self.components[self.local]
        self.local = l
        self.components[self.local] = vloc
    #end def promote_local


    # test needed
    # ensure that v_l=<l|vpp|l>==v, while v_l' remains unchanged
    def set_channel(self,l,v):
        self.assert_numeric('set_channel')
        if not self.has_local():
            self.error('cannot enforce channel matching via set_channel\nthe local potential is missing and must be present\nrequested channel: {0}'.format(l))
        #end if
        if l==self.local:
            self.promote_local()
        #end if
        vloc = self.get_local()
        vnl  = self.get_nonlocal()
        if self.has_L2():
            vL2 = self.get_L2()
        else:
            vL2 = 0*vloc
        #end if
        li = self.channel_indices[l]
        vl = vloc + li*(li+1)*vL2
        self.components[l] = v-vl
    #end def set_channel


    # test needed
    def expand_L2(self,lmax):
        self.assert_numeric('expand_L2')
        if lmax not in self.channel_indices:
            self.error('cannot expand L2 up to angular momentum "{0}"\nvalid options for lmax are: {1}'.format(lmax,self.l_channels))
        #end if
        limax = self.channel_indices[lmax]
        vps = obj()
        for li in range(limax+1):
            l = self.l_channels[li]
            vps[l] = self.evaluate_channel(l=l,rpow=1)
        #end for
        if self.has_L2():
            self.remove_L2()
        #end if
        self.components[self.local] = vps[self.local]
        for l,v in vps.items():
            if l!=self.local:
                self.set_channel(l,v)
            #end if
        #end for
    #end def expand_L2


    # test needed
    def angular_channels(self):
        channels = []
        for l in self.l_channels:
            if l in self.components:
                channels.append(l)
            #end if
        #end for
        return channels
    #end def angular_channels
        

    # evaluate r*potential based on a potential component object
    #  component representation is specific to each derived class
    def evaluate_comp_rV(self,r,l,vcomp):
        self.not_implemented()
    #end def evaluate_comp_rV


    def find_r_rng(self,r,rmin):
        if r is None and self.numeric and not self.interpolatable:
            r = self.r
        #end if
        rng = None
        if rmin>1e-12:
            rng = r>rmin
            r = r[rng]
        #end if
        return r,rng
    #end def find_r_rng


    # evaluate potential based on a potential component object
    #  local, nonlocal, and L2 are all represented by separate component objects
    def evaluate_comp(self,r,l,vcomp,rpow=0,rmin=0,rret=True):
        v = self.evaluate_comp_rV(r,l,vcomp)
        r,rng = self.find_r_rng(r,rmin)
        if rng is not None:
            v = v[rng]
        #end if
        if rpow!=1:
            v = r**(rpow-1)*v
        #end if
        if rret:
            return r,v
        else:
            return v
        #end if
    #end def evaluate_comp


    # evaluate the local component potential only
    def evaluate_local(self,r=None,rpow=0,rmin=0,rret=False):
        l = self.local
        if not self.has_component(l):
            self.error('cannot evaluate local potential\nlocal potential is not present')
        #end if
        vcomp = self.get_component(l)
        ret = self.evaluate_comp(r,l,vcomp,rpow,rmin,rret)
        return ret
    #end def evaluate_local


    # evaluate a nonlocal component potential
    def evaluate_nonlocal(self,r=None,l=None,rpow=0,rmin=0,rret=False):
        if l==self.local:
            self.error('called evaluate_nonlocal requesting local potential\nthe l index of the local potential is: {0}'.format(self.local))
        elif l=='L2':
            self.error('called evaluate_nonlocal requesting L2 potential')
        elif l is None:
            self.error('called evaluate_nonlocal without specifying the angular channel')
        elif not self.has_component(l):
            self.error('cannot evaluate non-local potential\nlocal potential is not present\nrequested potential: {0}'.format(l))
        #end if
        vcomp = self.get_component(l)
        ret = self.evaluate_comp(r,l,vcomp,rpow,rmin,rret)
        return ret
    #end def evaluate_nonlocal


    # test needed
    # evaluate the L2 component potential
    def evaluate_L2(self,r=None,rpow=0,rmin=0,rret=False):
        l = 'L2'
        if not self.has_component(l):
            self.error('cannot evaluate L2 potential\nL2 potential is not present')
        #end if
        vcomp = self.get_component(l)
        ret = self.evaluate_comp(r,l,vcomp,rpow,rmin,rret)
        return ret
    #end def evaluate_L2


    # evaluate semilocal potential components in isolation
    def evaluate_component(self,r=None,l=None,rpow=0,rmin=0,rret=False,optional=False):
        vcomp = self.get_component(l)
        if vcomp is not None:
            return self.evaluate_comp(r,l,vcomp,rpow,rmin,rret)
        elif optional:
            z = 0*self.evaluate_local(r,rpow,rmin,rret=False)
            if rret:
                return z,z
            else:
                return z
            #end if
        else:
            self.error('requested evaluation of non-existent component\ncomponent requested: {0}'.format(l))
        #end if
    #end def evaluate_component


    # evaluate angular momentum channel of full potential
    def evaluate_channel(self,r=None,l=None,rpow=0,rmin=0,rret=False,with_local=True,with_L2=True):
        if l not in self.l_channels:
            self.error('evaluate_channel must be called with a valid angular momentum label\nvalid options are l=s,p,d,f,...\nyou provided: l={0}'.format(l))
        #end if
        eval_any = False
        loc_present = self.has_component(self.local)
        l_present   = self.has_component(l)
        L2_present  = self.has_component('L2')
        if (l==self.local or with_local) and loc_present:
            vloc = self.evaluate_local(r,rpow,rmin)
            eval_any = True
        else:
            vloc = 0
        #end if
        if with_L2 and L2_present:
            l_int = self.channel_indices[l]
            vL2 = self.evaluate_L2(r,rpow,rmin)*l_int*(l_int+1)
            eval_any = True
        else:
            vL2 = 0
        #end if
        if l!=self.local and l_present:
            vnonloc = self.evaluate_nonlocal(r,l,rpow,rmin)
            eval_any = True
        else:
            vnonloc = 0
        #end if
        if rret or not eval_any:
            r,rng = self.find_r_rng(r,rmin)
        #end if
        if eval_any:
            v = vloc+vL2+vnonloc
        else:
            v = 0*r
        #end if
        if rret:
            return r,v
        else:
            return v
        #end if
    #end def evaluate_channel


    # similar to evaluate_channel but with defaults appropriate for QMCPACK
    def numeric_channel(self,l=None,rmin=0.,rmax=10.,npts=10001,rpow=0,with_local=True,with_L2=True):
        if self.numeric and not self.interpolatable:
            r = None
        else:
            r = linspace(rmin,rmax,npts)
        #end if
        if l=='L2':
            r,v = self.evaluate_L2(r,rpow,rmin,rret=True)
        else:
            r,v = self.evaluate_channel(r,l,rpow,rmin,rret=True,with_local=with_local,with_L2=with_L2)
        #end if
        return r,v
    #end def numeric_channel


    def update_rcut(self,tol=1e-5,optional=False):
        if not optional or self.rcut is None:
            self.rcut = self.find_rcut(tol=tol,with_L2=False)
        #end if
        has_vL2 = self.has_component('L2')
        if has_vL2 and (not optional or self.rcut_L2 is None):
            self.rcut_L2 = self.find_rcut(tol=tol,with_L2=True)
        #end if
        return self.rcut
    #end def update_rcut


    def find_rcut(self,tol=1e-5,with_L2=False):
        vnl = self.get_nonlocal()
        if with_L2 and self.has_L2():
            vnl.L2 = self.get_L2()
        #end if
        if len(vnl)==0:
            rmax = 10.
            return rmax
        #end if
        channels = list(vnl.keys())
        rv = obj()
        # add a zero potential
        l = channels[0]
        rvl = self.numeric_channel(l,rmin=0.01,with_local=False,with_L2=l=='L2')
        rv[l] = rvl
        rv['0'] = rvl[0],0*rvl[1]
        for l in channels[1:]:
            rv[l] = self.numeric_channel(l,rmin=0.01,with_local=False,with_L2=l=='L2')
        #end for
        r    = None
        vmin = None
        vmax = None
        for l,(rc,vc) in rv.items():
            if r is None:
                r = rc
                vmin = array(vc)
                vmax = array(vc)
            elif len(rc)!=len(r):
                self.error('numeric representation of channels do not match in length')
            else:
                vmin = minimum(vmin,vc)
                vmax = maximum(vmax,vc)
            #end if
        #end for
        vspread = vmax-vmin
        rcut = r[-1]
        nr = len(r)
        for i in range(nr):
            n = nr-1-i
            if vspread[n]>tol:
                rcut = r[n]
                break
            #end if
        #end for

        #figure()
        #plot(r,vspread,'k-')
        #plot([rcut,rcut],[0,vspread[1:].max()],'k--')
        #title('rcut = {0}'.format(rcut))
        #show()

        return rcut
    #end def find_rcut


    def plot(self,r=None,show=True,fig=True,linestyle='-',channels=None,with_local=False,rmin=0.01,rmax=5.0,title=None,metric=None,color=None):
        if channels is None:
            channels = self.l_channels
        #end if
        if fig:
            figure()
        #end if
        if r is None and self.numeric:
            r = self.r
        elif r is None:
            r = linspace(rmin,rmax,1000)
        #end if
        rin = r
        color_in = color
        for c in channels:
            r = rin
            color = color_in
            if c in self.components:
                if color is None:
                    color = self.channel_colors[c]
                #end if
                lab = c
                if c==self.local:
                    lab = self.local+' loc'
                #end if
                if self.name!=None:
                    lab = self.name+' '+lab
                #end if
                v = self.evaluate_channel(r,c,with_local=with_local,rmin=rmin-1e-12)
                rng = r>rmin-1e-12
                r = r[rng]
                if metric=='r2':
                    v = r**2*v
                    if c==self.local:
                        v += self.Zval*r
                    #end if
                elif metric!=None:
                    self.error('invalid metric for plotting: {0}\nvalid options are: r2'.format(metric))
                #end if
                plot(r,v,color+linestyle,label=lab)
            #end for
        #end for
        if fig:
            if title is None:
                title = 'Semilocal {0} PP ({1} core)'.format(self.element,self.core)
            #end if
            set_title(title)
            ylabel('channel potentials (Ha)')
            xlabel('r (Bohr)')
            legend()
        #end if
        if show:
            show_plots()
        #end if
    #end def plot


    def plot_components(self,r=None,show=True,fig=True,linestyle='-',rmin=0.01,rmax=5.0,title=None,metric=None,color=None,rpow=0):
        channels = list(self.l_channels)+['L2']
        if fig:
            plt.figure(tight_layout=True)
        #end if
        if r is None and self.numeric:
            r = self.r
            #r = None
        elif r is None:
            r = linspace(rmin,rmax,1000)
        #end if
        rin = r
        color_in = color
        for c in channels:
            r = rin
            color = color_in
            channel = self.get_component(c)
            if channel is not None:
                if color is None:
                    color = self.channel_colors[c]
                #end if
                lab = c
                if c==self.local:
                    lab += ' loc'
                elif c!='L2':
                    lab += '-'+self.local
                #end if
                if self.name!=None:
                    lab = self.name+' '+lab
                #end if
                v = self.evaluate_component(r,c,rpow,rmin-1e-12)
                rng = r>rmin-1e-12
                r = r[rng]
                if metric=='r2':
                    v = r**2*v
                    if c==self.local:
                        v += self.Zval*r
                    #end if
                elif metric!=None:
                    self.error('invalid metric for plotting: {0}\nvalid options are: r2'.format(metric))
                #end if
                plt.plot(r,v,color+linestyle,label=lab)
            #end for
        #end for
        if fig:
            if title is None:
                title = 'Semilocal {0} PP ({1} core)'.format(self.element,self.core)
            #end if
            plt.title(title)
            plt.ylabel('component potentials (Ha)')
            plt.xlabel('r (Bohr)')
            plt.legend()
        #end if
        if show:
            plt.show()
        #end if
    #end def plot_components


    def plot_channels(self,r=None,channels=None,show=True,fig=True,linestyle='-',rmin=0.01,rmax=5.0,title=None,metric=None,color=None,rpow=0,with_local=True,with_L2=True):
        if channels is None:
            channels = list(self.l_channels)
        #end if
        if fig:
            figure(tight_layout=True)
        #end if
        if r is None and self.numeric:
            r = self.r
            #r = None
        elif r is None:
            r = linspace(rmin,rmax,1000)
        #end if
        rin = r
        color_in = color
        if not self.has_component('L2'):
            loc_label = self.local
            for c in channels:
                if c not in self.components:
                    loc_label+=c
                #end if
            #end for
            for c in channels:
                r = rin
                color = color_in
                if self.has_component(c):
                    if color is None:
                        color = self.channel_colors[c]
                    #end if
                    lab = c
                    if c==self.local:
                        lab = loc_label
                    #end if
                    if self.name!=None:
                        lab = self.name+' '+lab
                    #end if
                    v = self.evaluate_channel(r,c,rpow,rmin-1e-12,False,with_local,with_L2)
                    rng = r>rmin-1e-12
                    r = r[rng]
                    if metric=='r2':
                        v = r**2*v
                    elif metric!=None:
                        self.error('invalid metric for plotting: {0}\nvalid options are: r2'.format(metric))
                    #end if
                    plot(r,v,color+linestyle,label=lab)
                #end for
            #end for
        else:
            for c in channels:
                r = rin
                color = color_in
                if color is None:
                    color = self.channel_colors[c]
                #end if
                lab = c
                if self.name!=None:
                    lab = self.name+' '+lab
                #end if
                v = self.evaluate_channel(r,c,rpow,rmin-1e-12,False,with_local,with_L2)
                rng = r>rmin-1e-12
                r = r[rng]
                if metric=='r2':
                    v = r**2*v
                elif metric!=None:
                    self.error('invalid metric for plotting: {0}\nvalid options are: r2'.format(metric))
                #end if
                plot(r,v,color+linestyle,label=lab)
            #end for
        #end if
        if fig:
            if title is None:
                title = 'Semilocal {0} PP angular channels ({1} core)'.format(self.element,self.core)
            #end if
            set_title(title)
            ylabel('channels')
            xlabel('r')
            legend()
        #end if
        if show:
            show_plots()
        #end if
    #end def plot_channels


    def plot_positive_definite(self,r=None,show=True,fig=True,linestyle='-',rmin=0.01,rmax=5.0,title=None,color='k'):
        if not self.has_L2():
            self.error('positive definite condition only applies to L2 potentials')
        #end if
        if fig:
            figure(tight_layout=True)
        #end if
        if r is None and self.numeric:
            r = self.r
        elif r is None:
            r = linspace(rmin,rmax,1000)
        #end if
        vL2 = self.evaluate_L2(r,0,rmin-1e-12)
        rng = r>rmin-1e-12
        r = r[rng] 
        b = vL2*(2*r**2)
        plot(r,1+b,color+linestyle,label='1+b')
        plot(r,0*r,'r-')
        if fig:
            if title is None:
                title = 'L2 positive definite condition {0} PP ({1} core)'.format(self.element,self.core)
            #end if
            set_title(title)
            ylabel('1+b > 0')
            xlabel('r (Bohr)')
            legend()
        #end if
        if show:
            show_plots()
        #end if
    #end def plot_positive_definite

                
    def plot_L2(self,show=True,fig=True,r=None,rmin=0.01,rmax=5.0,linestyle='-',title=None,color=None):
        if fig:
            plt.figure(tight_layout=True)
        #end if
        if r is None and self.numeric:
            r = self.r
        elif r is None:
            r = linspace(rmin,rmax,1000)
        #end if
        vs = self.evaluate_channel(r,'s',with_local=True,rmin=rmin-1e-12)
        for c in self.l_channels[1:]:
            if c in self.components:
                if color is None:
                    color = self.channel_colors[c]
                #end if
                v = self.evaluate_channel(r,c,with_L2=False,rmin=rmin-1e-12)
                rng = r>rmin-1e-12
                r = r[rng]
                l = self.channel_indices[c]
                vL2 = (v-vs)/(l*(l+1))
                plt.plot(r,vL2,color+linestyle,label='(v{0}-vs)/(l(l+1))'.format(c))
            #end if
        #end for
        if fig:
            plt.xlim([0,rmax])
            if title is None:
                title = 'Semilocal {0} PP ({1} core)'.format(self.element,self.core)
            #end if
            plt.title(title)
            plt.ylabel('vL2 for channels above s')
            plt.xlabel('r')
            plt.legend()
            if show:
                plt.show()
            #end if
        #end if 
    #end def plot_L2


    def plot_nonlocal_polar(self,show=True,lmax=10,rmin=0.01,rmax=2.0,nr=100,nt=100,levels=100,label=''):
        from scipy.special import eval_legendre as legendre

        tlabel = label
        rpow = 1

        # update rcut, if needed
        if self.rcut is None:
            self.update_rcut()
        #end if
        rc = self.rcut

        # get legendre polynomials
        th = np.linspace(0.0,2*np.pi,nt)
        cos = np.cos(th)
        sin = np.sin(th)
        Pl  = obj()
        Pli = obj()
        for li in range(lmax):
            pl_i =  (2*li+1)/(4*np.pi)*legendre(li,cos)
            Pli[li] = pl_i
            if li<len(self.l_channels):
                l = self.l_channels[li]
                Pl[l] = pl_i
            #end if
        #end for

        # set the colormap and centre the colorbar
        import matplotlib.colors as colors
        class MidNorm(colors.Normalize):
            def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
                self.midpoint = midpoint
                colors.Normalize.__init__(self, vmin, vmax, clip)
            #end def __init__
            def __call__(self, value, clip=None):
                x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
                return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
            #end def __call__
        #end class MidNorm

        # select a subset of the radial points
        r = None
        vl = obj()
        vnl = self.get_nonlocal()
        if self.numeric and not self.interpolatable:
            for l in vnl.keys():
                rf,vf = self.evaluate_nonlocal(l=l,rpow=rpow,rret=True)
                if r is None:
                    drf = rf[1]-rf[0]
                    dr  = rmax/nr
                    ndrf = int(np.round(np.ceil(dr/drf)))
                    dr  = drf*ndrf
                    rmax = nr*dr + 0.1*drf
                    rng = rf<rmax
                    r = rf[rng][::ndrf]
                #end if
                vl[l] = vf[rng][::ndrf]
            #end for
        else:
            self.error('plot_polar does not yet support non-numeric potentials')
        #end if

        # plot the radial potentials
        figure()
        vmin = 1e99
        vmax = -1e99
        for l in self.l_channels:
            if l in vl:
                color = self.channel_colors[l]
                v = vl[l]
                plot(r,v,color+'-',label='v'+l)
                vmin = min(v.min(),vmin)
                vmax = max(v.max(),vmax)
            #end if
        #end for
        plot([rc,rc],[vmin,vmax],'k--',lw=2)
        xlim([-0.1*rc,1.1*rc])
        xlabel('r (Bohr)')
        ylabel('V NL (Ha)')
        title((tlabel+'  NL channels').strip())
        legend()

        # function for a single polar plot
        def plot_V(V,label):
            vmin = V.min()
            vmax = V.max()
            vm = max(np.abs(vmin),np.abs(vmax))
            vmin = -vm
            vmax = vm
            lev = linspace(vmin,vmax,levels)

            fig = figure(tight_layout=True)
            ax = fig.add_subplot(111)
            ax.set_xlabel('x')
            ax.set_ylabel('z')
            ax.set_aspect('equal','box')
            fig.tight_layout()

            xlim(lim)
            ylim(lim)

            cmap = plt.cm.get_cmap('seismic')

            mid_norm = MidNorm(vmin,vmax,0.0)

            cs = ax.contourf(X,Z,V,levels=lev,cmap=cmap,clim=(vmin,vmax),norm=mid_norm)
            plot(rc*cos,rc*sin,'k--',lw=2)

            fig.colorbar(cs, ax=ax, shrink=0.9)

            title((tlabel+'  V {}'.format(label)).strip())
        #end def plot_V

        # make a polar plot of each non-local component
        R,COS = np.array(np.meshgrid(r,cos,indexing='ij'))
        Z = R*COS
        R,SIN = np.array(np.meshgrid(r,sin,indexing='ij'))
        X = R*SIN
        lim = [-1.1*rc,1.1*rc]
        VSUM = None
        for l in self.l_channels:
            if l in vl:
                VL,PL = np.array(np.meshgrid(vl[l],Pl[l],indexing='ij'))
                V = VL*PL

                if VSUM is None:
                    VSUM = V.copy()
                else:
                    VSUM += V
                #end if

                plot_V(V,'NL '+l)
            #end if
        #end for

        # plot the total non-local operator
        plot_V(VSUM,'NL sum')

        # plot approximation to L2 operator, if present
        if self.has_L2():
            vL2 = self.evaluate_L2(rpow=rpow)
            vL2 = vL2[rng][::ndrf]
            VL2SUM = None
            for li in range(1,lmax):
                #print('V L2 {} of {}'.format(li,lmax))
                v = li*(li+1)*vL2
                VL,PL = np.array(np.meshgrid(v,Pli[li],indexing='ij'))
                V = VL*PL
                if VL2SUM is None:
                    VL2SUM = V.copy()
                else:
                    VL2SUM += V
                #end if
                #plot_V(V,'L2 '+str(li))
            #end for
            plot_V(VL2SUM,'L2 sum (Lmax={})'.format(lmax))
        #end if

        if show:
            show_plots()
        #end if

    #end def plot_nonlocal_polar


    def write_qmcpack(self,filepath=None):
        self.update_rcut(tol=1e-5,optional=True)

        channels = self.angular_channels()

        symbol        = self.element
        atomic_number = self.Zcore+self.Zval
        zval          = self.Zval
        creator       = 'Nexus'
        npots_down    = len(channels)
        l_local       = self.channel_indices[self.local]
        if l_local == -1:
            self.error('Local channel, {}, not coded.'.format(self.local))
        #end if


        rmin = 1e99
        rmax = -1e99
        npts = 0
        vps = obj()
        for l in channels:
            r,v = self.numeric_channel(l,rpow=1,with_local=True,with_L2=False)
            rmin   = min(rmin,r.min())
            rmax   = max(rmax,r.max())
            npts   = len(r)
            vps[l] = v
        #end for

        header = '''<?xml version="1.0" encoding="UTF-8"?>
<pseudo version="0.5">
  <header symbol="{0}" atomic-number="{1}" zval="{2}" relativistic="unknown" 
   polarized="unknown" creator="{3}" flavor="unknown" 
   core-corrections="unknown" xc-functional-type="unknown" 
   xc-functional-parametrization="unknown"/>
'''.format(symbol,atomic_number,zval,creator)

        grid = '  <grid type="linear" units="bohr" ri="{0}" rf="{1}" npts="{2}"/>\n'.format(rmin,rmax,npts)
        L2 = ''
        if self.has_component('L2'):
            dpad = '\n      '
            L2 += '  <L2 units="hartree" format="r*V" cutoff="{0}">\n'.format(self.rcut_L2)
            L2 += '    <radfunc>\n'
            L2 += '    '+grid
            L2 += '      <data>'
            r,v = self.numeric_channel('L2',rpow=1)
            n=0
            for d in v:
                if n%3==0:
                    L2 += dpad
                #end if
                L2 += '{0:22.14e}'.format(d)
                n+=1
            #end for
            L2 = L2.rstrip()+'\n'
            L2 += '      </data>\n'
            L2 += '    </radfunc>\n'
            L2 += '  </L2>\n'
        #end if
        semilocal =   '  <semilocal units="hartree" format="r*V" npots-down="{0}" npots-up="0" l-local="{1}">\n'.format(npots_down,l_local)
        dpad = '\n        '
        for l in self.l_channels:
            if l in vps:
                semilocal+='    <vps principal-n="0" l="{0}" spin="-1" cutoff="{1}" occupation="unknown">\n'.format(l,self.rcut)
                semilocal+='      <radfunc>\n'
                semilocal+='      '+grid
                semilocal+='        <data>'
                v = vps[l]
                n=0
                for d in v:
                    if n%3==0:
                        semilocal+=dpad
                    #end if
                    semilocal+='{0:22.14e}'.format(d)
                    n+=1
                #end for
                semilocal = semilocal.rstrip()+'\n'
                semilocal+='        </data>\n'
                semilocal+='      </radfunc>\n'
                semilocal+='    </vps>\n'
            #end if
        #end for
        semilocal+='  </semilocal>\n'
        footer = '</pseudo>\n'

        text = header+grid+L2+semilocal+footer

        if filepath!=None:
            open(filepath,'w').write(text)
        #end if
        return text
    #end def write_qmcpack


    def write_casino(self,filepath=None):
        if self.has_component('L2'):
            self.error('cannot write potential in CASINO format\nan L2 term is present, but this is not supported by CASINO')
        #end if

        channels = self.angular_channels()

        name          = self.name
        symbol        = self.element
        atomic_number = self.Zcore+self.Zval
        zval          = float(self.Zval)
        l_local       = 'spdfgi'.find(self.local)

        if name is None:
            name = '{0} pseudopotential converted by Nexus'.format(symbol)
        #end if

        rmin = 1e99
        rmax = -1e99
        npts = 0
        vps = obj()
        for l in channels:
            r,v = self.numeric_channel(l,rpow=1,with_local=True,with_L2=False)
            rmin   = min(rmin,r.min())
            rmax   = max(rmax,r.max())
            npts   = len(r)
            vps[l] = v
        #end for

        header = '''{0}
Atomic number and pseudo-charge
  {1} {2}
Energy units (rydberg/hartree/ev):
  hartree
Angular momentum of local component (0=s,1=p,2=d..)
  {3}
NLRULE override (1) VMC/DMC (2) config gen (0 ==> input/default value)
  0 0
Number of grid points
  {4}
'''.format(name,atomic_number,zval,l_local,npts)

        grid = 'R(i) in atomic units\n'
        for d in r:
            grid += '  {0:20.14e}\n'.format(d)
        #end for

        channels = ''
        for l in self.l_channels:
            if l in vps:
                channels += 'r*potential (L={0}) in Ha\n'.format('spdfgi'.find(l))
                v = vps[l]
                for d in v:
                    channels += '  {0:20.14e}\n'.format(d)
                #end for
            #end if
        #end for
        text = header+grid+channels

        if filepath!=None:
            open(filepath,'w').write(text)
        #end if
        return text
    #end def write_casino

#end class SemilocalPP





class GaussianPP(SemilocalPP):
    requires_format = True
    formats = SemilocalPP.formats + 'gaussian gamess crystal numhf'.split()

    @staticmethod
    def process_float(s):
        return float(s.replace('D','e').replace('d','e'))
    #end def process_float

    def __init__(self,filepath=None,format=None,name=None,src=None):
        self.basis = None
        SemilocalPP.__init__(self,filepath,format,name,src)
    #end def __init__


    # test needed for gaussian and crystal
    def read_text(self,text,format=None,filepath=None):
        lines,basis_lines = process_gaussian_text(text,format)

        format=format.lower()
        channels = []
        basis    = None
        # need lmax, element, and Zcore
        if format=='gamess':
            i=0
            name,type,Zcore,lmax = lines[i].split(); i+=1
            Zcore = int(Zcore)
            lmax  = int(lmax)
            element = split_delims(name)[0]
            while i<len(lines):
                n = int(lines[i]); i+=1
                terms = []
                for j in range(n):
                    coeff,rpow,expon = lines[i].split(); i+=1
                    terms.append((float(coeff),int(rpow),float(expon)))
                #end for
                channels.append(terms)
            #end while
        elif format=='gaussian':
            i=0
            element,token = lines[i].split(); i+=1
            label,lmax,Zcore = lines[i].split(); i+=1
            lmax = int(lmax)
            Zcore = int(Zcore)
            while i<len(lines):
                i+=1 # skip comment line
                n = int(lines[i]); i+=1
                terms = []
                for j in range(n):
                    rpow,expon,coeff = lines[i].split(); i+=1
                    terms.append((float(coeff),int(rpow),float(expon)))
                #end for
                channels.append(terms)
            #end while
        elif format=='crystal':
            i = 0
            conv_atomic_number,nshells = lines[i].split(); i+=1
            if len(conv_atomic_number)==1:
                atomic_number = int(conv_atomic_number)
            else:
                atomic_number = int(conv_atomic_number[-2:])
            #end if
            element = pt.simple_elements[atomic_number].symbol
            if 'input' not in lines[i].lower():
                self.error('INPUT must be present for crystal pseudpotential read')
            #end if
            i+=1
            tokens = lines[i].split()
            Zval = int(float(tokens[0])); i+=1
            Zcore = atomic_number-Zval
            nterms = array(tokens[1:],dtype=int)
            lmax = 0
            if nterms[0]==0:
                lmax+=1
                channels.append([(0.0,2,1.0)])
                nterms = nterms[1:]
            #end if
            for nt in nterms:
                lmax += 1
                terms = []
                for n in range(nt):
                    expon,coeff,rpow = lines[i].split(); i+=1
                    terms.append((float(coeff),int(rpow)+2,float(expon)))
                #end for
                channels.append(terms)
            #end for
            lmax-=1
        elif format=='atomscf':
            #i=0
            #self.name = lines[i].strip(); i+=1
            i=1 # skip title line
            element = 'Rn' # text does not contain element (must be corrected downstream)
            lmax    = -1   
            Zcore = int(lines[i].strip()); i+=1
            while i<len(lines):
                n = int(lines[i]); i+=1
                terms = []
                for j in range(n):
                    rpow,expon,coeff = lines[i].split(); i+=1
                    terms.append((float(coeff),int(rpow),float(expon)))
                #end for
                channels.append(terms)
                lmax+=1
            #end while
            # PPs in atomscf input are s,p,d,f, etc
            #  rather than, e.g. f,s-f,p-f,d-f
            #  so rearrange into form similar to f,s-f,p-f,d-f
            loc = channels.pop()
            for l in range(lmax):
                c = channels[l]
                channels[l] = c[0:len(c)-len(loc)]
            #end for
            channels = [loc] + channels
        elif format=='numhf':
            name = None
            i=0
            Zval,lmax = lines[i].split(); i+=1;
            Zval = int(Zval)
            lmax = int(lmax)-1
            element = self.element
            Zcore = int(pt[element].atomic_number)-Zval
            ns =  [int(n) for n in lines[i].split()]; i+= 1;
            while i<len(lines):
                for n in ns:
                    terms = []
                    for j in range(n):
                        rpow,expon,coeff = lines[i].split(); i+=1
                        terms.append((float(coeff),int(rpow),float(expon)))
                    #end for
                    channels.append(terms)
                #end for
            #end while
            # Bring local channel to front
            channels.insert(0,channels.pop())
        else:
            self.error('ability to read file format {0} has not been implemented'.format(format))
        #end if

        if basis_lines!=None:
            bs = GaussianBasisSet()
            bs.read_lines(basis_lines,format)
            basis = bs.basis
        #end if

        if not element in pt:
            if not self.element in pt:
                self.error('cannot identify element for pseudopotential file '+filepath)
            #end if
        else:
            self.element = element
        #end if
        Zatom = pt[element].atomic_number
        Zval = Zatom-Zcore
        if Zcore==0:
            core = None
        else:
            core = pt.simple_elements[Zcore].symbol
        #end if
        self.set(
            core    = core,
            Zval    = Zval,
            Zcore   = Zcore,
            lmax    = lmax
            )
        for c in range(len(channels)):
            if c==0:
                #cname = 'loc'
                cname = self.l_channels[lmax]
                self.local = cname
            else:
                cname = self.l_channels[c-1]
            #end if
            channel = obj()
            terms = channels[c]
            for t in range(len(terms)):
                coeff,rpow,expon = terms[t]
                channel[t] = obj(coeff=coeff,rpow=rpow,expon=expon)
            #end for
            self.components[cname] = channel
        #end for
        self.basis = basis
        if len(self.components)!=self.lmax+1:
            self.error('number of channels is not lmax+1!')
        #end if
    #end def read_text


    # test needed for crystal
    def write_text(self,format=None,occ=None):
        text = ''
        format = format.lower()
        if format=='qmcpack':
            return self.write_qmcpack()
        elif format=='casino':
            return self.write_casino()
        #end if
        channel_order = [self.local]
        for c in self.l_channels:
            if c in self.components and c!=self.local:
                channel_order.append(c)
            #end if
        #end for
        basis = self.basis
        if basis!=None:
            bs = GaussianBasisSet()
            bs.basis = basis
            basis = bs
        #end if
        if format=='gamess':
            if basis!=None:
                text += '{0} {1} 0. 0. 0.\n'.format(self.element,self.Zcore+self.Zval)
                text += basis.write_text(format)
                text += '\n'
            #end if
            text += '{0}-PP GEN {1} {2}\n'.format(self.element,self.Zcore,self.lmax)
            for c in channel_order:
                channel = self.components[c]
                text += '{0}\n'.format(len(channel)) 
                for i in sorted(channel.keys()):
                    g = channel[i]
                    text += '{0:12.8f} {1} {2:12.8f}\n'.format(g.coeff,g.rpow,g.expon)
                #end for
            #end for
            text += '\n'
        elif format=='gaussian':
            if basis!=None:
                text += '{0} 0\n'.format(self.element)
                text += basis.write_text(format)
                text += '\n'
            #end if
            text += '{0} 0\n'.format(self.element)
            text += '{0}_PP {1} {2}\n'.format(self.element,self.lmax,self.Zcore)
            for c in channel_order:
                channel = self.components[c]
                text += '{0} channel\n'.format(c)
                text += '{0}\n'.format(len(channel)) 
                for i in sorted(channel.keys()):
                    g = channel[i]
                    text += '{0} {1:12.8f} {2:12.8f}\n'.format(g.rpow,g.expon,g.coeff)
                #end for
            #end for
            text += '\n'
        elif format=='crystal':
            if basis!=None:
                conv_atomic_number = 200 + pt[self.element].atomic_number
                text+='{0} {1}\n'.format(conv_atomic_number,basis.size())
                btext = basis.write_text(format,occ=occ)
            else:
                btext = ''
            #end if
            text += 'INPUT\n'
            tline = '{0}'.format(int(self.Zval))
            channels = []
            cloc = self.components[channel_order[0]]
            if len(cloc)==1 and abs(cloc[0].coeff)<1e-8:
                tline += ' 0'
            else:
                tline += ' {0}'.format(len(cloc))
                channels.append(cloc)
            #end if
            ccount = 1
            for c in channel_order[1:]:
                channel = self.components[c]
                tline += ' {0}'.format(len(channel))
                channels.append(channel)
                ccount += 1
            #end for
            for i in range(6-ccount): # crystal14 goes up to g (hence the 6)
                tline += ' 0'
            #end for
            text += tline+'\n'
            for channel in channels:
                for i in sorted(channel.keys()):
                    g = channel[i]
                    text += '{0} {1} {2}\n'.format(g.expon,g.coeff,g.rpow-2)
                #end for
            #end for
            text += btext
        elif format=='atomscf':
            text += '{0} core potential\n'.format(self.element)
            text += '{0}\n'.format(self.Zcore)
            local_channel = self.components[self.local]
            for c in self.l_channels:
                if c in self.components:
                    channel = self.components[c]
                    if c!=self.local:
                        text += '{0}\n'.format(len(channel)+len(local_channel)) 
                    else:
                        text += '{0}\n'.format(len(channel)) 
                    #end if
                    for i in sorted(channel.keys()):
                        g = channel[i]
                        text += '{0} {1:12.8f} {2:12.8f}\n'.format(g.rpow,g.expon,g.coeff)
                    #end for
                    if c!=self.local:
                        channel = local_channel
                        for i in sorted(channel.keys()):
                            g = channel[i]
                            text += '{0} {1:12.8f} {2:12.8f}\n'.format(g.rpow,g.expon,g.coeff)
                        #end for
                    #end if
                #end if
            #end for
            text += '\n'
        elif format=='numhf':
            channel_order = self.l_channels[:self.lmax+1]
            text += '{} {}\n'.format(self.Zval,len(self.components))
            for c in channel_order:
                text += '{} '.format(len(self.components[c]))
            #end for
            text = text[:-1]+'\n'
            for c in channel_order:
                comp = self.components[c]
                for i in sorted(comp.keys()):
                    g = comp[i]
                    text += '{0} {1:12.8f} {2:12.8f}\n'.format(g.rpow,g.expon,g.coeff)
                #end for
            #end for
        else:
            self.error('ability to write file format {0} has not been implemented'.format(format))
        #end if
        return text
    #end def write_text


    # test needed
    def get_basis(self):
        bs = None
        if self.basis!=None:
            bs = GaussianBasisSet()
            bs.basis = self.basis.copy()
        #end if
        return bs
    #end def get_basis


    def set_basis(self,bs):
        self.basis = bs.basis
    #end def set_basis


    # test needed
    def uncontract(self):
        if self.basis!=None:
            bs = GaussianBasisSet()
            bs.basis = self.basis.copy()
            bs.uncontract()
            self.basis = bs.basis
        #end if
    #end def uncontract


    # test needed
    def write_basis(self,filepath=None,format=None):
        basis = self.get_basis()
        text = ''
        if basis!=None:
            if format=='gamess':
                text += '{0} {1} 0. 0. 0.\n'.format(self.element,self.Zcore+self.Zval)
                text += basis.write_text(format)
                text += '\n'
            elif format=='gaussian':
                text += '{0} 0\n'.format(self.element)
                text += basis.write_text(format)
                text += '\n'
            else:
                self.error('ability to write basis for file format {0} has not been implemented'.format(format))
            #end if
        #end if
        if filepath!=None:
            fobj = open(filepath,'w')
            fobj.write(text)
            fobj.close()
        #end if
        return text
    #end def write_basis


    def evaluate_comp_rV(self,r,l,vcomp):
        r = array(r)
        v = zeros(r.shape)
        if l==self.local or l==None:
            v += -self.Zval
        #end if
        for g in vcomp:
            if g.rpow==1:
                v += g.coeff * exp(-g.expon*r**2)
            else:
                v += g.coeff * r**(g.rpow-1) * exp(-g.expon*r**2)
            #end if
        #end for
        return v
    #end def evaluate_comp_rV


    # test needed
    def ppconvert(self,outfile,ref,extra=None):
        of = outfile.lower()
        if of.endswith('.xml'):
            opts = '--xml'
        elif of.endswith('.upf'):
            opts = '--log_grid --upf'
        else:
            self.error('output file format unrecognized for {0}\nvalid extensions are .xml and .upf'.format(outfile))
        #end if
        tmpfile = 'tmp.gamess'
        self.write(tmpfile,'gamess')
        if extra is not None:
            command = 'ppconvert --gamess_pot {0} --s_ref "{1}" --p_ref "{1}" --d_ref "{1}" {2} {3} {4}'.format(tmpfile,ref,extra,opts,outfile)
        else:
            command = 'ppconvert --gamess_pot {0} --s_ref "{1}" --p_ref "{1}" --d_ref "{1}" {2} {3}'.format(tmpfile,ref,opts,outfile)
        execute(command,verbose=True)
        os.system('rm '+tmpfile)
    #end def ppconvert


    # test needed
    def append_to_component(self,l,coeff,expon,rpow):
        '''
        This function is used to append a term to a Gaussian ECP component.
        l: the angular ccomponent that the Gaussian term will be appended to 
        coeff, expon, rpow: the coefficient, exponent, and r-power of the Gaussian term
        '''
        if l>self.lmax:
            self.error('component {} not present in PP.'.format(l))
        #end if
        chan_labels = ['s','p','d','f','g','h','i','j']
        self.components[chan_labels[l]].append(obj(coeff=coeff,expon=expon,rpow=rpow))
    #end def append_to_component


    # test needed
    def scale_component(self,l,scale):
        '''
        This function is used to scale a Gaussian ECP component by a factor.
        l: the angular component that is scaled.
        scale: the scaling factor
        '''
        if l>self.lmax:
            self.error('component {} not present in PP.'.format(l))
        #end if
        chan_labels = ['s','p','d','f','g','h','i','j']
        for term in self.components[chan_labels[l]]:
            term.coeff*=scale
        #end for
    #end def scale_component


    # test needed
    def simplify(self):
        '''
        This function simplifies the Guassian ECP. The simplificactions are as follows:
        1. Remove all terms with coefficients that are equal to zero -- unless only one term exists.
        2. Within each component, look for terms that have matching exponents and r-powers, if any are
           present, then sum their coefficicents to make a single term. If the coefficients sum to
           zero, then remove the terms (unless that leaves the component empty)
        '''
        dec=16
        # Remove terms with coefficients equivalent to zero
        #chan_labels = ['s','p','d','f','g','h','i','k']
        chan_labels = list(self.l_channels)
        remove = []
        for l in np.arange(self.lmax+1):
            for term_idx,term in enumerate(self.components[chan_labels[l]]):
                if abs(term.coeff)<1e-12 and len(self.components[chan_labels[l]])>1:
                    remove.append((chan_labels[l],term_idx))
                #end if
            #end for
        #end for
        for r in remove:
            self.components[r[0]].delete(r[1])
        #end for
        comps = self.components.copy()
        for l in np.arange(self.lmax+1):
            comps[chan_labels[l]] = obj()
            for term_idx,term in enumerate(self.components[chan_labels[l]]):
                comps[chan_labels[l]].append(term)
            #end for
        #end for
        self.components = comps.copy()
        comps = self.components.copy()
        for l in np.arange(self.lmax+1):
            terms = []
            comps[chan_labels[l]] = obj()
            for term_idx,term in enumerate(self.components[chan_labels[l]]):
                terms.append(term.list())
            #end for
            terms = np.array(terms)
            like_terms = []
            if len(terms)>1:
                rpows  = terms[:,2]
                expons = terms[:,1]
                coeffs = terms[:,0]
                for ex_idx,ex in enumerate(expons.round(decimals=dec)):
                    if any(ex_idx in subl for subl in like_terms):
                        continue
                    #end if
                    match = np.argwhere(expons.round(decimals=dec)==ex)
                    if len(match)>1:
                        match = match.flatten()
                        unique_pows = np.unique(rpows[match].round(decimals=dec))
                        if len(unique_pows)==1:
                            like_terms.append(match.tolist())
                        else:
                            for uv in unique_pows:
                                uv_count = rpows[match].round(decimals=dec).tolist().count(uv)
                                if uv_count>1:
                                    m = match[rpows[match].round(decimals=dec)==uv][0]
                                    if any(m in subl for subl in like_terms):
                                        continue
                                    else:
                                        like_terms.append(match[rpows[match].round(decimals=dec)==uv].tolist())
                                    #end if
                                #end if
                            #end for
                        #end if
                    #end if
                #end for
            #end if
            # update comps
            comps[chan_labels[l]] = obj()
            added = []
            for term_idx,term in enumerate(self.components[chan_labels[l]]):
                if any(term_idx in subl for subl in like_terms):
                    if term_idx in added:
                        continue
                    else:
                        for mlist in like_terms:
                            if term_idx in mlist and not term_idx in added:
                                coeff = 0.0
                                mod_term = term.copy()
                                for ti in mlist: 
                                    coeff += self.components[chan_labels[l]][ti].coeff
                                #end for
                                if abs(coeff)>1e-12:
                                    mod_term.coeff = coeff
                                    comps[chan_labels[l]].append(mod_term)
                                #end if
                                added.extend(mlist)
                            #end if
                        #end for
                    #end if
                else:
                    comps[chan_labels[l]].append(term)
                    added.append(term_idx)
                #end if
            #end for
            if len(comps[chan_labels[l]])==0:
                # All terms cancelled. Add placeholder
                plcehldr = self.components['s'][0].copy()
                plcehldr.coeff = 0.0
                plcehldr.rpow = 2
                plcehldr.expon = 1.0
                comps[chan_labels[l]].append(plcehldr)
        #end for
        self.components = comps.copy()
    #end def simplify


    # test needed
    def is_truncated_L2(self):
        '''
        Determine if the Gaussian ECP's channels follow an L2 relationship.
        '''
        # CHECK IF THIS WORKS FOR lmax=1 !!!!!!
        # Only checked for lmax=2 and higher
        p1 = self.copy()
        p1.simplify()
        p2 = self.copy()
        p2.transform_to_truncated_L2(keep='s p',lmax=p2.lmax)
        p2.simplify()
        return object_eq(p2,p1)
    #end def is_truncated_L2


    # test needed
    def get_unboundedness(self,db,dbs):
        '''
        This function quantifies how unbounded a truncated L2 potential is.
        This is done by constructing a function that corrects VL2 in the unbounded region.
        Then integrating the difference between VL2 and the correcting function.
        '''
        if not self.is_truncated_L2():
            self.error('The PP must be in the truncated L2 form.')
        #end if
        import math
        def poly(x,c):
            val=0
            for ci,cv in enumerate(c):
                val+=cv*x**ci
            return val
        #end def
        
        def Rs(x,dx,s,c):
            if x+1-s<-dx:
                return 0-(1-s)
            elif x+1-s>dx:
                return x
            else:
                return poly(x+1-s,c)-(1-s)
        #end def
        A=[]
        for i in range(8):
            row=[]
            for j in range(8):
                if i<4:
                    if j-i<0:
                        dcoeff=0
                        dpower=0
                    else:
                        dcoeff=math.factorial(j)/math.factorial(j-i)
                        dpower=j-i
                    #end if
                    row.append(dcoeff*db**dpower)
                else:
                    if j-(i-4)<0:
                        dcoeff=0
                        dpower=0
                    else:
                        dcoeff=math.factorial(j)/math.factorial(j-(i-4))
                        dpower=j-(i-4)
                    #end if
                    row.append(dcoeff*(-db)**dpower)
                #end if
            #end for
            A.append(row)
        #end for
        
        A = np.array(A)
        b = np.array([db,1]+[0]*6)
        c = np.linalg.inv(A).dot(b)

        ng=3000
        gmin=0.02
        gmax=0.85
        r = np.linspace( gmin, gmax, ng )


        # PP
        v = []
        for l in ['s','p']:
            vtmp = []
            for ri in r:
                vtmp.append(self.evaluate_component(r=ri,l=l))
            #for
            v.append(vtmp)
        #for
        v=np.array(v)

        # 2*r^2*VL2 
        if self.lmax>1:
            f = r*r*(v[1]-v[0])
        elif self.lmax==1:
            f = -r*r*v[0]
        else:
            self.error('Not sure what to do with fully local potential.')
        #end if
            
        # 2*r^2*V'L2 
        fp = [Rs(fr,db,dbs,c) for fr in f]

        unboundedness = 0
        for fi,fx in enumerate(f):
            unboundedness+=(fp[fi]-fx)*(gmax-gmin)/ng 
        #end for

        return unboundedness
    #end def get_unboundedness


    # test needed
    def make_L2_bounded(self,db,dbs,exps0=None,plot=False):
        '''
        For a truncated L2 potential, this function constructs a correction to VL2 in the unbounded region.
        Then the correction is fit to a set of Gaussian primitives that are provided in the array 'exps0'.
        The fitted Gaussian primitives are then appended to the ECP, resulting in a bounded truncated L2 potential.
        '''
        if not self.is_truncated_L2():
            self.error('The PP must be in the truncated L2 form.')
        #end if
        if exps0 is None:
            self.error('Please provide a aet of exponents to be used for correction.')
        import math
        def poly(x,c):
            val=0
            for ci,cv in enumerate(c):
                val+=cv*x**ci
            return val
        #end def
        
        def Rs(x,dx,s,c):
            if x+1-s<-dx:
                return 0-(1-s)
            elif x+1-s>dx:
                return x
            else:
                return poly(x+1-s,c)-(1-s)
        #end def
        class fitClass:
        
            def __init__(self):
                pass
        
            def gauss_correction(self,x,c1,c2,c3):
                val = 0
                for ci,c in enumerate([c1,c2,c3]):
                    val+=x**2.*c*np.exp(-self.exps[ci]*x**2.)
                #end for
                return val
            #end def
        
            def gauss_correction_2_param(self,x,c1,c2):
                val = 0
                for ci,c in enumerate([c1,c2]):
                    val+=x**2.*c*np.exp(-self.exps[ci]*x**2.)
                #end for
                return val
            #end def
        
            def gauss_correction_1_param(self,x,c1):
                val = 0
                for ci,c in enumerate([c1]):
                    val+=x**2.*c*np.exp(-self.exps[ci]*x**2.)
                #end for
                return val
            #end def
        
        #end class

        A=[]
        for i in range(8):
            row=[]
            for j in range(8):
                if i<4:
                    if j-i<0:
                        dcoeff=0
                        dpower=0
                    else:
                        dcoeff=math.factorial(j)/math.factorial(j-i)
                        dpower=j-i
                    #end if
                    row.append(dcoeff*db**dpower)
                else:
                    if j-(i-4)<0:
                        dcoeff=0
                        dpower=0
                    else:
                        dcoeff=math.factorial(j)/math.factorial(j-(i-4))
                        dpower=j-(i-4)
                    #end if
                    row.append(dcoeff*(-db)**dpower)
                #end if
            #end for
            A.append(row)
        #end for
        
        A = np.array(A)
        b = np.array([db,1]+[0]*6)
        c = np.linalg.inv(A).dot(b)

        ng=3000
        gmin=0.02
        gmax=0.85
        r = np.linspace( gmin, gmax, ng )


        # PP
        v = []
        for l in ['s','p']:
            vtmp = []
            for ri in r:
                vtmp.append(self.evaluate_component(r=ri,l=l))
            #for
            v.append(vtmp)
        #for
        v=np.array(v)

        # 2*r^2*VL2 
        if self.lmax>1:
            f = r*r*(v[1]-v[0])
        elif self.lmax==1:
            f = -r*r*v[0]
        else:
            self.error('Not sure what to do with fully local potential.')
        #end if
            
        # 2*r^2*V'L2 
        fp = [Rs(fr,db,dbs,c) for fr in f]

        unboundedness = 0
        for fi,fx in enumerate(f):
            unboundedness+=(fp[fi]-fx)*(gmax-gmin)/ng 
        #end for
        #print('\npseudopotential undoundedness: ',undoundedness)
        from scipy.optimize import curve_fit

        fit_instance = fitClass()
        fit_instance.exps=exps0
        if len(exps0)==1:
            popt, pcov = curve_fit(fit_instance.gauss_correction_1_param, r, f-fp)
        elif len(exps0)==2:
            popt, pcov = curve_fit(fit_instance.gauss_correction_2_param, r, f-fp)
        elif len(exps0)==3:
            popt, pcov = curve_fit(fit_instance.gauss_correction, r, f-fp)
        else:
            self.error('Number of correction primitives not coded.')
        #end if

        if plot:
            import matplotlib.pyplot as plt

            plt.plot(r, fp, 'g-', label='fp')

            plt.xlabel('r (bohr)')
            plt.ylabel('$2r^2v_{L^2}$')

            if len(exps0)==1:
                plt.plot(r, f-fit_instance.gauss_correction_1_param(r, *popt), 'r-',label='f-corr')
            elif len(exps0)==2:
                plt.plot(r, f-fit_instance.gauss_correction_2_param(r, *popt), 'r-',label='f-corr')
            elif len(exps0)==3:
                plt.plot(r, f-fit_instance.gauss_correction(r, *popt), 'r-',label='f-corr')
            #end if

            plt.plot(r, f, 'b-',label='f')
            plt.plot(r, [-1]*len(r), 'k-',label=None)
            plt.legend()
            plt.show()
        #end if
        for expon_idx,expon in enumerate(exps0):
            self.components['s'].append(obj(coeff=1.0*popt[expon_idx],expon=expon,rpow=2))
        #end if
        self.transform_to_truncated_L2(keep='s p',lmax=self.lmax)
        self.simplify()

    #end def make_L2_bounded


    # test needed
    def transform_to_truncated_L2(self,keep=None,lmax=None,outfile=None,inplace=True):
        '''
        This function transforms a Gaussian ECP into a truncated L2 form, i.e., a form
        for which all channels follow an L2 relationship. For a semi-local ECP, this 
        transformation can have a significant negative impact on transferability. For
        an ECP that is already in a trucnated L2 form, the transformation has no affect.
        '''
        ##############################################################################
        # WARNING: ONLY PERFORM THIS TRANSFORMATION IF YOU KNOW WHAT YOU ARE DOING.
        #          TRANSFERABILITY IS GENERALLY REDUCED SEVERELY AFTER TRANSFORM.
        ##############################################################################
        comps = list(self.components.keys())
        if keep is None or lmax is None:
            self.error('parameters \'keep\' and \'lmax\' must be specified.')
        #end if
        chan_labels = ['s','p','d','f','g','h','i','j']
        keep_chans = keep.split()
        # Are the labels recognized?
        if keep_chans[0] not in chan_labels or keep_chans[1] not in chan_labels:
            slef.error('Requested channel to keep is not recognized')
        #end if
        # Does the original potential contain the requested channels?
        if keep_chans[0] not in comps or keep_chans[1] not in comps:
            self.error('Cannot keep channel that is not already present')
        #end if
        ## Are the requested 'keep' channels different?
        if chan_labels.index(keep_chans[0]) == chan_labels.index(keep_chans[1]):
            self.error('The two channels must be different.')
        #end if
        keep_l_vals = []
        keep_l_vals.append(chan_labels.index(keep_chans[0]))
        keep_l_vals.append(chan_labels.index(keep_chans[1]))
        keep_l_vals.sort()
        # Is one of the channels the local channel?
        if keep_l_vals[0]==self.lmax or keep_l_vals[1]==self.lmax:
            keep_local=True
        else:
            keep_local=False
        #end if
        old_lmax = self.lmax
        old_local = self.local
        self.lmax  = lmax
        self.local = chan_labels[lmax]
        if not keep_local:
            
            lm = keep_l_vals[0]
            ln = keep_l_vals[1]

            self.components[chan_labels[lmax]] = self.components[old_local]
            fctr = lm*(lm+1)-lmax*(lmax+1)
            fctr = float(fctr)/(lm*(lm+1)-ln*(ln+1))
            for term in self.components[chan_labels[ln]]:
                self.append_to_component(lmax,fctr*term.coeff,term.expon,term.rpow)
            #end for
            fctr = lmax*(lmax+1)-ln*(ln+1)
            fctr = float(fctr)/(lm*(lm+1)-ln*(ln+1))
            for term in self.components[chan_labels[lm]]:
                self.append_to_component(lmax,fctr*term.coeff,term.expon,term.rpow)
            #end for

            vm_comp = self.components[chan_labels[lm]].copy()
            vn_comp = self.components[chan_labels[ln]].copy()
            for l in np.arange(lmax):
                fctr = l*(l+1)-lmax*(lmax+1)
                fctr = float(fctr)/(lm*(lm+1)-ln*(ln+1))
                self.components[chan_labels[l]] = obj()
                for term_idx,term in enumerate(vm_comp):
                    self.append_to_component(l,coeff=fctr*term.coeff,expon=term.expon,rpow=term.rpow)
                #end for
                for term_idx,term in enumerate(vn_comp):
                    self.append_to_component(l,coeff=-fctr*term.coeff,expon=term.expon,rpow=term.rpow)
                #end for
            #end for

        else:
            
            lloc = keep_l_vals[1]
            lm = keep_l_vals[0]

            fctr = lmax*(lmax+1)-lloc*(lloc+1)
            fctr = float(fctr)/(lm*(lm+1)-lloc*(lloc+1))
            self.components[chan_labels[lmax]] = self.components[chan_labels[lloc]]
            for term in self.components[chan_labels[lm]]:
                self.append_to_component(lmax,fctr*term.coeff,term.expon,term.rpow)
            #end for

            vm_comp = self.components[chan_labels[lm]].copy()
            for l in np.arange(lmax):
                fctr = l*(l+1)-lmax*(lmax+1)
                fctr = float(fctr)/(lm*(lm+1)-lloc*(lloc+1))
                self.components[chan_labels[l]] = obj()
                for term_idx,term in enumerate(vm_comp):
                    self.append_to_component(l,coeff=fctr*term.coeff,expon=term.expon,rpow=term.rpow)
                #end for
            #end for

        #end if

        self.simplify()

    #end def transform_to_truncated_L2
#end class GaussianPP





class QmcpackPP(SemilocalPP):
    requires_format = False
    numeric         = True
    interpolatable  = False

    def read(self,filepath,format=None):
        if not os.path.exists(filepath):
            self.error('cannot read {0}, file does not exist'.format(filepath))
        #end if
        
        x = readxml(filepath,contract_names=True)
        x.convert_numeric()
        x.condense()
        x.remove_hidden()
        pp = x.pseudo
        
        h = pp.header
        self.element = h.symbol
        self.Zval    = h.zval
        self.Zcore   = h.atomic_number-h.zval
        if self.Zcore==0:
            self.core = '0'
        else:
            self.core = pt.simple_elements[self.Zcore].symbol
        #end if

        g = pp.grid
        if g.type=='linear':
            self.rmin = g.ri
            self.rmax = g.rf
            self.r = linspace(g.ri,g.rf,g.npts)
        else:
            self.error('functionality for '+g.type+' grids has not yet been implemented')
        #end if
        if 'l2' in pp:
            l2 = pp.l2
            if l2.format!='r*V':
                self.error('unrecognized potential format: {0}\nthe only supported format is r*V'.format(l2.format))
            #end if
            self.components.L2 = l2.radfunc.data.copy()
        #end if
        sl = pp.semilocal
        if sl.format!='r*V':
            self.error('unrecognized potential format: {0}\nthe only supported format is r*V'.format(sl.format))
        #end if
        lloc = self.l_channels[sl.l_local]
        self.local = lloc
        vps = sl.vps
        if not isinstance(vps,list):
            vps = [vps]
        #end if
        self.lmax = len(vps)-1
        for vp in vps:
            self.components[vp.l] = vp.radfunc.data.copy()
        #end for
        for l in self.angular_channels():
            if l!=self.local:
                self.components[l] -= self.components[self.local]
            #end if
        #end for
    #end def read


    def evaluate_comp_rV(self,r,l,vcomp):
        if r is not None:
            if len(r)==len(self.r) and abs( (r[1:]-self.r[1:])/self.r[1:] ).max()<1e-6:
                r = self.r
            else:
                self.error('ability to interpolate at arbitrary r has not been implemented\ncalling evaluate_channel() without specifying r will return the potential on a default grid')
            #end if
        else:
            r = self.r
        #end if
        v = vcomp.copy()
        return v
    #end def evaluate_comp_rV


    def v_at_zero(self,l):
        #r = self.r
        #v = self.get_component(l)/r
        #vz = (v[1]*r[2]**2-v[2]*r[1]**2)/(r[2]**2-r[1]**2)
        r = self.r[1:3]
        v = self.get_component(l)[1:3]/r
        vz = (v[0]*r[1]**2-v[1]*r[0]**2)/(r[1]**2-r[0]**2)
        return vz
    #end def v_at_zero
#end class QmcpackPP




class CasinoPP(SemilocalPP):
    requires_format = False
    numeric         = True
    interpolatable  = False

    unitmap = dict(rydberg='Ry',hartree='Ha',ev='eV')

    def read(self,filepath,format=None):
        if not os.path.exists(filepath):
            self.error('cannot read {0}, file does not exist'.format(filepath))
        #end if
        # open the file
        file = TextFile(filepath)
        # read scalar values at the top
        Zatom,Z = file.readtokensf('Atomic number and pseudo-charge',int,float)
        if Zatom>len(pt.simple_elements):
            self.error('element {0} is not in the periodic table')
        #end if
        element = pt.simple_elements[Zatom].symbol
        units = file.readtokensf('Energy units',str)
        if not units in self.unitmap:
            self.error('units {0} unrecognized from casino PP file {1}'.format(units,filepath))
        #end if
        lloc = file.readtokensf('Angular momentum of local component',int)
        lloc = self.l_channels[lloc]
        ngrid = file.readtokensf('Number of grid points',int)
        # read the radial grid
        file.seek('R(i)',1)
        file.readline()
        r = empty((ngrid,),dtype=float)
        for ir in range(ngrid):
            r[ir] = float(file.readline())
        #end for
        # read each channel, convert to hartree units
        lmax  = -1
        lvals = []
        lpots = []
        while(file.seek('pot',1)!=-1):
            lmax += 1
            potline = file.readline() # read the r*potential line
            eqloc = potline.find('=')
            if eqloc==-1:
                self.error('"=" not found in potential line\nline: {0}'.format(potline))
            #end if
            l = self.l_channels[int(potline[eqloc+1])] # get the l value
            lvals.append(l)
            v = empty((ngrid,),dtype=float)
            for ir in range(ngrid):
                v[ir] = float(file.readline())
            #end for
            lpots.append(convert(v,self.unitmap[units],'Ha'))
        #end while

        # fill in SemilocalPP class data obtained from read
        self.element = element
        self.Zval    = int(Z)
        self.Zcore   = Zatom-self.Zval
        if self.Zcore==0:
            self.core = '0'
        else:
            self.core = pt.simple_elements[self.Zcore].symbol
        #end if
        self.lmax  = lmax
        self.local = lloc
        self.r     = r
        for i in range(len(lpots)):
            self.components[lvals[i]] = lpots[i]
        #end for
        for l in self.angular_channels():
            if l!=self.local:
                self.components[l] -= self.components[self.local]
            #end if
        #end for
    #end def read_file


    def evaluate_comp_rV(self,r,l,vcomp):
        if r is not None:
            if len(r)==len(self.r) and abs( (r[1:]-self.r[1:])/self.r[1:] ).max()<1e-6:
                r = self.r
            else:
                self.error('ability to interpolate at arbitrary r has not been implemented\ncalling evaluate_channel() without specifying r will return the potential on a default grid')
            #end if
        else:
            r = self.r
        #end if
        v = vcomp.copy()
        return v
    #end def evaluate_comp_rV

#end class CasinoPP

