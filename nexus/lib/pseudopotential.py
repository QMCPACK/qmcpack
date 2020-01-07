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
    l_channels   = tuple('spdfghi')
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
        self.set_component('L2',guard=True)
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
                plot(r,v,color+linestyle,label=lab)
            #end for
        #end for
        if fig:
            if title is None:
                title = 'Semilocal {0} PP ({1} core)'.format(self.element,self.core)
            #end if
            set_title(title)
            ylabel('component potentials (Ha)')
            xlabel('r (Bohr)')
            legend()
        #end if
        if show:
            show_plots()
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
                    v = self.evaluate_channel(r,c,rpow,rmin-1e-12,with_local,with_L2)
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
                v = self.evaluate_channel(r,c,rpow,rmin-1e-12,with_local,with_L2)
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

                
    def plot_L2(self,show=True,fig=True,r=None,rmin=0.01,rmax=5.0,linestyle='-',title=None):
        if r is None and self.numeric:
            r = self.r
        elif r is None:
            r = linspace(rmin,rmax,1000)
        #end if
        vs = self.evaluate_channel(r,'s',with_local=True,rmin=rmin-1e-12)
        for c in self.l_channels[1:]:
            if c in self.components:
                color = self.channel_colors[c]
                v = self.evaluate_channel(r,c,with_L2=False,rmin=rmin-1e-12)
                rng = r>rmin-1e-12
                r = r[rng]
                l = self.channel_indices[c]
                vL2 = (v-vs)/(l*(l+1))
                plot(r,vL2,color+linestyle,label='(v{0}-vs)/(l(l+1))'.format(c))
            #end if
        #end for
        if fig:
            xlim([0,rmax])
            if title is None:
                title = 'Semilocal {0} PP ({1} core)'.format(self.element,self.core)
            #end if
            set_title(title)
            ylabel('vL2 for channels above s')
            xlabel('r')
            legend()
            if show:
                show_plots()
            #end if
        #end if 
    #end def plot_L2


    def write_qmcpack(self,filepath=None):
        self.update_rcut(tol=1e-5,optional=True)

        channels = self.angular_channels()

        symbol        = self.element
        atomic_number = self.Zcore+self.Zval
        zval          = self.Zval
        creator       = 'Nexus'
        npots_down    = len(channels)
        l_local       = 'spdfgi'.find(self.local)

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
    formats = SemilocalPP.formats + 'gaussian gamess crystal'.split()

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

