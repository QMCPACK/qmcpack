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


#! /usr/bin/env python

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
try:
    from scipy.optimize import curve_fit
except:
    curve_fit = unavailable('curve_fit')
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
            #elem_label = self.filename.split('.')[0]
            #is_elem,symbol = is_element(elem_label,symbol=True)
            if not is_elem:
                self.error('cannot determine element for pseudopotential file: {0}\npseudopotential file names must be prefixed by an atomic symbol or label\n(e.g. Si, Si1, etc)'.format(filepath))
            #end if
            self.element = symbol
            self.element_label = elem_label
            #self.element  = self.filename[0:2].strip('._-')
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
        print
        print '  Pseudopotentials'
        for filepath in ppfiles:
            print '    reading pp: ',filepath
            ext = filepath.split('.')[-1].lower()
            if ext=='gms':
                pp = gamessPPFile(filepath)
            else:
                pp = PseudoFile(filepath)
            #end if
            pps.append(pp)
        #end for
        print
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
        for code,pps in code_pps.iteritems():
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
    all_channels = ['loc']+list(l_channels)
    channel_colors = obj(s='g',p='r',d='b',f='m',g='c',h='k')

    numeric        = False
    interpolatable = True

    formats = ['qmcpack']

    channel_indices = obj()
    for i,c in enumerate(l_channels):
        channel_indices[c] = i
    #end for

    def __init__(self,filepath=None,format=None,name=None,src=None):
        self.name = name
        self.rcut = None
        self.lmax = None
        self.local = None
        if self.numeric:
            self.r = None
        #end if
        self.channels = obj()
        Pseudopotential.__init__(self,filepath,format)
        if src!=None:
            self.transfer_core_from(src)
        #end if
    #end def __init__

    
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


    def get_channel(self,l=None):
        if l is None:
            l = self.local
        #end if
        if not l in self.channels:
            self.error('cannot get invalid channel: {0}\n  valid options are: {1}'.format(l,self.channels.keys()))
        #end if
        return self.channels[l]
    #end def get_channel


    def remove_channel(self,l):
        if l in self.channels:
            del self.channels[l]
            if self.local==l:
                lmax = -1
                for lt in self.channels.keys():
                    li = self.l_int_from_text(lt)
                    if li>lmax:
                        lmax = li
                        self.local = lt
                    #end if
                #end for
                self.lmax = lmax
            #end if
        #end if
    #end def remove_channel


    def l_int_from_text(self,ltext):
        lint = 'spdfghi'.find(ltext)
        return lint
    #end def l_int_from_text
        

    def evaluate(self,r=None,l=None,rpow=0,with_local=False):
        if self.numeric and not self.interpolatable:
            r = None
        elif r is None and self.interpolatable:
            r = linspace(0.01,4.0,400)            
        #end if
        if not with_local:
            v = self.evaluate_rV(r,l)
        elif l==None or l==self.local:
            v = self.evaluate_rV(r,l)
        else:
            v = self.evaluate_rV(r,l)+self.evaluate_rV(r,self.local)
        #end if
        if self.numeric and not self.interpolatable:
            r = self.r
        #end if
        if rpow!=1:
            v = r**(rpow-1)*v
        #end if
        return v
    #end def evaluate


    def evaluate_rV(self,r=None,l=None): # returns r*V
        self.not_implemented()
    #end def evaluate_rV


    def numeric_channel(self,l=None,rmin=0.,rmax=10.,npts=10001,rpow=0,with_local=False):
        if self.numeric and not self.interpolatable:
            v = self.evaluate(None,l,rpow,with_local)
            r = self.r
        else:
            r = linspace(rmin,rmax,npts)
            v = self.evaluate(r,l,rpow,with_local)
        #end if
        return r,v
    #end def numeric_channel


    def update_rcut(self,tol=1e-5):
        self.rcut = self.find_rcut(tol=tol)
        return self.rcut
    #end def update_rcut


    def find_rcut(self,tol=1e-5):
        r    = None
        vmin = None
        vmax = None
        for l in self.channels.keys():
            rc,vc = self.numeric_channel(l,rmin=0.01,with_local=True)
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
        for i in xrange(nr):
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
        color_in = color
        if channels is None:
            channels = self.all_channels
        #end if
        if fig:
            figure()
        #end if
        if r is None and self.numeric:
            r = self.r
        elif r is None:
            r = linspace(rmin,rmax,1000)
        #end if
        for c in channels:
            if c in self.channels: 
                if c==self.local:
                    lab = self.local+' loc'
                    color = self.channel_colors[self.local]
                else:
                    lab = c
                    color = self.channel_colors[c]
                #end if
                if color_in!=None:
                    color = color_in
                #end if
                if self.name!=None:
                    lab = self.name+' '+lab
                #end if
                v = self.evaluate(r,c,with_local=with_local)
                if metric=='r2':
                    v = r**2*v
                    if c==self.local:
                        v += self.Zval*r
                    #end if
                elif metric!=None:
                    self.error('invalid metric for plotting: {0}\nvalid options are: r2'.format(metric))
                #end if
                #l = self.channel_indices[c]
                #if l==0:
                #    lmult = 1
                #else:
                #    lmult = (l*(l+1))
                #vs = self.evaluate(r,'s',with_local=with_local)
                #v-=vs
                #v/=lmult
                plot(r,v,color+linestyle,label=lab)
            #end for
        #end for
        if fig:
            if title is None:
                title = 'Semilocal {0} PP ({1} core)'.format(self.element,self.core)
            #end if
            set_title(title)
            ylabel('channels')
            xlabel('r')
            legend()
        #end if
        if show:
            show_plots()
        #end if
    #end def plot

                
    def plot_L2(self,show=True,fig=True,r=None,rmin=0.01,rmax=5.0,linestyle='-',title=None):
        if r is None and self.numeric:
            r = self.r
        elif r is None:
            r = linspace(rmin,rmax,1000)
        #end if
        vs = self.evaluate(r,'s',with_local=True)
        channels = self.all_channels
        for c in channels[2:]:
            if c in self.channels:
                color = self.channel_colors[c]
                v = self.evaluate(r,c,with_local=True)
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


    def gaussian_fit(self,r=None,pmax=3,maxfev=100000,verbose=False,filepath=None,format=None,offset=0):
        if verbose:
            self.log('\nfitting {0} pseudopotential to gaussians'.format(self.element))
        #end if
        gf.Z = self.Zval
        if r is None:
            r = self.r
        #end if
        #r2 = r**2
        channels = obj()
        for l in self.channels.keys():
            #channels[l] = self.evaluate(r,l)[offset:]*r2[offset:]
            channels[l] = self.evaluate(r,l,rpow=2)
        #end for
        #r = r[offset:]
        #del r2
        p0nl  = 10.0,0.3
        p0loc = 0.2,0.3,10.,0.3
        padd  = [0.0,0.3,.27]
        channel_params = obj()
        for l in self.channels.keys():
            if l==self.local:
                channel_params[l] = tuple(p0loc)
            else:
                channel_params[l] = tuple(p0nl)
            #end if
        #end for
        pmt = (1,pmax)
        if not pmt in gf.functions or not pmt in gf.loc_functions:
            self.error('cannot include {0} paired functions\nplease choose a smaller number for pmax'.format(pmax))
        #end if
        gfits = GFit()
        gfits.proto = GaussianPP()
        gfits.proto.transfer_core_from(self)
        channel_fits = obj()
        for l in self.channels.keys():
            if verbose:
                self.log('  fitting channel {0}'.format(l))
            #end if
            lfits = obj()
            fail = False
            for p in range(pmax+1):
                if verbose:
                    self.log('    performing fit of order {0}'.format(p))
                #end if
                try:

                    if l==self.local:
                        gf_fit = gf.loc_functions[1,p]
                        Zval   = self.Zval
                    else:
                        gf_fit = gf.functions[1,p]
                        Zval   = None
                    #end if
                    v  = channels[l]
                    p0 = channel_params[l]
                    vp,vc = curve_fit(gf_fit,r,v,p0,maxfev=maxfev)
                    gc = get_gf_channel(1,p,vp,Zval)
                    vf = gf_fit(r,*vp)
                    rms_dev = gf.rms_deviation(v,vf)
                    energy_dev = gf.energy_deviation(v,vf,r)
                    cfit = obj(
                        channel    = gc,
                        rms_dev    = rms_dev,
                        energy_dev = energy_dev
                        )
                    if verbose:
                        self.log('      rms    deviation: {0}'.format(rms_dev))
                        self.log('      energy deviation: {0}'.format(energy_dev))
                    #end if
                    lfits.append(cfit)

                except StandardError,e:
                    #print e
                    fail = True
                #end try
                if fail:
                    if verbose:
                        self.log('    order {0} fit failed, skipping higher orders'.format(p))
                    #end if
                    break
                #end if
                channel_params[l] = tuple(list(channel_params[l])+padd)
            #end for
            channel_fits[l] = lfits
        #end for
        gfits.channels = channel_fits
        if filepath!=None:
            for p,gpp in gfits.iteritems():
                fpath = filepath.format(p)
                if verbose:
                    self.log('  writing '+fpath)
                #end if
                gpp.write(fpath,format)
            #end for
        #end if
        if verbose:
            self.log('  fitting complete')
        #end if
        return gfits
    #end def gaussian_fit


    def write_qmcpack(self,filepath=None):
        if self.rcut is None:
            self.update_rcut(tol=1e-5)
        #end if

        symbol        = self.element
        atomic_number = self.Zcore+self.Zval
        zval          = self.Zval
        creator       = 'Nexus'
        npots_down    = len(self.channels)
        l_local       = 'spdfgi'.find(self.local)

        rmin = 1e99
        rmax = -1e99
        npts = 0
        vps = obj()
        for l in self.channels.keys():
            r,v = self.numeric_channel(l,rpow=1,with_local=True)
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
                semilocal+='        </data>\n'
                semilocal+='      </radfunc>\n'
                semilocal+='    </vps>\n'
            #end if
        #end for
        semilocal+='  </semilocal>\n'
        footer = '</pseudo>\n'

        text = header+grid+semilocal+footer

        if filepath!=None:
            open(filepath,'w').write(text)
        #end if
        return text
    #end def write_qmcpack


    def write_casino(self,filepath=None):
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
        for l in self.channels.keys():
            r,v = self.numeric_channel(l,rpow=1,with_local=True)
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



class GFit(DevBase):
    def __init__(self):
        self.proto    = None # SemilocalPP
        self.channels = None
    #end def __init__


    def best_channels(self):
        cfits = obj()
        for l,pfits in self.channels.iteritems():
            emin = 1e99
            cmin = None
            for p,cfit in pfits.iteritems():
                if abs(cfit.energy_dev)<emin:
                    emin = abs(cfit.energy_dev)
                    cmin = cfit
                #end if
            #end for
            cfits[l] = cfit
        #end for
        return cfits
    #end def best_channels

    def best(self):
        pp = self.proto.copy()
        for l,cfit in self.best_channels().iteritems():
            pp.channels[l] = cfit.channel
        #end for
        return pp
    #end def best


    def report(self,fits='best'):
        print
        print 'Gaussian fits for '+self.proto.element
        if fits=='best':
            cfits = self.best_channels()
            for l in self.proto.l_channels:
                if l in cfits:
                    cfit = cfits[l]
                    print '  {0}  {1}  {2}'.format(l,cfit.energy_dev,cfit.rms_dev)
                #end if
            #end for
        else:
            self.error('fits option {0} is not supported'.format(fits))
        #end if
    #end def report
#end class GFit





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
            self.channels[cname] = channel
        #end for
        self.basis = basis
        if len(self.channels)!=self.lmax+1:
            self.error('number of channels is not lmax+1!')
        #end if
    #end def read_text


    def write_text(self,format=None,occ=None):
        text = ''
        format = format.lower()
        if format=='qmcpack':
            return self.write_qmcpack()
        #end if
        channel_order = [self.local]
        for c in self.all_channels:
            if c in self.channels and c!=self.local:
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
                channel = self.channels[c]
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
                channel = self.channels[c]
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
            cloc = self.channels[channel_order[0]]
            if len(cloc)==1 and abs(cloc[0].coeff)<1e-8:
                tline += ' 0'
            else:
                tline += ' {0}'.format(len(cloc))
                channels.append(cloc)
            #end if
            ccount = 1
            for c in channel_order[1:]:
                channel = self.channels[c]
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
            local_channel = self.channels[self.local]
            for c in self.all_channels:
                if c in self.channels:
                    channel = self.channels[c]
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


    def uncontract(self):
        if self.basis!=None:
            bs = GaussianBasisSet()
            bs.basis = self.basis.copy()
            bs.uncontract()
            self.basis = bs.basis
        #end if
    #end def uncontract


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


    def evaluate_rV(self,r,l=None):
        r = array(r)
        v = zeros(r.shape)
        if l==self.local or l==None:
            v += -self.Zval
        #end if
        for g in self.get_channel(l):
            if g.rpow==1:
                v += g.coeff * exp(-g.expon*r**2)
            else:
                v += g.coeff * r**(g.rpow-1) * exp(-g.expon*r**2)
            #end if
        #end for
        return v
    #end def evaluate_rV


    def ppconvert(self,outfile,ref):
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
            self.channels[vp.l] = vp.radfunc.data.copy()
        #end for
        for l in self.channels.keys():
            if l!=self.local:
                self.channels[l] -= self.channels[self.local]
            #end if
        #end for
    #end def read


    def evaluate_rV(self,r=None,l=None):
        if r!=None:
            if len(r)==len(self.r) and abs( (r[1:]-self.r[1:])/self.r[1:] ).max()<1e-6:
                r = self.r
            else:
                self.error('ability to interpolate at arbitrary r has not been implemented\ncalling evaluate() without specifying r will return the potential on a default grid')
            #end if
        else:
            r = self.r
        #end if
        v = self.get_channel(l)
        return v
    #end def evaluate_rV


    def v_at_zero(self,l):
        #r = self.r
        #v = self.get_channel(l)/r
        #vz = (v[1]*r[2]**2-v[2]*r[1]**2)/(r[2]**2-r[1]**2)
        r = self.r[1:3]
        v = self.get_channel(l)[1:3]/r
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
        for ir in xrange(ngrid):
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
            for ir in xrange(ngrid):
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
            self.channels[lvals[i]] = lpots[i]
        #end for
        for l in self.channels.keys():
            if l!=self.local:
                self.channels[l] -= self.channels[self.local]
            #end if
        #end for
    #end def read_file
#end class CasinoPP







# functions for gaussing fitting
def get_gf_channel(ns,np,params,Zval=None):
    loc = Zval!=None
    g = obj(
        coeff = [],
        rpow  = [],
        expon = []
        )
    nparam = 2*ns+3*np
    if loc:
        nparam += 2
    #end if
    if nparam!=len(params):
        raise RuntimeError('wrong number of parameters given to get_gf_coeff_rpow_expon\nparameters expected: {0}\nparameters given: {1}'.format(nparam,len(params)))
    #end if
    pcur = 0
    if loc:
        a,b = params[pcur:pcur+2]
        g.coeff += [ -Zval*0, Zval    , .5*Zval/a**2 ]
        g.rpow  += [ -1   , -1      ,  1          ]
        g.expon += [  0   , .5/a**2 , .5/b**2     ]
        pcur+=2
    #end if
    for s in range(ns):
        b,bs = params[pcur:pcur+2]
        g.coeff.append(b/bs)
        g.rpow.append( 0)
        g.expon.append(.5/bs**2) 
        pcur+=2
    #end for
    for p in range(np):
        b,bs,bt = params[pcur:pcur+3]
        pf = 2*b/(bs+bt)
        g.coeff += [ pf/bs    , -pf/bt   ]
        g.rpow  += [  1       ,  1       ]
        g.expon += [ .5/bs**2 , .5/bt**2 ]
        pcur+=3
    #end for
    g.coeff = array(g.coeff)
    g.rpow  = array(g.rpow, dtype=int)
    g.expon = array(g.expon)
    g.rpow += 2
    channel = obj()
    for t in range(len(g.coeff)):
        channel[t] = obj(
            coeff = g.coeff[t],
            rpow  = g.rpow[t],
            expon = g.expon[t]
            )
    return channel
#end def get_gf_channel

def gf_s1_p0(r,b1,bs1):
    r2 = r**2/2 
    return 2*r2*(b1/bs1*exp(-r2/bs1**2))
#end def gf_s1_p0

def gf_s1_p1(r,b1,bs1,b2,bs2,bt2):
    r2 = r**2/2
    return 2*r2*(  b1/bs1*exp(-r2/bs1**2) \
                 + b2*(2*r/(bs2+bt2))*(exp(-r2/bs2**2)/bs2-exp(-r2/bt2**2)/bt2) )
#end def gf_s1_p1

def gf_s1_p2(r,b1,bs1,b2,bs2,bt2,b3,bs3,bt3):
    r2 = r**2/2
    return 2*r2*(  b1/bs1*exp(-r2/bs1**2) \
                 + b2*(2*r/(bs2+bt2))*(exp(-r2/bs2**2)/bs2-exp(-r2/bt2**2)/bt2) \
                 + b3*(2*r/(bs3+bt3))*(exp(-r2/bs3**2)/bs3-exp(-r2/bt3**2)/bt3) )
#end def gf_s1_p2

def gf_s1_p3(r,b1,bs1,b2,bs2,bt2,b3,bs3,bt3,b4,bs4,bt4):
    r2 = r**2/2
    return 2*r2*(  b1/bs1*exp(-r2/bs1**2) \
                 + b2*(2*r/(bs2+bt2))*(exp(-r2/bs2**2)/bs2-exp(-r2/bt2**2)/bt2) \
                 + b3*(2*r/(bs3+bt3))*(exp(-r2/bs3**2)/bs3-exp(-r2/bt3**2)/bt3) \
                 + b4*(2*r/(bs4+bt4))*(exp(-r2/bs4**2)/bs4-exp(-r2/bt4**2)/bt4) )
#end def gf_s1_p3



def gf_loc_s0_p0(r,a,b):
    r2 = r**2/2
    return -gf.Z*r*( 1.0 -exp(-r2/a**2) -r2/a**2*exp(-r2/b**2) )
#end def gf_loc_s0_p0

def gf_loc_s1_p0(r,a,b,b1,bs1):
    r2 = r**2/2
    return -gf.Z*r*( 1.0 -exp(-r2/a**2) -r2/a**2*exp(-r2/b**2) ) +2*r2*(b1/bs1*exp(-r2/bs1**2))
#end def gf_loc_s1_p0

def gf_loc_s1_p1(r,a,b,b1,bs1,b2,bs2,bt2):
    r2 = r**2/2
    return -gf.Z*r*( 1.0 -exp(-r2/a**2) -r2/a**2*exp(-r2/b**2) ) \
           +2*r2*(  b1/bs1*exp(-r2/bs1**2) \
                  + b2*(2*r/(bs2+bt2))*(exp(-r2/bs2**2)/bs2-exp(-r2/bt2**2)/bt2) )
#end def gf_loc_s1_p1

def gf_loc_s1_p2(r,a,b,b1,bs1,b2,bs2,bt2,b3,bs3,bt3):
    r2 = r**2/2
    return -gf.Z*r*( 1.0 -exp(-r2/a**2) -r2/a**2*exp(-r2/b**2) ) \
           +2*r2*(  b1/bs1*exp(-r2/bs1**2) \
                  + b2*(2*r/(bs2+bt2))*(exp(-r2/bs2**2)/bs2-exp(-r2/bt2**2)/bt2) \
                  + b3*(2*r/(bs3+bt3))*(exp(-r2/bs3**2)/bs3-exp(-r2/bt3**2)/bt3) )
#end def gf_loc_s1_p2

def gf_loc_s1_p3(r,a,b,b1,bs1,b2,bs2,bt2,b3,bs3,bt3,b4,bs4,bt4):
    r2 = r**2/2
    return -gf.Z*r*( 1.0 -exp(-r2/a**2) -r2/a**2*exp(-r2/b**2) ) \
           +2*r2*(  b1/bs1*exp(-r2/bs1**2) \
                  + b2*(2*r/(bs2+bt2))*(exp(-r2/bs2**2)/bs2-exp(-r2/bt2**2)/bt2) \
                  + b3*(2*r/(bs3+bt3))*(exp(-r2/bs3**2)/bs3-exp(-r2/bt3**2)/bt3) \
                  + b4*(2*r/(bs4+bt4))*(exp(-r2/bs4**2)/bs4-exp(-r2/bt4**2)/bt4) )
#end def gf_loc_s1_p3

def gf_rms_deviation(v1,v2):
    return sqrt(((v1-v2)**2).mean())
#end def gf_rms_deviation

def gf_energy_deviation(v1,v2,r):
    return (v1-v2).sum()*(r[1]-r[0])
#end def gf_energy_deviation


gf = obj(
    Z = None,
    functions = {
        (1,0):gf_s1_p0,
        (1,1):gf_s1_p1,
        (1,2):gf_s1_p2,
        (1,3):gf_s1_p3
        },
    loc_functions = {
        (0,0):gf_loc_s0_p0,
        (1,0):gf_loc_s1_p0,
        (1,1):gf_loc_s1_p1,
        (1,2):gf_loc_s1_p2,
        (1,3):gf_loc_s1_p3
        },
    rms_deviation = gf_rms_deviation,
    energy_deviation = gf_energy_deviation
    )











# older classes, kept in case contents are useful later

#class PseudoFile(DevBase):
#    colors = dict(s='k',p='r',d='b',f='m')
#    lcolors= ['k','r','b','m']
#    ldict = {0:'s',1:'p',2:'d',3:'f'}
#
#    format = 'generic'
#
#    conv_table = {('upf','fsatom'):('--upf_pot','--xml')}
#    extensions = dict(upf='upf',fsatom='xml',gamess='gamess',casino='data')
#
#
#    standard_energy_units   = 'eV'
#    standard_distance_units = 'B'
#
#    def readfile(self,filepath):
#        self.not_implemented()
#    #end def read
#
#
#    def write(self,filepath):
#        self.not_implemented()
#    #end def write
#
#
#    def __init__(self,filepath=None,energy_units=None):
#        self.initialize = False
#        if energy_units is None:
#            self.energy_units = self.standard_energy_units
#        else:
#            self.energy_units = energy_units
#        #end if
#        self.set(
#            filename         = None,
#            location         = None,
#            element          = None,
#            type             = None,
#            Z                = None,
#            r                = None,
#            potentials       = None,
#            rcut             = None,
#            potential_spread = None,
#            pp               = None
#            )
#        if filepath!=None:
#            self.read(filepath)
#        #end if
#    #end def __init__
#
#
#    def read(self,filepath):
#        self.filename = os.path.basename(filepath)
#        self.location = os.path.abspath(filepath)
#        self.readfile(filepath)
#        if self.initialize:
#            self.rcut = self.get_rcut(1e-16)
#            pot = array(self.potentials.values())
#            self.potential_spread = pot.max(0)-pot.min(0)
#        #end if
#    #end def read
#
#
#    def convert(self,pptype):
#        oform = pptype.format
#        forms = self.format,oform
#        have_conv = forms in self.conv_table 
#        if not have_conv:
#            self.error('pseudopotential conversion from '+self.format+' to '+pptype.format+' format has not yet been implemented')
#        #end if
#        sflags,oflags = self.conv_table[forms]
#        sfilepath = self.location
#        path,sfilename = os.path.split(sfilepath)
#        if not oform in self.extensions:
#            self.error('attempted to convert to an unknown format: '+oform)
#        #end if
#        ofilename = sfilename.rsplit('.',1)[0]+'.'+self.extensions[oform]
#        ofilepath = os.path.join(path,ofilename)
#        command = 'ppconvert {0} {1} {2} {3}\n'.format(sflags,sfilename,oflags,ofilename)
#        cwd = os.getcwd()
#        os.chdir(path)
#        ppcin = open(os.path.join(path,'ppconvert.in'),'w')
#        ppcin.write(command)
#        ppcin.close()
#        out = open('./ppconvert.out','w')
#        err = open('./ppconvert.err','w')
#        p = Popen(command,stdout=out,stderr=err,shell=True)
#        p.communicate()
#        out.close()
#        err.close()
#        os.chdir(cwd)
#        pp = pptype(ofilepath)
#        return pp
#    #end def convert
#
#
#    def get_rcut(self,Etol=1e-16,rlower=.5):
#        v = self.potentials.values()
#        rcut = 0
#        ir = self.r>rlower
#        r = self.r[ir]
#        for i in range(1,len(v)):
#            vd = abs(v[i][ir]-v[0][ir])
#            rc = r[abs(vd-Etol).argmin()]
#            rcut = max(rcut,rc)
#        #end for
#        return rcut
#    #end def get_rcut
#
#        
#    def plot(self,style='',rmax=4,vcoul=False,ptitle=None,mult=1,units='eV',lw=2):
#        colors = self.lcolors
#        r = self.r
#        rc= self.rcut
#        Z = self.Z
#        if mult=='r':
#            mult = r
#        #end if
#        vmin = 1e55
#        vmax = -1e55
#        for l,vpot in self.potentials.iteritems():
#            pot = convert(mult*vpot,self.energy_units,units)
#            vmin = min(vmin,pot.min())
#            vmax = max(vmax,pot.max())
#            plot(r,pot,colors[l]+style,lw=lw,label=self.ldict[l])
#        #end for
#        if vcoul:
#            vcoul = convert(-Z/r,'Ha',units)
#            vin = vcoul>vmin
#            plot(r[vin],mult*vcoul[vin],'k-.')
#        #dnd if
#        plot([rc,rc],[vmin,vmax],'k-.',lw=lw)
#        if rmax!=None:
#            xlim([0,rmax])
#        #end if
#        if ptitle is None:
#            title('Channel potentials for '+self.element+' '+self.type+' pseudopotential')
#        else:
#            title(ptitle)
#        #end if
#        ylabel('Potential Energy ({0})'.format(units))
#        xlabel('Radius (bohr)')
#        legend()
#    #end def plot
#
#        
#    def plot_spread(self,style='k',label='',rmax=5):
#        r = self.r
#        rc= self.rcut
#        vs = self.potential_spread
#        semilogy(r,vs,style,lw=2,label=label)
#        semilogy([rc,rc],[vs.min(),vs.max()],'k-.',lw=2)
#        if rmax!=None:
#            xlim([0,rmax])
#        #end if
#        grid()
#        title("Potential spread $(\max_{\ell m}\, |v_\ell(r)-v_m(r)|)$ for "+self.element+' '+self.type+' pseudopotential')
#        ylabel('Potential spread ({0})',self.energy_units)
#        xlabel('Radius (bohr)')
#        if label!='':
#            legend()
#        #end if
#    #end def plot_spread
##end class PseudoFile
#
#
#
#class fsatomPP(PseudoFile):
#    format = 'fsatom'
#    def readfile(self,filepath):
#        None
#        #x = readxml(filepath,contract_names=True)
#        #x.convert_numeric()
#        #x.condense()
#        #x.remove_hidden()
#        #pp = x.pseudo
#        #self.pp = pp
#        #
#        #h = pp.header
#        #self.element = h.symbol
#        #self.type = h.flavor        
#        #self.Z    = h.zval
#        #vps = self.pp.semilocal.vps
#        #if not isinstance(vps,list):
#        #    vps = [vps]
#        ##end if
#        #g = vps[0].radfunc.grid
#        #if g.type=='linear':
#        #    r = linspace(g.ri,g.rf,g.npts)
#        #    self.r = r[1:]
#        #else:
#        #    self.error('functionality for '+g.type+' grids has not yet been implemented')
#        ##end if
#        #p = obj()
#        #r = self.r
#        #ldict = dict(s=0,p=1,d=2,f=3)
#        #for vp in vps:
#        #    l = vp.l
#        #    v = 1./r*vp.radfunc.data[1:]
#        #    p[ldict[l]]= convert(v,'Ha',self.energy_units)
#        ##end for
#        #self.potentials = p
#    #end def readfile
##end class fsatomPP
#
#
#class ncppPP(PseudoFile):
#    format = 'ncpp'
#    def readfile(self,filepath):
#        text = open(filepath,'r').read()
#        
#        lines = text.splitlines()
#        
#        functional = lines[0].split()[0].strip("'")
#        es,zs = lines[1].split(',')[0:2]
#        element = es.strip("'")
#        Zval  = int(float(zs)+.5)
#        
#        self.set(
#            element = element,
#            type    = functional,
#            Z       = Zval,
#            r          = [],
#            potentials = obj(),
#            pp         = obj(),
#            initialize = False
#            )
#    #end def readfile
##end class ncppPP
#
#
#class upfPP(PseudoFile):
#    format='upf'
#    def readfile(self,filepath):
#        text = open(filepath,'r').read()
#        if '<UPF' not in text:
#            upf_format = 'old'
#            lines = text.split('\n')
#            xml = '<upf>\n'
#            for l in lines:
#                if l.find('/>')!=-1:
#                    ln = l.replace('/>','>').replace('<','</')
#                else:
#                    ln = l
#                #end if
#                xml+=ln+'\n'
#            #end for
#            xml += '</upf>\n'
#            tmppath = filepath+'_tmp'
#            open(tmppath,'w').write(xml)
#            x = readxml(tmppath,contract_names=True,strip_prefix='pp_')
#            os.system('rm '+tmppath)
#            x.convert_numeric()
#            x.condense()
#            x.remove_hidden()
#            pp = x.upf
#        else:
#            upf_format = 'new'
#            #x = readxml(filepath,contract_names=True,strip_prefix='pp_')
#            #x.convert_numeric()
#            #x.condense()
#            #x.remove_hidden()
#            #pp = x.upf
#            pp = obj()
#            pp_contents = open(filepath,'r').read()
#        #end if
#        if upf_format=='old':
#            lines = pp.header.split('\n')
#            h = obj()
#            i=0
#            h.version = string2val(lines[i].split()[0]); i+=1
#            h.element = string2val(lines[i].split()[0]); i+=1
#            h.type = string2val(lines[i].split()[0]); i+=1
#            ncc = string2val(lines[i].split()[0]); i+=1
#            h.nonlinear_cc = dict(T=True,F=False)[ncc]
#            h.functional = lines[i].split()[0:4]; i+=1
#            h.Z = string2val(lines[i].split()[0]); i+=1
#            h.total_energy = string2val(lines[i].split()[0]); i+=1
#            h.wfc_cutoff,h.rho_cutoff = array(lines[i].split()[0:2],dtype=float); i+=1
#            h.lmax = string2val(lines[i].split()[0]); i+=1
#            h.npts = string2val(lines[i].split()[0]); i+=1
#            h.nwfc = string2val(lines[i].split()[0]); i+=1
#            pp.header = h
#            if 'nonlocal' in pp:
#                if 'beta' in pp.nonlocal:
#                    beta = pp.nonlocal.beta
#                    if isinstance(beta,str):
#                        beta = [beta]
#                    #end if
#                    b = obj()
#                    for i in range(len(beta)):
#                        sections = beta[i].split('\n',2)
#                        p = string2val(sections[0].split()[0])
#                        v = string2val(sections[2])
#                        b[p]=v
#                    #end for
#                    pp.nonlocal.beta = b
#                #end if
#                dij = pp.nonlocal.dij
#                d = obj()
#                lines = dij.split('\n')[1:]
#                for l in lines:
#                    t = l.split()
#                    d[int(t[0]),int(t[1])] = string2val(t[2])
#                #end for
#                pp.nonlocal.dij = d
#            #end if
#            pswfc = pp.pswfc
#            tokens = pswfc.split()
#            nwfc= pp.header.nwfc
#            npts= pp.header.npts
#            wf = []
#            for n in range(nwfc):
#                i = n*(4+npts)
#                label = tokens[i]
#                l = int(tokens[i+1])
#                ws = []
#                for v in tokens[i+4:i+4+npts]:
#                    ws.append(float(v))
#                #end for
#                orb = obj()
#                orb.label = label
#                orb.l = l
#                orb.wfc = array(ws)
#                wf.append(orb)
#            #end for
#            pp.pswfc = wf
#        
#            #fill in standard fields
#            self.pp = pp
#            self.r = pp.mesh.r[1:]
#            if 'local' in pp:
#                self.local = convert(pp.local,'Ry',self.energy_units)[1:]
#                if 'nonlocal' in pp:
#                    nl = obj()
#                    vnl = zeros(self.local.shape)
#                    if 'beta' in pp.nonlocal:
#                        beta = pp.nonlocal.beta
#                        for t,d in pp.nonlocal.dij.iteritems():
#                            bi = beta[t[0]]
#                            bj = beta[t[1]]
#                            if not isinstance(bi,str) and not isinstance(bj,str):
#                                bb = d*bi*bj
#                            else: # the file is being misread, fix later
#                                bb  = 0*pp.mesh.r
#                            #end if
#                            naftcut = len(pp.mesh.r)-len(bb)
#                            if naftcut>0:
#                                bb=append(bb,zeros((naftcut,)))
#                            #end if
#                            vnl += bb[1:]/self.r**2
#                        #end for
#                    #end if
#                    vnl = convert(vnl,'Ry',self.energy_units)
#                    nl[0] = vnl
#                    self.nonlocal = nl
#                    h = pp.header
#                    p = obj()
#                    p[0] = self.local
#                    p[1] = self.local+self.nonlocal[0]
#                    self.potentials = p
#                #end if
#            #end if
#            self.element = h.element
#            self.type = h.type
#            self.Z  = h.Z
#        else:
#            header_start = pp_contents.find('<PP_HEADER')
#            if header_start==-1:
#                self.error('could not find <PP_HEADER>')
#            #end if
#            header_end   = pp_contents.find('/>')
#            if header_end==-1:
#                self.error('could not find </PP_HEADER>')
#            #end if
#            header = pp_contents[header_start:header_end]
#            tokens = header.split()[1:]
#            for token in tokens:
#                if '=' in token:
#                    name,value = token.split('=',1)
#                    value = value.strip('"')
#                    if name=='element':
#                        self.element = value
#                    elif name=='functional':
#                        self.type = value
#                    elif name=='z_valence':
#                        self.Z = int(round(float(value)))
#                    #end if
#                #end if
#            #end for
#        
#            #self.error('ability to read new UPF format has not yet been implemented\n  attempted to read '+filepath)
#            #self.warn('ability to read new UPF format has not yet been implemented\n  attempted to read '+filepath)
#            self.initialize = False
#        #end if
#    #end def readfile
##end class upfPP
#
#
#class gamessPP(PseudoFile):
#    format = 'gamess'
#    def readfile(self,filepath):
#        text = open(filepath,'r').read()
#        lines = text.splitlines()
#        pp = obj()
#        i=0
#        pp.name,pp.type,Zcore,lmax = lines[i].split(); i+=1
#        Zcore = int(Zcore)
#        lmax  = int(lmax)
#        pp.Zcore = Zcore
#        pp.lmax  = lmax
#        
#        element = split_delims(pp.name)[0]
#        if not element in pt:
#            element = split_delims(self.filename)[0]
#            if not element in pt:
#                self.error('cannot identify element for pseudopotential file '+filepath)
#            #end if
#        #end if
#        Zatom = pt[element].atomic_number
#        Z = Zatom-Zcore
#        
#        r = mgrid[1.e-10:150.00001:.005]
#        p = obj()
#        vlocal = None
#        for index in range(lmax+1):
#            l = (index+lmax)%(lmax+1)
#            ngpot = int(lines[i].strip()); i+=1
#            coeffs    = empty((ngpot,),dtype=float)
#            powers    = empty((ngpot,),dtype=int)
#            exponents = empty((ngpot,),dtype=float)
#            for ig in range(ngpot):
#                coef,power,exponent = lines[i].split(); i+=1
#                coeffs[ig]    = float(coef)
#                powers[ig]    = int(power) 
#                exponents[ig] = float(exponent)
#            #end for
#            pp[index] = obj(coeffs=coeffs,powers=powers,exponents=exponents)
#            v = 0*r
#            for ig in range(ngpot):
#                v += coeffs[ig]*r**(powers[ig]-2)*exp(-exponents[ig]*r**2)
#            #end for
#            if index==0:
#                vlocal = v - Z/r
#                p[l] = vlocal.copy()
#            else:
#                p[l] = v + vlocal
#            #end if
#        #end for
#        for l in p.keys():
#            p[l] = convert(p[l],'Ha',self.energy_units)
#        #end for
#        
#        self.set(
#            element    = element,
#            type       = pp.type,
#            Z          = Z,
#            r          = r,
#            potentials = p
#            )
#    #end def readfile
##end class gamessPP
#
#
#class TextFile(DevBase):
#    def __init__(self,filepath):
#        if not os.path.exists(filepath):
#            self.error('text file '+filepath+' does not exist')
#        #end if
#        text = open(filepath,'r').read()
#        lines= text.splitlines()
#        self.text = text
#        self.lines = lines
#        self.filepath = filepath
#    #end def __init__
#
#    def find_line(self,text,exit=False):
#        for i in range(len(self.lines)):
#            if text in self.lines[i]:
#                return i
#            #end if
#        #end for
#        if exit:
#            self.error('text "'+text+'" not found in file '+self.filepath)
#        else:
#            return None
#        #end if
#    #end def find_line
#
#    def read_tokens(self,iline,*formats):
#        if isinstance(iline,str):
#            iline = self.find_line(iline,exit=True)+1
#        #end if
#        tokens = []
#        stokens = self.lines[iline].split()
#        if len(formats)==1 and len(stokens)>1:
#            formats = len(stokens)*formats
#        elif len(formats)>len(stokens):
#            self.error('line {0} only has {1} tokens, you requested {2}'.format(iline,len(stokens),len(formats)))
#        #end if
#        for i in range(len(formats)):
#            tokens.append(formats[i](stokens[i]))
#        #end for
#        if len(tokens)==1:
#            return tokens[0]
#        else:
#            return tokens
#        #end if
#    #end def read_tokens
##end def TextFile
#
#
#class casinoPP(PseudoFile):
#    format = 'casino'
#    unitmap = dict(rydberg='Ry',hartree='Ha',ev='eV')
#    def readfile(self,filepath):
#        text = TextFile(filepath)
#        Zatom,Z = text.read_tokens('Atomic number and pseudo-charge',int,float)
#        if Zatom>len(pt.simple_elements):
#            self.error('element {0} is not in the periodic table')
#        #end if
#        element = pt.simple_elements[Zatom].symbol
#        units = text.read_tokens('Energy units',str)
#        if not units in self.unitmap:
#            self.error('units {0} unrecognized from casino PP file {1}'.format(units,filepath))
#        #end if
#        lloc = text.read_tokens('Angular momentum of local component',int)
#        ngrid = text.read_tokens('Number of grid points',int)
#        i = text.find_line('R(i)',exit=True)+1
#        r = empty((ngrid,),dtype=float)
#        for ir in xrange(ngrid):
#            r[ir] = float(text.lines[i])
#            i+=1
#        #end for
#        r=r[1:]
#        p = obj()
#        while i<len(text.lines):
#            line = text.lines[i]
#            if 'potential' in line:
#                eqloc = line.find('=')
#                if eqloc==-1:
#                    self.error('"=" not found in potential line')
#                #end if
#                l = int(line[eqloc+1])
#                i+=1
#                if i+ngrid>len(text.lines):
#                    self.error('potentials in file {0} are not the right length'.format(filepath))
#                #end if
#                v = empty((ngrid,),dtype=float)
#                for ir in xrange(ngrid):
#                    v[ir] = float(text.lines[i])
#                    i+=1
#                #end for
#                p[l] = v[1:]/r
#            #end if
#        #end while
#        
#        for l in p.keys():
#            p[l] = convert(p[l],self.unitmap[units],self.energy_units)
#        #end for
#        
#        self.set(
#            element = element,
#            type = 'Trail-Needs',
#            Z = Z,
#            r = r,
#            potentials = p,
#            pp = obj(
#                Zatom = Zatom,
#                Z     = Z,
#                units = units,
#                lloc  = lloc,
#                ngrid = ngrid
#                )
#            )
#    #end def read_file
##end class casinoPP
#
#
#class Pseudopotentials(DevBase):
#    def __init__(self,*pseudopotentials):
#        if len(pseudopotentials)==1 and isinstance(pseudopotentials[0],list):
#            pseudopotentials = pseudopotentials[0]
#        #end if
#        ppfiles = []
#        pps     = []
#        errors = False
#        print
#        for pp in pseudopotentials:
#            if isinstance(pp,PseudoFile):
#                pps.append(pp)
#            elif isinstance(pp,str):
#                ppfiles.append(pp)
#            else:
#                self.error('expected PseudoFile type or filepath, got '+str(type(pp)),exit=False)
#                errors = True
#            #end if
#        #end for
#        if errors:
#            self.error('cannot create Pseudopotentials object')
#        #end if
#
#        if len(pps)>0:
#            self.addpp(pps)
#        #end if
#        if len(ppfiles)>0:
#            self.readpp(ppfiles)
#        #end if
#    #end def __init__
#
#
#    def addpp(self,*pseudopotentials):
#        if len(pseudopotentials)==1 and isinstance(pseudopotentials[0],list):
#            pseudopotentials = pseudopotentials[0]
#        #end if
#        for pp in pseudopotentials:
#            self[pp.filename] = pp
#        #end for
#    #end def addpp
#
#        
#    def readpp(self,*ppfiles):
#        if len(ppfiles)==1 and isinstance(ppfiles[0],list):
#            ppfiles = ppfiles[0]
#        #end if
#        pps = []
#        errors = False
#        print '  Pseudopotentials'
#        for filepath in ppfiles:
#            print '    reading pp: ',filepath
#            ext = filepath.split('.')[-1].lower()
#            if ext=='upf':
#                pp = upfPP(filepath)
#            elif ext=='xml':
#                pp = fsatomPP(filepath)
#            elif ext=='ncpp':
#                pp = ncppPP(filepath)
#            elif ext=='gms':
#                pp = gamessPP(filepath)
#            elif ext=='data':
#                pp = casinoPP(filepath)
#            else:
#                self.error('cannot determine pseudopotential type from file extension '+ext+' ('+os.path.basename(filepath)+')')
#                errors = True
#                pp = None
#            #end if
#            pps.append(pp)
#        #end for
#        if errors:
#            self.error('cannot read all pseudopotentials')
#        #end if
#        self.addpp(pps)
#    #end def readpp
#
#
#    def read_type(self,pptype,*ppfiles):
#        if len(ppfiles)==1 and isinstance(ppfiles[0],list):
#            ppfiles = ppfiles[0]
#        #end if
#        pps = []
#        for filepath in ppfiles:
#            pp = pptype(filepath)
#            pps.append(pp)
#        #end for
#        self.addpp(pps)
#    #end def read_type
#    
#    
#    def read_upf(self,*ppfiles):
#        self.read_type(upfPP,*ppfiles)
#    #end def read_upf
#    
#    
#    def read_fsatom(self,*ppfiles):
#        self.read_type(fsatomPP,*ppfiles)
#    #end def read_fsatom
#
##end class Pseudopotentials

