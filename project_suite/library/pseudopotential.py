#! /usr/bin/env python

import os
from subprocess import Popen
from execute import execute
from numpy import linspace,array,zeros,append,mgrid,empty,exp
from xmlreader import readxml
from superstring import string2val,split_delims
from periodic_table import pt
from unit_converter import convert
from generic import obj
from developer import DevBase,unavailable
from debug import *
try:
    from matplotlib.pyplot import figure,plot,xlabel,ylabel,title,show,ylim,legend,xlim,rcParams,savefig,bar,xticks,subplot,grid,setp,errorbar,loglog,semilogx,semilogy

    params = {'legend.fontsize':14,'figure.facecolor':'white','figure.subplot.hspace':0.,
          'axes.labelsize':16,'xtick.labelsize':14,'ytick.labelsize':14}
    rcParams.update(params)
except (ImportError,RuntimeError):
   figure,plot,xlabel,ylabel,title,show,ylim,legend,xlim,rcParams,savefig,bar,xticks,subplot,grid,setp,errorbar,loglog,semilogx,semilogy = unavailable('matplotlib.pyplot','figure','plot','xlabel','ylabel','title','show','ylim','legend','xlim','rcParams','savefig','bar','xticks','subplot','grid','setp','errorbar','loglog','semilogx','semilogy')
#end try


# basic interface for project suite, only gamess really needs this for now
class PseudoFile(DevBase):
    def __init__(self,filepath=None):
        self.element  = None
        self.filename = None
        self.location = None
        if filepath!=None:
            self.filename = os.path.basename(filepath)
            self.location = os.path.abspath(filepath)
            self.element  = self.filename[0:2].strip('._-')
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
        print
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
        self.addpp(pps)
    #end def readpp


    def pseudos_by_atom(self,*ppfiles):
        pps = obj()
        for ppfile in ppfiles:
            if ppfile in self:
                pp = self[ppfile]
                pps[pp.element] = pp
            #end if
        #end for
        return pps
    #end def pseudos_by_atom
#end class Pseudopotentials















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
        self.Zeff    = None
        self.Zcore   = None
        if filepath!=None:
            self.read(filepath,format)
        #end if
    #end def __init__


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
        self.read_text(text,format)
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


    def read_text(self,text,format=None):
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
    l_channels   = tuple('spdfg')
    all_channels = ['loc']+list(l_channels)
    channel_colors = obj(loc='k',s='g',p='r',d='b',f='m',g='c')

    def __init__(self,filepath=None,format=None,name=None):
        self.name = name
        self.lmax = None
        self.local = None
        self.channels = obj()
        Pseudopotential.__init__(self,filepath,format)
    #end def __init__


    def read(self,filepath,format=None):
        Pseudopotential.read(self,filepath,format)
        lc = self.l_channels
        c = self.channels
        self.local = list(set(lc[0:len(c)])-set(c.keys()))[0]
    #end def read


    def get_channel(self,l):
        if l==self.local:
            l = 'loc'
        #end if
        if not l in self.channels:
            self.error('cannot get invalid channel: {0}\n  valid options are: {1}'.format(l,self.channels.keys()))
        #end if
        return self.channels[l]
    #end def get_channel

        
    def evaluate(self,r=None,l='loc'):
        self.not_implemented()
    #end def evaluate


    def evaluate_with_local(self,r=None,l='loc'):
        if l=='loc':
            return self.evaluate(r,l)
        else:
            return self.evaluate(r,l)+self.evaluate(r,'loc')
        #end if
    #end def evaluate_with_local


    def rcut(self,tol=1e-5):
        self.not_implemented()
    #end def rcut


    def plot(self,r=None,show=True,fig=True,linestyle='-',channels=None,with_local=False,rmin=0.01,rmax=5.0,title=None):
        if channels is None:
            channels = self.all_channels
        #end if
        if fig:
            figure()
        #end if
        for c in channels:
            if c in self.channels:
                if c=='loc':
                    lab = self.local+' loc'
                    color = self.channel_colors[self.local]
                else:
                    lab = c
                    color = self.channel_colors[c]
                #end if
                if self.name!=None:
                    lab = self.name+' '+lab
                #end if
                if with_local:                    
                    v = self.evaluate_with_local(r,c)
                else:
                    v = self.evaluate(r,c)
                #end if
                if r is None and 'r' in self:
                    plot(self.r,v,color+linestyle,label=lab)
                else:
                    plot(r,v,color+linestyle,label=lab)
                #end if
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
            if show:
                show_plots()
            #end if
        #end if
    #end def plot
#end class SemilocalPP

 

class GaussianPP(SemilocalPP):
    requires_format = True
    formats = 'gaussian gamess'.split()


    def read_text(self,text,format=None):
        rawlines = text.splitlines()
        lines = []
        for rline in rawlines:
            line = rline.strip()
            if (not line.startswith('!')) and (not line.startswith('#')) and len(line)>0:
                lines.append(line)
            #end if
        #end for

        format=format.lower()
        channels = []
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
        else:
            self.error('ability to read file format {0} has not been implemented'.format(format))
        #end if
        if not element in pt:
            if not self.element in pt:
                self.error('cannot identify element for pseudopotential file '+filepath)
            #end if
        else:
            self.element = element
        #end if
        Zatom = pt[element].atomic_number
        Zeff = Zatom-Zcore
        core = pt.simple_elements[Zcore].symbol
        self.set(
            core    = core,
            Zeff    = Zeff,
            Zcore   = Zcore,
            lmax    = lmax
            )
        for c in range(len(channels)):
            if c==0:
                cname = 'loc'
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
        if len(self.channels)!=self.lmax+1:
            self.error('number of channels is not lmax+1!')
        #end if
    #end def read_text


    def write_text(self,format=None):
        text = ''
        format = format.lower()
        if format=='gamess':
            text += '{0}-PP GEN {1} {2}\n'.format(self.element,self.Zcore,self.lmax)
            for c in self.all_channels:
                if c in self.channels:
                    channel = self.channels[c]
                    text += '{0}\n'.format(len(channel)) 
                    for i in sorted(channel.keys()):
                        g = channel[i]
                        text += '{0:12.8f} {1} {2:12.8f}\n'.format(g.coeff,g.rpow,g.expon)
                    #end for
                #end if
            #end for
            text += '\n'
        elif format=='gaussian':
            text += '{0} 0\n'.format(self.element)
            text += '{0}_PP {1} {2}\n'.format(self.element,self.lmax,self.Zcore)
            for c in self.all_channels:
                if c in self.channels:
                    channel = self.channels[c]
                    text += '{0} channel\n'.format(c)
                    text += '{0}\n'.format(len(channel)) 
                    for i in sorted(channel.keys()):
                        g = channel[i]
                        text += '{0} {1:12.8f} {2:12.8f}\n'.format(g.rpow,g.expon,g.coeff)
                    #end for
                #end if
            #end for
            text += '\n'
        else:
            self.error('ability to write file format {0} has not been implemented'.format(format))
        #end if
        return text
    #end def write_text


    def evaluate(self,r=None,l='loc'):
        if not l in self.channels:
            self.error('channel {0} is not present\nvalid options are: {1}'.format(self.channels.keys()))
        #end if
        if r is None:
            r = linspace(0.01,4.0,400)
        #end if
        r = array(r)
        v = zeros(r.shape)
        if l=='loc':
            v += -self.Zeff/r
        #end if
        for g in self.channels[l]:
            v += g.coeff * r**(g.rpow-2) * exp(-g.expon*r**2)
        #end for
        return v
    #end def evaluate


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

    def __init__(self,filepath=None,format=None,name=None):
        self.r = None
        SemilocalPP.__init__(self,filepath,format,name)
    #end def __init__


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
        self.Zeff    = h.zval
        self.Zcore   = h.atomic_number-h.zval
        self.core    = pt.simple_elements[self.Zcore].symbol

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
        vps = sl.vps
        if not isinstance(vps,list):
            vps = [vps]
        #end if
        for vp in vps:
            l = vp.l
            if l==lloc:
                l='loc'
            #end if
            self.channels[l] = vp.radfunc.data.copy()
        #end for
        for l in self.channels.keys():
            if l!='loc':
                self.channels[l] -= self.channels.loc
            #end if
        #end for
        lc = self.l_channels
        c = self.channels
        self.local = list(set(lc[0:len(c)])-set(c.keys()))[0]
    #end def read


    def evaluate(self,r=None,l='loc'):
        if r!=None:
            self.error('ability to interpolate at arbitrary r has not been implemented\ncalling evaluate() without specifying r will return the potential on a default grid')
        else:
            r = self.r
        #end if
        c = self.get_channel(l)
        v = c/r
        return v
    #end def evaluate


    def v_at_zero(self,l):
        r = self.r
        v = self.get_channel(l)/r
        vz = (v[1]*r[2]**2-v[2]*r[1]**2)/(r[2]**2-r[1]**2)
        return vz
    #end def v_at_zero
#end class QmcpackPP





















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
#        Zeff  = int(float(zs)+.5)
#        
#        self.set(
#            element = element,
#            type    = functional,
#            Z       = Zeff,
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

