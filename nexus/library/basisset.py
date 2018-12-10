##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


import os
from numpy import array,sqrt,arange,linspace,zeros,exp
from generic import obj
from periodic_table import is_element
from developer import DevBase,error
from fileio import TextFile
from plotting import *
from debug import *



def show_plots():
    show()
#end def show_plots

# container class for available basis set files
class BasisSets(DevBase):
    def __init__(self,*basissets):
        if len(basissets)==1 and isinstance(basissets[0],list):
            basissets = basissets[0]
        #end if
        bsfiles = []
        bss     = []
        errors = False
        for bs in basissets:
            if isinstance(bs,BasisFile):
                bss.append(bs)
            elif isinstance(bs,str):
                bsfiles.append(bs)
            else:
                self.error('expected BasisFile type or filepath, got '+str(type(bs)),exit=False)
                errors = True
            #end if
        #end for
        if errors:
            self.error('cannot create Basissets object')
        #end if

        if len(bss)>0:
            self.addbs(bss)
        #end if
        if len(bsfiles)>0:
            self.readbs(bsfiles)
        #end if
    #end def __init__


    def addbs(self,*basissets):
        if len(basissets)==1 and isinstance(basissets[0],list):
            basissets = basissets[0]
        #end if
        for bs in basissets:
            self[bs.filename] = bs
        #end for
    #end def addbs

        
    def readbs(self,*bsfiles):
        if len(bsfiles)==1 and isinstance(bsfiles[0],list):
            bsfiles = bsfiles[0]
        #end if
        bss = []
        print
        print '  Basissets'
        for filepath in bsfiles:
            print '    reading basis: ',filepath
            ext = filepath.split('.')[-1].lower()
            if ext=='gms_bas':
                bs = gamessBasisFile(filepath)
            else:
                bs = BasisFile(filepath)
            #end if
            bss.append(bs)
        #end for
        print
        self.addbs(bss)
    #end def readbs


    def bases_by_atom(self,*bsfiles):
        bss = obj()
        for bsfile in bsfiles:
            if bsfile in self:
                bs = self[bsfile]
                bss[bs.element_label] = bs
            else:
                self.error('basis file not found\nmissing file: {0}'.format(bsfile))
            #end if
        #end for
        return bss
    #end def bases_by_atom
#end class BasisSets


class BasisFile(DevBase):
    def __init__(self,filepath=None):
        self.element       = None
        self.element_label = None
        self.filename      = None
        self.location      = None
        if filepath!=None:
            self.filename = os.path.basename(filepath)
            self.location = os.path.abspath(filepath)
            elem_label = self.filename.split('.')[0]
            is_elem,symbol = is_element(elem_label,symbol=True)
            if not is_elem:
                self.error('cannot determine element for basis file: {0}\nbasis file names must be prefixed by an atomic symbol or label\n(e.g. Si, Si1, etc)'.format(filepath))
            #end if
            self.element = symbol
            self.element_label = elem_label
        #end if
    #end def __init__

    def cleaned_text(self):
        self.not_implemented()
    #end def cleaned_text
#end class BasisFile




class gaussBasisFile(BasisFile):
    angular_terms = 'spdfghiklmn'

    def __init__(self,filepath=None):
        BasisFile.__init__(self,filepath)
        self.text = None
        if filepath!=None:
            self.read(filepath)
        #end if
    #end def __init__

    def cleaned_text(self):
        if self.text is None:
            self.error('text requested prior to read\nfile: {0}'.format(self.location))
        #end if
        return self.text
    #end def cleaned_text

    def read(self,filepath=None):
        if filepath is None:
            filepath = self.location
        #end if
        if not os.path.exists(filepath):
            self.error('file does not exist: {0}'.format(filepath))
        #end if
        file = TextFile(filepath)
        self.read_file(file)
    #end def read

    def read_file(self,file):
        self.not_implemented()
    #end def read_file
#end class gaussBasisFile



class gamessBasisFile(gaussBasisFile):
    def read_file(self,file):
        dstart = file.find('$DATA')
        estart = file.find('$END')
        if dstart!=-1:
            file.seek(dstart)
            file.readline()
            cur = file.tell()
            tokens = file.readtokens()
            if len(tokens)==2 and tokens[0].lower() in self.angular_terms and tokens[1].isdigit():
                file.seek(cur) # skip atom header line
            #end if
            if estart==-1:
                estart = file.size()
            #end if
            self.text = file[file.tell():estart].strip() 
        else:
            lines = file[:].splitlines()
            blines = []
            for line in lines:
                ls = line.strip()
                if len(ls)>0 and ls[0]!='!' and ls[0]!='#':
                    blines.append(line)
                #end if
            #end for
            if len(blines)>0:
                tokens = blines[0].split()
                if not (len(tokens)==2 and tokens[0].lower() in self.angular_terms and tokens[1].isdigit()):
                    blines = blines[1:] # skip atom header line
                #end if
            #end if
            text = ''
            for line in blines:
                text += line + '\n'
            #end for
            self.text = text.strip()
        #end if
    #end def read_file
#end class gamessBasisFile





 

def process_gaussian_text(text,format,pp=True,basis=True,preserve_spacing=False):
    if format=='gamess' or format=='gaussian' or format=='atomscf':
        rawlines = text.splitlines()
        sections = []
        last_empty = True
        for rline in rawlines:
            line = rline.strip()
            if (not line.startswith('!')) and (not line.startswith('#')) and len(line)>0:
                if last_empty:
                    lines = []
                    sections.append(lines)
                #end if
                if preserve_spacing:
                    lines.append(rline)
                else:
                    lines.append(line)
                #end if
                last_empty = False
            else:
                last_empty = True
            #end if
        #end for
        del lines
        if len(sections)==2:
            basis_lines = sections[0]
            pp_lines    = sections[1]
        elif pp:
            basis_lines = None
            pp_lines    = sections[0]
        elif basis:
            basis_lines = sections[0]
            pp_lines    = None
        #end if
    elif format=='crystal':
        rawlines = text.splitlines()
        pp_lines = []
        basis_lines = []
        foundpp = False
        for line in rawlines:
            if not foundpp:
                foundpp = len(line.split())==5
            #end if
            if not foundpp:
                pp_lines.append(line)
            else:
                basis_lines.append(line)
            #end if
        #end for
        if len(pp_lines)==0:
            pp_lines = None
        #end if
        if len(basis_lines)==0:
            basis_lines = None
        #end if
    else:
        error('{0} format is unknown'.format(format),'process_gaussian_text')
    #end if
    if pp and basis:
        return pp_lines,basis_lines
    elif pp:
        return pp_lines
    elif basis:
        return basis_lines
    else:
        error('must request pp or basis')
    #end if
#end def process_gaussian_text




class GaussianBasisSet(DevBase):
    lset_full = tuple('spdfghijk')
    lstyles = obj(s='g-',p='r-',d='b-',f='m-',g='c-',h='k-',i='g-.',j='r-.',k='b-.')
    formats = 'gaussian gamess'.split()

    crystal_lmap = {0:'s',1:'sp',2:'p',3:'d',4:'f'}
    crystal_lmap_reverse = dict(s=0,sp=1,p=2,d=3,f=4)

    @staticmethod
    def process_float(s):
        return float(s.replace('D','e').replace('d','e'))
    #end def process_float


    def __init__(self,filepath=None,format=None):
        self.name  = None
        self.basis = obj()
        if filepath!=None:
            self.read(filepath,format)
        #end if
    #end def __init__


    def read(self,filepath,format=None):
        if format is None:
            self.error('format keyword must be specified to read file {0}\nvalid options are: {1}'.format(filepath,self.formats))
        elif not format in self.formats:
            self.error('incorrect format requested: {0}\nvalid options are: {1}'.format(format,self.formats))
        #end if
        if not os.path.exists(filepath):
            self.error('cannot read {0}, file does not exist'.format(filepath))
        #end if
        #self.name = split_delims(os.path.split(filepath)[1])[0]
        self.name = os.path.split(filepath)[1].split('.')[0]
        text = open(filepath,'r').read()
        self.read_text(text,format)
    #end def read

        
    def write(self,filepath=None,format=None):
        if format is None:
            self.error('format keyword must be specified to write file {0}\nvalid options are: {1}'.format(filepath,self.formats))
        elif not format in self.formats:
            self.error('incorrect format requested: {0}\nvalid options are: {1}'.format(format,self.formats))
        #end if
        text = self.write_text(format)
        if filepath!=None:
            open(filepath,'w').write(text)
        #end if
        return text
    #end def write


    def read_text(self,text,format=None):
        basis_lines = process_gaussian_text(text,format,pp=False)
        self.read_lines(basis_lines,format)
    #end def read_text


    def read_lines(self,basis_lines,format=None):
        basis = self.basis
        basis.clear()
        if format=='gamess':
            i=1
            while i<len(basis_lines):
                tokens = basis_lines[i].split(); i+=1
                ltext = tokens[0].lower()
                ngauss = int(tokens[1])
                scale  = array(tokens[2:],dtype=float)
                bterms = obj()
                for j in xrange(ngauss):
                    index,expon,coeff = basis_lines[i].split(); i+=1
                    expon = GaussianBasisSet.process_float(expon)
                    coeff = GaussianBasisSet.process_float(coeff)
                    bterms.append(obj(expon=expon,coeff=coeff))
                #end for
                basis.append(obj(l=ltext,scale=scale,terms=bterms))
            #end while
        #end if
        elif format=='gaussian':
            i=1
            while i<len(basis_lines):
                tokens = basis_lines[i].split(); i+=1
                ltext = tokens[0].lower()
                ngauss = int(tokens[1])
                scale  = array(tokens[2:],dtype=float)
                bterms = obj()
                for j in xrange(ngauss):
                    expon,coeff = basis_lines[i].split(); i+=1
                    expon = GaussianBasisSet.process_float(expon)
                    coeff = GaussianBasisSet.process_float(coeff)
                    bterms.append(obj(expon=expon,coeff=coeff))
                #end for
                basis.append(obj(l=ltext,scale=scale,terms=bterms))
            #end while
        elif format=='crystal':
            i=0
            while i<len(basis_lines):
                tokens = basis_lines[i].split(); i+=1
                if len(tokens)!=5:
                    self.error('could not parse crystal basisset, input may be misformatted')
                #end if
                basis_type    =   int(tokens[0])
                l_type        =   int(tokens[1])
                ngauss        =   int(tokens[2])
                formal_charge = float(tokens[3])
                scale         = array([tokens[4]],dtype=float)
                ltext = GaussianBasisSet.crystal_lmap[l_type]
                if ltext!='sp':
                    bterms = obj()
                    for j in xrange(ngauss):
                        expon,coeff = basis_lines[i].split(); i+=1
                        expon = GaussianBasisSet.process_float(expon)
                        coeff = GaussianBasisSet.process_float(coeff)
                        bterms.append(obj(expon=expon,coeff=coeff))
                    #end for
                    basis.append(obj(l=ltext,scale=scale,terms=bterms))
                else: # sp has shared exponent for s and p, split them now
                    sterms = obj()
                    pterms = obj()
                    for j in xrange(ngauss):
                        expon,scoeff,pcoeff = basis_lines[i].split(); i+=1
                        expon = GaussianBasisSet.process_float(expon)
                        scoeff = GaussianBasisSet.process_float(scoeff)
                        pcoeff = GaussianBasisSet.process_float(pcoeff)
                        sterms.append(obj(expon=expon,coeff=scoeff))
                        pterms.append(obj(expon=expon,coeff=pcoeff))
                    #end for
                    basis.append(obj(l='s',scale=scale,terms=sterms))
                    basis.append(obj(l='p',scale=scale,terms=pterms))
                #end if
            #end while
        else:
            self.error('ability to read file format {0} has not been implemented'.format(format))
        #end if
        # sort the basis in s,p,d,f,... order
        self.lsort()
    #end def read_lines


    
    def write_text(self,format=None,occ=None):
        text = ''
        format = format.lower()
        if format=='gamess':
            #text += '{0} {1} 0. 0. 0.\n'.format(self.element,self.Zcore+self.Zval)
            for ib in xrange(len(self.basis)):
                b = self.basis[ib]
                line = '{0} {1}'.format(b.l,len(b.terms))
                for s in b.scale:
                    line += ' {0}'.format(s)
                #end for
                text += line + '\n'
                for it in xrange(len(b.terms)):
                    t = b.terms[it]
                    text += '{0} {1:12.8f} {2: 12.8f}\n'.format(it+1,t.expon,t.coeff)
                #end for
            #end for
        elif format=='gaussian':
            #text += '{0} 0\n'.format(self.element)
            for ib in xrange(len(self.basis)):
                b = self.basis[ib]
                line = '{0} {1}'.format(b.l,len(b.terms))
                for s in b.scale:
                    line += ' {0}'.format(s)
                #end for
                text += line + '\n'
                for it in xrange(len(b.terms)):
                    t = b.terms[it]
                    text += '{0:12.8f}{1: 12.8f}\n'.format(t.expon,t.coeff)
                #end for
            #end for
        elif format=='crystal':
            if occ is not None:
                lcounts = dict(s=0,p=0,d=0,f=0)
            #end if
            for ib in xrange(len(self.basis)):
                b = self.basis[ib]
                if b.l not in self.crystal_lmap_reverse:
                    self.error('{0} channels cannot be handled by crystal'.format(b.l))
                #end if
                Zf = 0
                if occ is not None and b.l in occ and lcounts[b.l]<len(occ[b.l]):
                    Zf = occ[b.l][lcounts[b.l]]
                    lcounts[b.l]+=1
                #end if
                lnum = self.crystal_lmap_reverse[b.l]
                line = '0 {0} {1} {2} {3}'.format(lnum,len(b.terms),Zf,b.scale[0])
                text += line + '\n'
                for it in xrange(len(b.terms)):
                    t = b.terms[it]
                    text += '{0:12.8f}{1: 12.8f}\n'.format(t.expon,t.coeff)
                #end for
            #end for
        else:
            self.error('ability to write file format {0} has not been implemented'.format(format))
        #end if
        return text
    #end def write_text


    def size(self):
        return len(self.basis)
    #end def size


    def lset(self):
        lset = set()
        for bf in self.basis:
            lset.add(bf.l)
        #end for
        return lset
    #end def lset


    def lcount(self):
        return len(self.lset())
    #end def lcount


    def lbasis(self):
        lbasis = obj()
        for n in range(len(self.basis)):
            bf = self.basis[n]
            l  = bf.l
            if l not in lbasis:
                lbasis[l] = obj()
            #end if
            lbasis[l].append(bf)
        #end for
        return lbasis
    #end def lbasis

    
    def lsort(self):
        lbasis = self.lbasis()
        self.basis.clear()
        for l in self.lset_full:
            if l in lbasis:
                lbas = lbasis[l]
                for n in range(len(lbas)):
                    bf = lbas[n]
                    self.basis.append(bf)
                #end for
            #end if
        #end for
    #end def lsort


    def uncontracted(self):
        all_uncon = True
        for bf in self.basis:
            all_uncon &= len(bf.terms)==1
        #end for
        return all_uncon
    #end def uncontracted


    def contracted(self):
        return not self.uncontracted()
    #end def contracted


    def uncontract(self,tol=1e-3):
        if self.uncontracted():
            return
        #end if
        lbasis = self.lbasis()
        self.basis.clear()
        for l in self.lset_full:
            if l in lbasis:
                exponents = []
                lbas = lbasis[l]
                for n in xrange(len(lbas)):
                    uterms = lbas[n].terms
                    for i in xrange(len(uterms)):
                        expon = uterms[i].expon
                        if len(exponents)==0:
                            exponents = array([expon],dtype=float)
                        elif abs(exponents-expon).min()>tol:
                            exponents = array(list(exponents)+[expon],dtype=float)
                        #end if
                    #end for
                #end for
                for expon in exponents:
                    cterms = obj()
                    cterms.append(obj(expon=expon,coeff=1.0))
                    bf = obj(l=l,scale=array([1.0]),terms=cterms)
                    self.basis.append(bf)
                #end for
            #end if
        #end for
    #end def uncontract


    def contracted_basis_size(self):
        bcount = obj()
        for bf in self.basis:
            l = bf.l
            if l not in bcount:
                bcount[l]=0
            #end if
            bcount[l] += 1
        #end for
        bs = ''
        for l in self.lset_full:
            if l in bcount:
                bs += str(bcount[l])+l
            #end if
        #end for
        return bs
    #end def contracted_basis_size


    def uncontracted_basis_size(self):
        if self.uncontracted():
            return self.contracted_basis_size()
        #end if
        uc = self.copy()
        uc.uncontract()
        return uc.contracted_basis_size()
    #end def uncontracted_basis_size


    def basis_size(self):
        us = self.uncontracted_basis_size()
        cs = self.contracted_basis_size()
        return '({0})/[{1}]'.format(us,cs)
    #end def basis_size


    def prim_expons(self):
        if self.contracted():
            self.error('cannot find primitive gaussian expons because basis is contracted')
        #end if
        lbasis = self.lbasis()
        gexpon = obj()
        for l,lbas in lbasis.iteritems():
            e = []
            for n in range(len(lbas)):
                e.append(lbas[n].terms[0].expon)
            #end for
            gexpon[l] = array(e,dtype=float)
        #end for
        return gexpon
    #end def prim_expons


    def prim_widths(self):
        if self.contracted():
            self.error('cannot find primitive gaussian widths because basis is contracted')
        #end if
        lbasis = self.lbasis()
        gwidth = obj()
        for l,lbas in lbasis.iteritems():
            w = []
            for n in range(len(lbas)):
                w.append(1./sqrt(2.*lbas[n].terms[0].expon))
            #end for
            gwidth[l] = array(w,dtype=float)
        #end for
        return gwidth
    #end def prim_widths

    
    def remove_prims(self,comp=None,keep=None,**lselectors):
        lbasis = self.lbasis()
        if comp!=None:
            gwidths = self.prim_widths()
        #end if
        for l,lsel in lselectors.iteritems():
            if l not in lbasis:
                self.error('cannot remove basis functions from channel {0}, channel not present'.format(l))
            #end if
            lbas = lbasis[l]
            if isinstance(lsel,float):
                rcut = lsel
                less = False
                if comp=='<':
                    less = True
                elif comp=='>':
                    less = False
                elif comp is None:
                    self.error('comp argument must be provided (< or >)')
                else:
                    self.error('comp must be < or >, you provided: {0}'.format(comp))
                #end if
                gw = gwidths[l]
                iw = arange(len(gw))
                nkeep = 0
                if keep!=None and l in keep:
                    nkeep = keep[l]
                #end if
                if less:
                    rem = iw[gw<rcut]
                    for i in xrange(len(rem)-nkeep):
                        del lbas[rem[i]]
                    #end for
                else:
                    rem = iw[gw>rcut]
                    for i in xrange(nkeep,len(rem)):
                        del lbas[rem[i]]
                    #end for
                #end if
            elif isinstance(lsel,int):                
                if comp=='<':
                    if lsel>len(lbas):
                        self.error('cannot remove {0} basis functions from channel {1} as it only has {2}'.format(lsel,l,len(lbas)))
                    #end if
                    for i in xrange(lsel):
                        del lbas[i]
                    #end for
                elif comp=='>':
                    if lsel>len(lbas):
                        self.error('cannot remove {0} basis functions from channel {1} as it only has {2}'.format(lsel,l,len(lbas)))
                    #end if
                    for i in xrange(len(lbas)-lsel,len(lbas)):
                        del lbas[i]
                    #end for
                else:
                    if lsel>=len(lbas):
                        self.error('cannot remove basis function {0} from channel {1} as it only has {2}'.format(lsel,l,len(lbas)))
                    #end if
                    del lbas[lsel]
                #end if
            else:
                for ind in lsel:
                    del lbas[ind]
                #end for
            #end if
        #end for
        self.basis.clear()
        for l in self.lset_full:
            if l in lbasis:
                lbas = lbasis[l]
                for k in sorted(lbas.keys()):
                    self.basis.append(lbas[k])
                #end for
            #end if
        #end for
    #end def remove_prims


    def remove_small_prims(self,**keep):
        lsel = obj()
        for l,lbas in self.lbasis().iteritems():
            if l in keep:
                lsel[l] = len(lbas)-keep[l]
            #end if
        #end for
        self.remove_prims(comp='<',**lsel)
    #end def remove_small_prims


    def remove_large_prims(self,**keep):
        lsel = obj()
        for l,lbas in self.lbasis().iteritems():
            if l in keep:
                lsel[l] = len(lbas)-keep[l]
            #end if
        #end for
        self.remove_prims(comp='>',**lsel)
    #end def remove_large_prims


    def remove_small_prims_rel(self,other,**keep):
        gwidths = other.prim_widths()
        lsel = obj()
        for l,gw in gwidths.iteritems():
            lsel[l] = gw.min()
        #end for
        self.remove_prims(comp='<',keep=keep,**lsel)
    #end def remove_small_prims_rel


    def remove_large_prims_rel(self,other,**keep):
        gwidths = other.prim_widths()
        lsel = obj()
        for l,gw in gwidths.iteritems():
            lsel[l] = gw.max()
        #end for
        self.remove_prims(comp='>',keep=keep,**lsel)
    #end def remove_large_prims_rel


    def remove_channels(self,llist):
        lbasis = self.lbasis()
        for l in llist:
            if l in lbasis:
                del lbasis[l]
            #end if
        #end for
        self.basis.clear()
        for l in self.lset_full:
            if l in lbasis:
                lbas = lbasis[l]
                for n in range(len(lbas)):
                    bf = lbas[n]
                    self.basis.append(bf)
                #end for
            #end if
        #end for
    #end def remove_channels
                

    def incorporate(self,other,tol=1e-3,unique=False):
        uncontracted = self.uncontracted() and other.uncontracted()
        lbasis       = self.lbasis()
        lbasis_other = other.lbasis()
        if uncontracted:
            gwidths       = self.prim_widths()
            gwidths_other = other.prim_widths()
        #end if
        self.basis.clear()
        if not unique: # simple, direct merge of basis sets
            for l in self.lset_full:
                if l in lbasis:
                    lbas = lbasis[l]
                    for n in range(len(lbas)):
                        bf = lbas[n]
                        self.basis.append(bf)
                    #end for
                #end if
                if l in lbasis_other:
                    lbas = lbasis_other[l]
                    for n in range(len(lbas)):
                        bf = lbas[n]
                        self.basis.append(bf)
                    #end for
                #end if
            #end for
        else: # merge uncontracted basis sets preserving order
            for l in self.lset_full:
                primitives = []
                widths     = []
                orig_widths = array([])
                if l in lbasis:
                    primitives.extend(lbasis[l].list())
                    widths.extend(gwidths[l])
                    orig_widths = gwidths[l]
                #end if
                if l in lbasis_other:
                    prims = lbasis_other[l].list()
                    owidths = gwidths_other[l]
                    for n in range(len(prims)):
                        w = owidths[n]
                        if len(orig_widths)==0 or abs(orig_widths-w).min()>tol:
                            primitives.append(prims[n])
                            widths.append(w)
                        #end if
                    #end if
                #end if
                primitives = array(primitives,dtype=object)[array(widths).argsort()]
                for bf in primitives:
                    self.basis.append(bf)
                #end for
            #end for
        #end if
    #end def incorporate


    def plot(self,r=None,rmin=0.01,rmax=8.0,show=True,fig=True,sep=False,prim=False,style=None,fmt=None,nsub=None):
        if r is None:
            r = linspace(rmin,rmax,1000)
        #end if
        if not prim:
            ptitle = '{0} {1} basis'.format(self.name,self.basis_size())
        else:
            ptitle = '{0} {1} primitives'.format(self.name,self.basis_size())
        #end if
        if fig:
            figure()
        #end if
        r2 = r**2
        lcount = self.lcount()
        if nsub!=None:
            lcount = max(lcount,nsub)
        #end if
        lbasis = self.lbasis()
        lc = 0
        for l in self.lset_full:
            if l in lbasis:
                lc+=1
                if sep:
                    subplot(lcount,1,lc)
                    ylabel(l)
                    if lc==1:
                        title(ptitle)
                    #end if
                #end if
                lstyle=self.lstyles[l]
                if style!=None:
                    lstyle = lstyle[0]+style
                #end if
                if fmt!=None:
                    lstyle=fmt
                #end if
                lbas = lbasis[l]
                for n in range(len(lbas)):
                    bf = lbas[n]
                    br = zeros(r.shape)
                    s  = bf.scale[0]
                    for pf in bf.terms:
                        c =  pf.coeff
                        a = -pf.expon*s**2
                        pr = exp(a*r2)
                        if not prim:
                            br += c*pr
                        else:
                            plot(r,pr,lstyle,label=l)
                        #end if
                    #end for
                    if not prim:
                        plot(r,br,lstyle,label=l)
                    #end if
                #end for
            #end if
        #end for
        if fig:
            if not sep:
                if self.name!=None:
                    title(ptitle)
                #end if
                ylabel('basis functions')
                legend()
            #end if
            xlabel('r')
            if show:
                show_plots()
            #end if
        #end if
    #end def plot


    def plot_primitives(self):
        None
    #end def plot_primitives


    def plot_prim_widths(self,show=True,fig=True,sep=False,style='o',fmt=None,nsub=None,semilog=True,label=True):
        if self.contracted():
            self.error('cannot plot primitive gaussian widths because basis is contracted')
        #end if
        ptitle = '{0} {1} primitive widths'.format(self.name,self.basis_size())
        if fig:
            figure()
        #end if
        pwidths = self.prim_widths()
        lcount = self.lcount()
        if nsub!=None:
            lcount = max(lcount,nsub)
        #end if
        lbasis = self.lbasis()
        lc = 0
        for l in self.lset_full:
            if l in lbasis:
                lwidths = pwidths[l]
                lc+=1
                if sep:
                    subplot(lcount,1,lc)
                    ylabel(l)
                    if lc==1:
                        title(ptitle)
                    #end if
                #end if
                lstyle=self.lstyles[l]
                if style!=None:
                    lstyle = lstyle[0]+style
                #end if
                if fmt!=None:
                    lstyle=fmt
                #end if
                if semilog:
                    plotter = semilogy
                else:
                    plotter = plot
                #end if
                plotter(range(len(lwidths)),lwidths,lstyle,label=l)
            #end if
        #end for
        if label:
            if not sep:
                if self.name!=None:
                    title(ptitle)
                #end if
                ylabel('primitive widths')
                legend()
            #end if
            xlabel('primitive index')
        #end if
        if show:
            show_plots()
        #end if
    #end def plot_prim_widths
#end class GaussianBasisSet
