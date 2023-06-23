##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  fileio.py                                                         #
#    Support for I/O with various file formats.  Currently this only #
#    contains a generic file I/O class for XSF files.  In the future #
#    generic XML and HDF5 support should go here.  Input only        #
#    interfaces to these formats can be found in hdfreader.py and    #
#    xmlreader.py.                                                   #
#                                                                    #
#  Content summary:                                                  #
#    XsfFile                                                         #
#      Represents generic XSF, AXSF, and BXSF files.                 #
#      Can read/write arbitrary files of these formats.              #
#      Useful for atomic structure and electronic density I/O.       #       
#                                                                    #
#====================================================================#


import os
import mmap
import numpy as np
from numpy import array,zeros,ndarray,around,arange,dot,savetxt,empty,reshape
from numpy.linalg import det,norm
from generic import obj
from developer import DevBase,error,to_str
from periodic_table import pt as ptable,is_element
from unit_converter import convert
from debug import *


class TextFile(DevBase):
    # interface to mmap files
    # see Python 2 documentation for mmap

    def __init__(self,filepath=None):
        self.mm = None
        self.f  = None
        if filepath!=None:
            self.open(filepath)
        #end if
    #end def __init__

    def open(self,filepath):
        if not os.path.exists(filepath):
            self.error('cannot open non-existent file: {0}'.format(filepath))
        #end if
        f = open(filepath,'r')
        fno = f.fileno()
        #fno = os.open(filepath,os.O_RDONLY)
        self.f = f
        self.mm = mmap.mmap(fno,0,prot=mmap.PROT_READ)
    #end def open

    def __iter__(self):
        for line in self.f:
            yield line
        #end for
    #end def __iter__

    def __getitem__(self,slc):
        return self.mm[slc]
    #end def __getitem__

    def lines(self):
        return self.read().splitlines()
    #end def lines

    def tokens(self):
        return self.read().split()
    #end def tokens

    def readtokens(self,s=None):
        return self.readline(s).split()
    #end def readtokens

    def readtokensf(self,s=None,*formats):
        if s!=None:
            self.seek(s)
        #end if
        self.mm.readline()
        line = to_str(self.mm.readline())
        stokens = line.split()
        all_same = False
        if len(formats)==1 and len(stokens)>1:
            format = formats[0]
            all_same = True
        elif len(formats)>len(stokens):
            self.error('formatted line read failed\nnumber of tokens and provided number of formats do not match\nline: {0}\nnumber of tokens: {1}\nnumber of formats provided: {2}'.format(line,len(stokens),len(formats)))
        #end if
        tokens = []
        if all_same:
            for stoken in stokens:
                tokens.append(format(stoken))
            #end for
        else:
            for format,stoken in zip(formats,stokens):
                tokens.append(format(stoken))
            #end for
        #end if
        if len(tokens)==1:
            return tokens[0]
        else:
            return tokens
        #end if
    #end def readtokensf

    # extended mmap interface below
    def close(self):
        r = self.mm.close()
        self.f.close()
        return r
    #end def close

    def seek(self,pos,whence=0,start=None,end=None):
        if isinstance(pos,str):
            pos = pos.encode('ASCII')
            if whence!=2 and start is None:
                if whence==0:
                    start = 0
                elif whence==1:
                    start = self.mm.tell()
                else:
                    self.error('relative positioning must be either 0 (begin), 1 (current), or 2 (end)\nyou provided: {0}'.format(whence))
                #end if
            #end if
            if whence!=2:
                if end!=None:
                    pos = self.mm.find(pos,start,end)
                else:
                    pos = self.mm.find(pos,start)
                #end if
            else:
                if end!=None:
                    pos = self.mm.rfind(pos,start,end)
                else:
                    pos = self.mm.rfind(pos,start)
                #end if
            #end if
            if pos!=-1:
                return self.mm.seek(pos,0)
            else:
                return -1
            #end if
        else:
            return self.mm.seek(pos,whence)
        #end if
    #end def seek

    def readline(self,s=None):
        if s!=None:
            self.seek(s)
        #end if
        return to_str(self.mm.readline())
    #end def readline

    def read(self,num=None):
        if num is None:
            return to_str(self.mm[:])
        else:
            return to_str(self.mm.read(num))
        #end if
    #end def read


    # unchanged mmap interface below
    def find(self,*a,**kw):
        args = []
        for v in a:
            if isinstance(v,str):
                args.append(v.encode('ASCII'))
            else:
                args.append(a)
            #end if
        #end for
        return self.mm.find(*args,**kw)
    #end def find

    def flush(self,*a,**kw):
        return self.mm(*a,**kw)
    #end def flush

    def move(self,dest,src,count):
        return self.mm.move(dest,src,count)
    #end def move

    def read_byte(self):
        return self.mm.read_byte()
    #end def read_byte

    def resize(self,newsize):
        return self.mm.resize(newsize)
    #end def resize

    def rfind(self,*a,**kw):
        args = []
        for v in a:
            if isinstance(v,str):
                args.append(v.encode('ASCII'))
            else:
                args.append(a)
            #end if
        #end for
        return self.mm.rfind(*args,**kw)
    #end def rfind

    def size(self):
        return self.mm.size()
    #end def size

    def tell(self):
        return self.mm.tell()
    #end def tell

    def write(self,string):
        return self.mm.write(string)
    #end def write

    def write_byte(self,byte):
        return self.mm.write_byte(byte)
    #end def write_byte
#end class TextFile



class StandardFile(DevBase):

    sftype = ''

    def __init__(self,filepath=None):
        if filepath is None:
            None
        elif isinstance(filepath,str):
            self.read(filepath)
        else:
            self.error('unsupported input: {0}'.format(filepath))
        #end if
    #end def __init__


    def read(self,filepath):
        if not os.path.exists(filepath):
            self.error('read failed\nfile does not exist: {0}'.format(filepath))
        #end if
        self.read_text(open(filepath,'r').read())
        self.check_valid('read failed')
    #end def read


    def write(self,filepath=None):
        self.check_valid('write failed')
        text = self.write_text()
        if filepath!=None:
            open(filepath,'w').write(text)
        #end if
        return text
    #end def write


    def is_valid(self):
        return len(self.validity_checks())==0
    #end def is_valid


    def check_valid(self,header=None):
        messages = self.validity_checks()
        if len(messages)>0:
            msg = ''
            if header is not None:
                msg += header+'\n'
            #end if
            msg += 'not a valid {0} file, see below for details\n'.format(self.sftype)
            for m in messages:
                msg+=m+'\n'
            #end for
            self.error(msg)
        #end if
    #end def check_valid


    def validity_checks(self):
        messages = []
        return messages
    #end def validity_checks


    def read_text(self,text):
        self.not_implemented()
    #end def read_text


    def write_text(self):
        self.not_implemented()
    #end def write_text

#end class StandardFile



class XsfFile(StandardFile):

    sftype = 'xsf'

    filetypes     = set(['xsf','axsf','bxsf'])
    periodicities = set(['molecule','polymer','slab','crystal']) 
    dimensions    = obj(molecule=0,polymer=1,slab=2,crystal=3)

    # ATOMS  are in units of Angstrom, only provided for 'molecule'
    # forces are in units of Hatree/Angstrom
    # each section should be followed by a blank line

    def __init__(self,filepath=None):
        self.filetype    = None
        self.periodicity = None
        StandardFile.__init__(self,filepath)
    #end def __init__


    def add_to_image(self,image,name,value):
        if image is None:
            self[name] = value
        else:
            if 'images' not in self:
                self.images = obj()
            #end if
            if not image in self.images:
                self.images[image] = obj()
            #end if
            self.images[image][name] = value
        #end if
    #end def add_to_image


    # test needed for axsf and bxsf
    def read_text(self,text):
        lines = text.splitlines()
        i=0
        self.filetype = 'xsf'
        while(i<len(lines)):
            line = lines[i].strip().lower()
            if len(line)>0 and line[0]!='#':
                tokens = line.split()
                keyword = tokens[0]
                image = None
                if len(tokens)==2:
                    image = int(tokens[1])
                #end if
                if keyword in self.periodicities:
                    self.periodicity = keyword
                elif keyword=='animsteps':
                    self.animsteps = int(tokens[1])
                    self.filetype = 'axsf'
                elif keyword=='primvec':
                    primvec = array((lines[i+1]+' '+
                                     lines[i+2]+' '+
                                     lines[i+3]).split(),dtype=float)
                    primvec.shape = 3,3
                    self.add_to_image(image,'primvec',primvec)
                    i+=3
                elif keyword=='convvec':
                    convvec = array((lines[i+1]+' '+
                                     lines[i+2]+' '+
                                     lines[i+3]).split(),dtype=float)
                    convvec.shape = 3,3
                    self.add_to_image(image,'convvec',convvec)
                    i+=3
                elif keyword=='atoms':
                    if self.periodicity is None:
                        self.periodicity='molecule'
                    #end if
                    i+=1
                    tokens = lines[i].strip().split()
                    elem  = []
                    pos   = []
                    force = []
                    natoms = 0
                    while len(tokens)==4 or len(tokens)==7:
                        natoms+=1
                        elem.append(tokens[0])
                        pos.extend(tokens[1:4])
                        if len(tokens)==7:
                            force.extend(tokens[4:7])
                        #end if
                        i+=1
                        tokens = lines[i].strip().split()
                    #end while
                    elem = array(elem,dtype=int)
                    pos  = array(pos,dtype=float)
                    pos.shape = natoms,3
                    self.add_to_image(image,'elem',elem)
                    self.add_to_image(image,'pos',pos)
                    if len(force)>0:
                        force = array(force,dtype=float)
                        force.shape = natoms,3
                        self.add_to_image(image,'force',force)
                    #end if
                    i-=1
                elif keyword=='primcoord':
                    natoms = int(lines[i+1].split()[0])
                    elem  = []
                    pos   = []
                    force = []
                    for iat in range(natoms):
                        tokens = lines[i+2+iat].split()
                        elem.append(tokens[0])
                        pos.extend(tokens[1:4])
                        if len(tokens)==7:
                            force.extend(tokens[4:7])
                        #end if
                    #end for
                    try:
                        elem = array(elem,dtype=int)
                    except:
                        elem = array(elem,dtype=str)
                    #end try
                    pos  = array(pos,dtype=float)
                    pos.shape = natoms,3
                    self.add_to_image(image,'elem',elem)
                    self.add_to_image(image,'pos',pos)
                    if len(force)>0:
                        force = array(force,dtype=float)
                        force.shape = natoms,3
                        self.add_to_image(image,'force',force)
                    #end if
                    i+=natoms+1
                elif keyword.startswith('begin_block_datagrid'):
                    if keyword.endswith('2d'):
                        d=2
                    elif keyword.endswith('3d'):
                        d=3
                    else:
                        self.error('dimension of datagrid could not be identified: '+line)
                    #end if
                    i+=1
                    block_identifier = lines[i].strip().lower()
                    if not 'data' in self:
                        self.data = obj()
                    #end if
                    if not d in self.data:
                        self.data[d] = obj()
                    #end if
                    if not block_identifier in self.data[d]:
                        self.data[d][block_identifier]=obj()
                    #end if
                    data = self.data[d][block_identifier]

                    line = ''
                    while not line.startswith('end_block_datagrid'):
                        line = lines[i].strip().lower()
                        if line.startswith('begin_datagrid') or line.startswith('datagrid_'):
                            grid_identifier = line.replace('begin_datagrid_{0}d_'.format(d),'')
                            grid   = array(lines[i+1].split(),dtype=int)[:d]
                            corner = array(lines[i+2].split(),dtype=float)
                            if d==2:
                                cell   = array((lines[i+3]+' '+
                                                lines[i+4]).split(),dtype=float)
                                i+=5
                            elif d==3:
                                cell   = array((lines[i+3]+' '+
                                                lines[i+4]+' '+
                                                lines[i+5]).split(),dtype=float)
                                i+=6
                            #end if
                            cell.shape = d,3
                            dtokens = []
                            line = lines[i].strip().lower()
                            while not line.startswith('end_datagrid'):
                                dtokens.extend(line.split())
                                i+=1
                                line = lines[i].strip().lower()
                            #end while
                            grid_data = array(dtokens,dtype=float)
                            grid_data=reshape(grid_data,grid,order='F')
                            data[grid_identifier] = obj(
                                grid   = grid,
                                corner = corner,
                                cell   = cell,
                                values = grid_data
                                )
                        #end if
                        i+=1
                    #end while
                elif keyword=='begin_info':
                    self.info = obj()
                    while line.lower()!='end_info':
                        line = lines[i].strip()
                        if len(line)>0 and line[0]!='#' and ':' in line:
                            k,v = line.split(':')
                            self.info[k.strip()] = v.strip()
                        #end if
                        i+=1
                    #end while
                elif keyword.startswith('begin_block_bandgrid'):
                    self.filetype = 'bxsf'
                    if keyword.endswith('2d'):
                        d=2
                    elif keyword.endswith('3d'):
                        d=3
                    else:
                        self.error('dimension of bandgrid could not be identified: '+line)
                    #end if
                    i+=1
                    block_identifier = lines[i].strip().lower()
                    if not 'band' in self:
                        self.band = obj()
                    #end if
                    if not d in self.band:
                        self.band[d] = obj()
                    #end if
                    if not block_identifier in self.band[d]:
                        self.band[d][block_identifier]=obj()
                    #end if
                    band = self.band[d][block_identifier]

                    line = ''
                    while not line.startswith('end_block_bandgrid'):
                        line = lines[i].strip().lower()
                        if line.startswith('begin_bandgrid'):
                            grid_identifier = line.replace('begin_bandgrid_{0}d_'.format(d),'')
                            nbands = int(lines[i+1].strip())
                            grid   = array(lines[i+2].split(),dtype=int)[:d]
                            corner = array(lines[i+3].split(),dtype=float)
                            if d==2:
                                cell   = array((lines[i+4]+' '+
                                                lines[i+5]).split(),dtype=float)
                                i+=6
                            elif d==3:
                                cell   = array((lines[i+4]+' '+
                                                lines[i+5]+' '+
                                                lines[i+6]).split(),dtype=float)
                                i+=7
                            #end if
                            cell.shape = d,3
                            bands = obj()
                            line = lines[i].strip().lower()
                            while not line.startswith('end_bandgrid'):
                                if line.startswith('band'):
                                    band_index = int(line.split(':')[1].strip())
                                    bands[band_index] = []
                                else:
                                    bands[band_index].extend(line.split())
                                #end if
                                i+=1
                                line = lines[i].strip().lower()
                            #end while
                            for bi,bv in bands.items():
                                bands[bi] = array(bv,dtype=float)
                                bands[bi].shape = tuple(grid)
                            #end for
                            band[grid_identifier] = obj(
                                grid   = grid,
                                corner = corner,
                                cell   = cell,
                                bands  = bands
                                )
                        #end if
                        i+=1
                    #end while
                else:
                    self.error('invalid keyword encountered: {0}'.format(keyword))
                #end if
            #end if
            i+=1
        #end while
    #end def read_text


    # test needed for axsf and bxsf
    def write_text(self):
        c=''
        if self.filetype=='xsf':    # only write structure/datagrid if present
            if self.periodicity=='molecule' and 'elem' in self:
                c += self.write_coord()
            elif 'primvec' in self:
                c += ' {0}\n'.format(self.periodicity.upper())
                c += self.write_vec('primvec',self.primvec)
                if 'convvec' in self:
                    c += self.write_vec('convvec',self.convvec)
                #end if
                if 'elem' in self:
                    c+= self.write_coord()
                #end if
            #end if
            if 'data' in self:
                c += self.write_data()
            #end if
        elif self.filetype=='axsf': # only write image structures
            c += ' ANIMSTEPS {0}\n'.format(self.animsteps)
            if self.periodicity!='molecule':
                c += ' {0}\n'.format(self.periodicity.upper())
            #end if
            if 'primvec' in self:
                c += self.write_vec('primvec',self.primvec)
            #end if
            if 'convvec' in self:
                c += self.write_vec('convvec',self.convvec)
            #end if
            for i in range(1,len(self.images)+1):
                image = self.images[i]
                if 'primvec' in image:
                    c += self.write_vec('primvec',image.primvec,i)
                #end if
                if 'convvec' in image:
                    c += self.write_vec('convvec',image.convvec,i)
                #end if
                c += self.write_coord(image,i)
            #end for
        elif self.filetype=='bxsf': # only write bandgrid
            c += self.write_band()
        #end if
        return c
    #end def write_text


    def write_coord(self,image=None,index=''):
        if image is None:
            s = self
        else:
            s = image
        #end if
        c = ''
        if self.periodicity=='molecule':
            c += ' ATOMS {0}\n'.format(index)
        else:
            c += ' PRIMCOORD {0}\n'.format(index)
            c += '   {0} 1\n'.format(len(s.elem))
        if not 'force' in s:
            for i in range(len(s.elem)):
                r = s.pos[i]
                c += '   {0:>3} {1:12.8f} {2:12.8f} {3:12.8f}\n'.format(s.elem[i],r[0],r[1],r[2])
            #end for
        else:
            for i in range(len(s.elem)):
                r = s.pos[i]
                f = s.force[i]
                c += '   {0:>3} {1:12.8f} {2:12.8f} {3:12.8f}  {4:12.8f} {5:12.8f} {6:12.8f}\n'.format(s.elem[i],r[0],r[1],r[2],f[0],f[1],f[2])
            #end for
        #end if
        return c
    #end def write_coord


    def write_vec(self,name,vec,index=''):
        c = ' {0} {1}\n'.format(name.upper(),index)
        for v in vec:
            c += '   {0:12.8f} {1:12.8f} {2:12.8f}\n'.format(v[0],v[1],v[2])
        #end for
        return c
    #end def write_vec


    def write_data(self):
        c = ''
        ncols = 4
        data = self.data
        for d in sorted(data.keys()):
            bdg_xd = data[d]       # all block datagrids 2 or 3 D
            for bdgk in sorted(bdg_xd.keys()):
                c += ' BEGIN_BLOCK_DATAGRID_{0}D\n'.format(d)
                c += '   {0}\n'.format(bdgk)
                bdg = bdg_xd[bdgk] # single named block data grid
                for dgk in sorted(bdg.keys()):
                    c += '   BEGIN_DATAGRID_{0}D_{1}\n'.format(d,dgk)
                    dg = bdg[dgk]  # single named data grid
                    if d==2:
                        c += '     {0} {1}\n'.format(*dg.grid)
                    elif d==3:
                        c += '     {0} {1} {2}\n'.format(*dg.grid)
                    #end if
                    c += '   {0:12.8f} {1:12.8f} {2:12.8f}\n'.format(*dg.corner)
                    for v in dg.cell:
                        c += '   {0:12.8f} {1:12.8f} {2:12.8f}\n'.format(*v)
                    #end for
                    c = c[:-1]
                    n=0
                    for v in dg.values.ravel(order='F'):
                        if n%ncols==0:
                            c += '\n    '
                        #end if
                        c += ' {0:12.8f}'.format(v)
                        n+=1
                    #end for
                    c += '\n   END_DATAGRID_{0}D_{1}\n'.format(d,dgk)
                #end for
                c += ' END_BLOCK_DATAGRID_{0}D\n'.format(d)
            #end for
        #end for                    
        return c
    #end def write_data


    def write_band(self):
        c = ''
        ncols = 4
        band = self.band
        for d in sorted(band.keys()):
            bdg_xd = band[d]       # all block bandgrids 2 or 3 D
            for bdgk in sorted(bdg_xd.keys()):
                c += ' BEGIN_BLOCK_BANDGRID_{0}D\n'.format(d)
                c += '   {0}\n'.format(bdgk)
                bdg = bdg_xd[bdgk] # single named block band grid
                for dgk in sorted(bdg.keys()):
                    c += '   BEGIN_BANDGRID_{0}D_{1}\n'.format(d,dgk)
                    dg = bdg[dgk]  # single named band grid
                    if d==2:
                        c += '     {0} {1}\n'.format(*dg.grid)
                    elif d==3:
                        c += '     {0} {1} {2}\n'.format(*dg.grid)
                    #end if
                    c += '   {0:12.8f} {1:12.8f} {2:12.8f}\n'.format(*dg.corner)
                    for v in dg.cell:
                        c += '   {0:12.8f} {1:12.8f} {2:12.8f}\n'.format(*v)
                    #end for
                    for bi in sorted(dg.bands.keys()):
                        c += '   BAND:  {0}'.format(bi)
                        n=0
                        for v in dg.bands[bi].ravel():
                            if n%ncols==0:
                                c += '\n    '
                            #end if
                            c += ' {0:12.8f}'.format(v)
                            n+=1
                        #end for
                        c += '\n'
                    #end for
                    c += '   END_BANDGRID_{0}D_{1}\n'.format(d,dgk)
                #end for
                c += ' END_BLOCK_BANDGRID_{0}D\n'.format(d)
            #end for
        #end for
        return c
    #end def write_band


    def dimension(self):
        if self.periodicity in self.dimensions:
            return self.dimensions[self.periodicity]
        else:
            return None
        #end if
    #end def dimension


    def initialized(self):
        return self.filetype!=None
    #end def initialized


    def has_animation(self):
        return self.filetype=='axsf' and 'animsteps' in self
    #end def has_animation


    def has_bands(self):
        return self.filetype=='bxsf' and 'band' in self and 'info' in self
    #end def has_bands


    def has_structure(self):
        hs = self.filetype=='xsf'
        hs &= 'elem' in self and 'pos' in self
        d = self.dimension()
        if d!=0:
            hs &= 'primvec' in self
        #end if
        return hs
    #end def has_structure


    def has_data(self):
        return self.filetype=='xsf' and 'data' in self
    #end def has_data


    def validity_checks(self):
        ha = self.has_animation() 
        hb = self.has_bands()
        hs = self.has_structure()
        hd = self.has_data()
        v = ha or hb or hs or hd
        if v:
            return []
        else:
            return ['xsf file must have animation, bands, structure, or data\nthe current file is missing all of these']
        #end if
    #end def validity_checks


    # test needed
    def incorporate_structure(self,structure):
        s = structure.copy()
        s.change_units('A')
        s.recenter()
        elem = []
        for e in s.elem:
            ne = len(e)
            if ne>1:
                if ne==2 and not e[1].isalpha():
                    e = e[0]
                elif ne>2:
                    e = e[0:2]
                #end if
            #end if
            if is_element(e):
                elem.append(ptable.elements[e].atomic_number)
            else:
                elem.append(0)
            #end if
        #end for
        self.filetype    = 'xsf'
        self.periodicity = 'crystal' # assumed
        self.primvec     = s.axes
        self.elem        = array(elem,dtype=int)
        self.pos         = s.pos
    #end def incorporate_structure


    def add_density(self,cell,density,name='density',corner=None,grid=None,centered=False,add_ghost=False):
        if corner is None:
            corner = zeros((3,),dtype=float)
        #end if
        if grid is None:
            grid = density.shape
        #end if
        grid    = array(grid,dtype=int)
        corner  = array(corner,dtype=float)
        cell    = array(cell  ,dtype=float)
        density = array(density,dtype=float)
        density.shape = tuple(grid)
        
        if centered: # shift corner by half a grid cell to center it
            dc = 0.5/grid     
            dc = dot(dc,cell)
            corner += dc
        #end if

        if add_ghost: # add ghost points to make a 'general' xsf grid
            g = grid  # this is an extra shell of points in PBC
            d = density
            grid = g+1
            density = zeros(tuple(grid),dtype=float)
            density[:g[0],:g[1],:g[2]] = d[:,:,:] # volume copy
            density[   -1,:g[1],:g[2]] = d[0,:,:] # face copies
            density[:g[0],   -1,:g[2]] = d[:,0,:]
            density[:g[0],:g[1],   -1] = d[:,:,0]
            density[   -1,   -1,:g[2]] = d[0,0,:] # edge copies
            density[   -1,:g[1],   -1] = d[0,:,0] 
            density[:g[0],   -1,   -1] = d[:,0,0] 
            density[   -1,   -1,   -1] = d[0,0,0] # corner copy
            density.shape = tuple(grid)
        #end if

        self.data = obj()     
        self.data[3] = obj()
        self.data[3][name] = obj()
        self.data[3][name][name] = obj(
            grid   = grid,
            corner = corner,
            cell   = cell,
            values = density
            )
    #end def add_density


    def get_density(self):
        return self.data.first().first().first()
    #end def get_density


    # test needed
    def change_units(self,in_unit,out_unit):
        fac = 1.0/convert(1.0,in_unit,out_unit)**3
        density = self.get_density()
        density.values *= fac
        if 'values_noghost' in density:
            density.values_noghost *= fac
        #end if
    #end def change_units


    # test needed
    def remove_ghost(self,density=None):
        if density is None:
            density = self.get_density()
        #end if
        if 'values_noghost' in density:
            return density.values_noghost
        #end if
        data = density.values

        # remove the ghost cells
        d = data
        g = array(d.shape,dtype=int)-1
        data = zeros(tuple(g),dtype=float)
        data[:,:,:] = d[:g[0],:g[1],:g[2]]
        density.values_noghost = data
        return data
    #end def remove_ghost

    
    # test needed
    def norm(self,density=None,vnorm=True):
        if density is None:
            density = self.get_density()
        #end if
        if 'values_noghost' not in density:
            self.remove_ghost(density)
        #end if
        data = density.values_noghost
        if vnorm:
            dV = det(density.cell)/data.size
        else:
            dV = 1.0
        #end if
        return data.ravel().sum()*dV
    #end def norm


    # test needed
    def line_data(self,dim,density=None):
        if density is None:
            density = self.get_density()
        #end if
        if 'values_noghost' not in density:
            self.remove_ghost(density)
        #end if
        data = density.values_noghost
        dV = det(density.cell)/data.size
        dr = norm(density.cell[dim])/data.shape[dim]
        ndim = 3
        permute = dim!=0
        if permute:
            r = list(range(0,ndim))
            r.pop(dim)
            permutation = tuple([dim]+r)
            data = data.transpose(permutation)
        #end if
        s = data.shape
        data.shape = s[0],s[1]*s[2]
        line_data = data.sum(1)*dV/dr
        r_data = density.corner[dim] + dr*arange(len(line_data),dtype=float)
        return r_data,line_data
    #end def line_data


    def line_plot(self,dim,filepath):
        r,d = self.line_data(dim)
        savetxt(filepath,array(list(zip(r,d))))
    #end def line_plot

    # test needed
    def interpolate_plane(self,r1,r2,r3,density=None,meshsize=50,fill_value=0):
        if density is None:
            density = self.get_density()
        #end if
   
        dens_values = np.array(density.values)

        # Construct crystal meshgrid for dens
        da = 1./(density.grid[0]-1)
        db = 1./(density.grid[1]-1)
        dc = 1./(density.grid[2]-1)
    
        cry_corner = np.matmul(density.corner,np.linalg.inv(density.cell))
        a0  = cry_corner[0]
        b0  = cry_corner[1]
        c0  = cry_corner[2]
        
        ra = np.arange(a0, density.grid[0]*da, da)
        rb = np.arange(b0, density.grid[1]*db, db)
        rc = np.arange(c0, density.grid[2]*dc, dc)
    
        [mra, mrb, mrc] = np.meshgrid(ra, rb, rc)
    
        # 3d Interpolation on crystal coordinates
        from scipy.interpolate import RegularGridInterpolator
        g = RegularGridInterpolator((ra,rb,rc), dens_values, bounds_error=False,fill_value=fill_value)
    
        # Construct cartesian meshgrid for dens
        mrx,mry,mrz = np.array([mra,mrb,mrc]).T.dot(density.cell).T
     
        # First construct a basis (x'^,y'^,z'^) where z'^ is normal to the plane formed from ra, rb, and rc
        zph = np.cross((r2-r3),(r1-r3))
        zph = zph/np.linalg.norm(zph)
        yph = r2-r3
        yph = yph/np.linalg.norm(yph)
        xph = np.cross(yph,zph)
    
        # Positions in (x'^,y'^,z'^) basis
        rp1 = np.dot(r1,np.linalg.inv((xph,yph,zph))) 
        rp2 = np.dot(r2,np.linalg.inv((xph,yph,zph))) 
        rp3 = np.dot(r3,np.linalg.inv((xph,yph,zph)))  
    
        # Meshgrid in (x'^,y'^,z'^) basis
        mrxp,mryp,mrzp = np.array([mrx,mry,mrz]).T.dot(np.linalg.inv([xph,yph,zph])).T
    
        # Generate mesh in (x'^,y'^,z'^) basis. Ensure all points are in cell.
        xp_min = np.amin(mrxp)
        xp_max = np.amax(mrxp)
        yp_min = np.amin(mryp)
        yp_max = np.amax(mryp)
    
    
        rpx = np.arange(xp_min,xp_max,(xp_max-xp_min)/meshsize)
        rpy = np.arange(yp_min,yp_max,(yp_max-yp_min)/meshsize)
        mrpx, mrpy = np.meshgrid(rpx,rpy)
    
        slice_dens = []
        for xpi in np.arange(xp_min,xp_max,(xp_max-xp_min)/meshsize):
            yline = []
            for ypi in np.arange(yp_min,yp_max,(yp_max-yp_min)/meshsize):
                # xpi,ypi,rp1[2] to crystal coords
                rcry = np.matmul( np.dot((xpi,ypi,rp1[2]),(xph,yph,zph)) , np.linalg.inv(density.cell))
                yline.extend(g(rcry))
                #end if
            #end for
            slice_dens.append(yline)
        #end for
        slice_dens = np.array(slice_dens).T
       
        # return the following...
        # slice_dens: density on slice
        # mrpx, mrpy: meshgrid for x',y' coordinates parallel to slice, i.e., (x'^,y'^) basis
        # rp1, rp2, rp3: Input positions in (x'^,y'^,z'^) basis
        return slice_dens, mrpx, mrpy, rp1, rp2, rp3
    #end def coordinatesToSlice
#end class XsfFile



class PoscarFile(StandardFile):

    sftype = 'POSCAR'

    def __init__(self,filepath=None):
        self.description = None
        self.scale       = None
        self.axes        = None
        self.elem        = None
        self.elem_count  = None
        self.coord       = None
        self.pos         = None
        self.dynamic     = None
        self.vel_coord   = None
        self.vel         = None
        StandardFile.__init__(self,filepath)
    #end def __init__


    def assign_defaults(self):
        if self.description is None:
            self.description = 'System cell and coordinates'
        #end if
    #end def assign_defaults


    def validity_checks(self):
        msgs = []
        if self.description is None:
            msgs.append('description is missing')
        elif not isinstance(self.description,str):
            msgs.append('description must be text')
        #end if
        if self.scale is None:
            msgs.append('scale is missing')
        elif not isinstance(self.scale,(float,int)):
            msgs.append('scale must be a real number')
        elif self.scale<0:
            msgs.append('scale must be greater than zero')
        #end if
        if self.axes is None:
            msgs.append('axes is missing')
        elif not isinstance(self.axes,ndarray):
            msgs.append('axes must be an array')
        elif self.axes.shape!=(3,3):
            msgs.append('axes must be a 3x3 array, shape provided is {0}'.format(self.axes.shape))
        elif not isinstance(self.axes[0,0],float):
            msgs.append('axes must be an array of real numbers')
        #end if
        natoms = -1
        if self.elem_count is None:
            msgs.append('elem_count is missing')
        elif not isinstance(self.elem_count,ndarray):
            msgs.append('elem_count must be an array')
        elif len(self.elem_count)==0:
            msgs.append('elem_count array must contain at least one entry')
        elif not isinstance(self.elem_count[0],(int,np.int_)):
            msgs.append('elem_count must be an array of integers')
        else:
            if (self.elem_count<1).sum()>0:
                msgs.append('all elem_count entries must be greater than zero')
            #end if
            natoms = self.elem_count.sum()
        #end if
        if self.elem is not None: # presence depends on vasp version
            if not isinstance(self.elem,ndarray):
                msgs.append('elem must be an array')
            elif isinstance(self.elem_count,ndarray) and len(self.elem)!=len(self.elem_count):
                msgs.append('elem and elem_count arrays must be the same length')
            elif not isinstance(self.elem[0],str):
                msgs.append('elem must be an array of text')
            else:
                for e in self.elem:
                    iselem,symbol = is_element(e,symbol=True)
                    if not iselem:
                        msgs.append('elem entry "{0}" is not an element'.format(e))
                    #end if
                #end for
            #end for
        #end if
        if self.coord is None:
            msgs.append('coord is missing')
        elif not isinstance(self.coord,str):
            msgs.append('coord must be text')
        #end if
        if self.pos is None:
            msgs.append('pos is missing')
        elif not isinstance(self.pos,ndarray):
            msgs.append('pos must be an array')
        elif natoms>0 and self.pos.shape!=(natoms,3):
            msgs.append('pos must be a {0}x3 array, shape provided is {1}'.format(natoms),self.pos.shape)
        elif natoms>0 and not isinstance(self.pos[0,0],float):
            msgs.append('pos must be an array of real numbers')
        #end if
        if self.dynamic is not None: # dynamic is optional
            if not isinstance(self.dynamic,ndarray):
                msgs.append('dynamic must be an array')
            elif natoms>0 and self.dynamic.shape!=(natoms,3):
                msgs.append('dynamic must be a {0}x3 array, shape provided is {1}'.format(natoms),self.dynamic.shape)
            elif natoms>0 and not isinstance(self.dynamic[0,0],bool):
                msgs.append('dynamic must be an array of booleans (true/false)')
            #end if
        #end if
        if self.vel_coord is not None: # velocities are optional
            if not isinstance(self.vel_coord,str):
                msgs.append('vel_coord must be text')
            #end if
        #end if
        if self.vel is not None: # velocities are optional
            if not isinstance(self.vel,ndarray):
                msgs.append('vel must be an array')
            elif natoms>0 and self.vel.shape!=(natoms,3):
                msgs.append('vel must be a {0}x3 array, shape provided is {1}'.format(natoms),self.vel.shape)
            elif natoms>0 and not isinstance(self.vel[0,0],float):
                msgs.append('vel must be an array of real numbers')
            #end if
        #end if
        return msgs
    #end def validity_checks


    def read_text(self,text):
        read_poscar_chgcar(self,text)
    #end def read_text


    def write_text(self):
        text = ''
        if self.description is None:
            text += 'System cell and coordinates\n'
        else:
            text += self.description+'\n'
        #end if
        text += ' {0}\n'.format(self.scale)
        for a in self.axes:
            text += ' {0:20.14f} {1:20.14f} {2:20.14f}\n'.format(*a)
        #end for
        if self.elem is not None:
            for e in self.elem:
                iselem,symbol = is_element(e,symbol=True)
                if not iselem:
                    self.error('{0} is not an element'.format(e))
                #end if
                text += e+' '
            #end for
            text += '\n'
        #end if
        for ec in self.elem_count:
            text += ' {0}'.format(ec)
        #end for
        text += '\n'
        if self.dynamic!=None:
            text += 'selective dynamics\n'
        #end if
        text += self.coord+'\n'
        if self.dynamic is None:
            for p in self.pos:
                text += ' {0:20.14f} {1:20.14f} {2:20.14f}\n'.format(*p)
            #end for
        else:
            bm = self.bool_map
            for i in range(len(self.pos)):
                p = self.pos[i]
                d = self.dynamic[i]
                text += ' {0:20.14f} {1:20.14f} {2:20.14f}  {3}  {4}  {5}\n'.format(p[0],p[1],p[2],bm[d[0]],bm[d[1]],bm[d[2]])
            #end for
        #end if
        if self.vel!=None:
            text += self.vel_coord+'\n'
            for v in self.vel:
                text += ' {0:20.14f} {1:20.14f} {2:20.14f}\n'.format(*v)
            #end for
        #end if
        return text
    #end def write_text


    def incorporate_xsf(self,xsf):
        if 'primvec' in xsf:
            axes = xsf.primvec.copy()
        #end if
        if 'convvec' in xsf:
            axes = xsf.convvec.copy()
        #end if
        elem = xsf.elem.copy()
        pos  = xsf.pos.copy()

        species        = []
        species_counts = []
        elem_indices   = []

        spec_set = set()
        for i in range(len(elem)):
            e = elem[i]
            if not e in spec_set:
                spec_set.add(e)
                species.append(e)
                species_counts.append(0)
                elem_indices.append([])
            #end if
            sindex = species.index(e)
            species_counts[sindex] += 1
            elem_indices[sindex].append(i)
        #end for

        elem_order = []
        for elem_inds in elem_indices:
            elem_order.extend(elem_inds)
        #end for

        pos = pos[elem_order]

        species_ind = species
        species = []
        for i in species_ind:
            species.append(ptable.simple_elements[i].symbol)
        #end for

        self.scale      = 1.0
        self.axes       = axes
        self.elem       = array(species,dtype=str)
        self.elem_count = array(species_counts,dtype=int)
        self.coord      = 'cartesian'
        self.pos        = pos

        self.assign_defaults()
    #end def incorporate_xsf
#end class PoscarFile



class ChgcarFile(StandardFile):

    sftype = 'CHGCAR'

    def __init__(self,filepath=None):
        self.poscar         = None
        self.grid           = None
        self.charge_density = None
        self.spin_density   = None
        StandardFile.__init__(self,filepath)
    #end def __init__


    def validity_checks(self):
        msgs = []
        if self.poscar is None:
            msgs.append('poscar elements are missing')
        elif not isinstance(self.poscar,PoscarFile):
            msgs.append('poscar is not an instance of PoscarFile')
        else:
            msgs.extend(self.poscar.validity_checks())
        #end if
        if self.grid is None:
            msgs.append('grid is missing')
        elif not isinstance(self.grid,ndarray):
            msgs.append('grid must be an array')
        elif len(self.grid)!=3 or self.grid.size!=3:
            msgs.append('grid must have 3 entries')
        elif not isinstance(self.grid[0],(int,np.int_)):
            msgs.append('grid must be an array of integers')
        elif (self.grid<1).sum()>0:
            msgs.append('all grid entries must be greater than zero')
        #end if
        if self.grid is not None:
            ng = self.grid.prod()
        #end if
        if self.charge_density is None:
            msgs.append('charge_density is missing')
        elif not isinstance(self.charge_density,ndarray):
            msgs.append('charge_density must be an array')
        elif len(self.charge_density)!=ng:
            msgs.append('charge_density must have {0} entries ({1} present by length)'.format(ng,len(self.charge_density)))
        elif self.charge_density.size!=ng:
            msgs.append('charge_density must have {0} entries ({1} present by size)'.format(ng,self.charge_density.size))
        elif not isinstance(self.charge_density[0],float):
            msgs.append('charge_density must be an array of real numbers')
        #end if
        if self.spin_density is not None: # spin density is optional
            if not isinstance(self.spin_density,ndarray):
                msgs.append('spin_density must be an array')
            elif len(self.spin_density)!=ng:
                msgs.append('spin_density must have {0} entries ({1} present)'.format(ng,len(self.spin_density)))
            elif self.spin_density.size!=ng and self.spin_density.shape!=(ng,3):
                msgs.append('non-collinear spin_density must be a {0}x3 array, shape provided: {1}'.format(ng,self.spin_density.shape))
            elif not isinstance(self.spin_density.ravel()[0],float):
                msgs.append('spin_density must be an array of real numbers')
            #end if
        #end if
        return msgs
    #end def validity_checks


    def read_text(self,text):
        read_poscar_chgcar(self,text)
    #end def read_text


    def write_text(self):
        text = self.poscar.write_text()
        text+= '\n {0} {1} {2}\n'.format(*self.grid)
        densities = [self.charge_density]
        if self.spin_density is not None:
            if self.spin_density.size==self.charge_density.size:
                densities.append(self.spin_density)
            else:
                for i in range(3):
                    densities.append(self.spin_density[:,i])
                #end for
            #end if
        #end if
        n=0
        for dens in densities:
            for d in dens:
                text += '{0:20.12E}'.format(d)
                n+=1
                if n%5==0:
                    text+='\n'
                #end if
            #end for
        #end for
        return text
    #end def write_text


    def incorporate_xsf(self,xsf):
        poscar = PoscarFile()
        poscar.incorporate_xsf(xsf)
        density = xsf.remove_ghost().copy()
        self.poscar         = poscar
        self.grid           = array(density.shape,dtype=int)
        self.charge_density = density.ravel(order='F')
        self.check_valid()
    #end def incorporate_xsf
#end class ChgcarFile



def read_poscar_chgcar(host,text):
    is_poscar = isinstance(host,PoscarFile)
    is_chgcar = isinstance(host,ChgcarFile)
    if not is_poscar and not is_chgcar:
        error('read_poscar_chgcar must be used in conjunction with PoscarFile or ChgcarFile objects only\nencountered object of type: {0}'.format(host.__class__.__name__))
    #end if

    # read lines and remove fortran comments
    raw_lines = text.splitlines()
    lines = []
    for line in raw_lines:
        # remove fortran comments
        cloc1 = line.find('!')
        cloc2 = line.find('#')
        has1  = cloc1!=-1
        has2  = cloc2!=-1
        if has1 or has2:
            if has1 and has2:
                cloc = min(cloc1,cloc2)
            elif has1:
                cloc = cloc1
            else:
                cloc = cloc2
            #end if
            line = line[:cloc]
        #end if
        lines.append(line.strip())
    #end for

    # extract file information
    nlines = len(lines)
    min_lines = 8
    if nlines<min_lines:
        host.error('file {0} must have at least {1} lines\nonly {2} lines found'.format(host.filepath,min_lines,nlines))
    #end if
    description = lines[0]
    dim = 3
    scale = float(lines[1].strip())
    axes = empty((dim,dim))
    axes[0] = array(lines[2].split(),dtype=float)
    axes[1] = array(lines[3].split(),dtype=float)
    axes[2] = array(lines[4].split(),dtype=float)
    tokens = lines[5].split()
    if tokens[0].isdigit():
        counts = array(tokens,dtype=int)
        elem   = None
        lcur   = 6
    else:
        elem   = array(tokens,dtype=str)
        counts = array(lines[6].split(),dtype=int)
        lcur   = 7
    #end if

    if lcur<len(lines) and len(lines[lcur])>0:
        c = lines[lcur].lower()[0]
        lcur+=1
    else:
        host.error('file {0} is incomplete (missing positions)'.format(host.filepath))
    #end if
    selective_dynamics = c=='s'
    if selective_dynamics: # Selective dynamics
        if lcur<len(lines) and len(lines[lcur])>0:
            c = lines[lcur].lower()[0]
            lcur+=1
        else:
            host.error('file {0} is incomplete (missing positions)'.format(host.filepath))
        #end if
    #end if
    cartesian = c=='c' or c=='k'
    if cartesian:
        coord = 'cartesian'
    else:
        coord = 'direct'
    #end if
    npos = counts.sum()
    if lcur+npos>len(lines):
        host.error('file {0} is incomplete (missing positions)'.format(host.filepath))
    #end if
    spos = []
    for i in range(npos):
        spos.append(lines[lcur+i].split())
    #end for
    lcur += npos
    spos = array(spos)
    pos  = array(spos[:,0:3],dtype=float)
    if selective_dynamics:
        dynamic = array(spos[:,3:6],dtype=str)
        dynamic = dynamic=='T'
    else:
        dynamic = None
    #end if

    def is_empty(lines,start=None,end=None):
        if start is None:
            start = 0
        #end if
        if end is None:
            end = len(lines)
        #end if
        is_empty = True
        for line in lines[start:end]:
            is_empty &= len(line)==0
        #end for
        return is_empty
    #end def is_empty

    # velocities may be present for poscar
    #   assume they are not for chgcar
    if is_poscar and lcur<len(lines) and not is_empty(lines,lcur):
        cline = lines[lcur].lower()
        lcur+=1
        if lcur+npos>len(lines):
            host.error('file {0} is incomplete (missing velocities)'.format(host.filepath))
        #end if
        cartesian = len(cline)>0 and (cline[0]=='c' or cline[0]=='k')
        if cartesian:
            vel_coord = 'cartesian'
        else:
            vel_coord = 'direct'
        #end if
        svel = []
        for i in range(npos):
            svel.append(lines[lcur+i].split())
        #end for
        lcur += npos
        vel = array(svel,dtype=float)
    else:
        vel_coord = None
        vel = None
    #end if

    # grid data is present for chgcar
    if is_chgcar:
        lcur+=1
        if lcur<len(lines) and len(lines[lcur])>0:
            grid = array(lines[lcur].split(),dtype=int)
            lcur+=1
        else:
            host.error('file {0} is incomplete (missing grid)'.format(host.filepath))
        #end if
        if lcur<len(lines):
            ng = grid.prod()
            density = []
            for line in lines[lcur:]:
                density.extend(line.split())
            #end for
            if len(density)>0:
                def is_float(val):
                    try:
                        v = float(val)
                        return True
                    except:
                        return False
                    #end try
                #end def is_float
                # remove anything but the densities (e.g. augmentation charges)
                n=0
                while is_float(density[n]):
                    n+=ng
                    if n+ng>=len(density):
                        break
                    #end if
                #end while
                density = array(density[:n],dtype=float)
            else:
                host.error('file {0} is incomplete (missing density)'.format(host.filepath))
            #end if
            if density.size%ng!=0:
                host.error('number of density data entries is not a multiple of the grid\ngrid shape: {0}\ngrid size: {1}\ndensity size: {2}'.format(grid,ng,density.size))
            #end if
            ndens = density.size//ng
            if ndens==1:
                charge_density = density
                spin_density   = None
            elif ndens==2:
                charge_density = density[:ng]
                spin_density   = density[ng:]
            elif ndens==4:
                charge_density = density[:ng]
                spin_density   = empty((ng,3),dtype=float)
                for i in range(3):
                    spin_density[:,i] = density[(i+1)*ng:(i+2)*ng]
                #end for
            else:
                host.error('density data must be present for one of the following situations\n  1) charge density only (1 density)\n  2) charge and collinear spin densities (2 densities)\n  3) charge and non-collinear spin densities (4 densities)\nnumber of densities found: {0}'.format(ndens))
            #end if
        else:
            host.error('file {0} is incomplete (missing density)'.format(host.filepath))
        #end if
    #end if

    if is_poscar:
        poscar = host
    elif is_chgcar:
        poscar = PoscarFile()
    #end if

    poscar.set(
        description = description,
        scale       = scale,
        axes        = axes,
        elem        = elem,
        elem_count  = counts,
        coord       = coord,
        pos         = pos,
        dynamic     = dynamic,
        vel_coord   = vel_coord,
        vel         = vel
        )

    if is_chgcar:
        host.set(
            poscar         = poscar,
            grid           = grid,
            charge_density = charge_density,
            spin_density   = spin_density,
            )
    #end if
#end def read_poscar_chgcar
