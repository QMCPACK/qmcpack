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
from numpy import array,zeros,ndarray,around,arange,dot,savetxt
from numpy.linalg import det,norm
from generic import obj
from developer import DevBase
from periodic_table import pt as ptable
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
        line = self.mm.readline()
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
            for i in xrange(len(formats)):
                tokens.append(formats[i](stokens[i]))
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
        return self.mm.readline()
    #end def readline

    def read(self,num=None):
        if num is None:
            return self.mm[:]
        else:
            return self.mm.read(num)
        #end if
    #end def read


    # unchanged mmap interface below
    def find(self,*a,**kw):
        return self.mm.find(*a,**kw)
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
        return self.mm.rfind(*a,**kw)
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



class XsfFile(DevBase):
    filetypes     = set(['xsf','axsf','bxsf'])
    periodicities = set(['molecule','polymer','slab','crystal']) 
    dimensions    = obj(molecule=0,polymer=1,slab=2,crystal=3)

    # ATOMS  are in units of Angstrom, only provided for 'molecule'
    # forces are in units of Hatree/Angstrom
    # each section should be followed by a blank line

    def __init__(self,arg0=None):
        self.filetype    = None
        self.periodicity = None
        # primvec convec elem pos force data
        # animsteps images
        # band info
        if arg0!=None:
            if isinstance(arg0,str):
                filepath = arg0
                self.read(filepath)
            else:
                self.error('unsupported input: {0}'.format(arg0))
            #end if
        #end if
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


    def read(self,filepath):
        if not os.path.exists(filepath):
            self.error('read failed, file does not exist: {0}'.format(filepath))
        #end if
        self.read_text(open(filepath,'r').read(),check=False)
        if not self.is_valid():
            self.error('read failed,file {0} is not a valid xsf file'.format(filepath))
        #end if
    #end def read


    def write(self,filepath=None):
        text = self.write_text()
        if filepath!=None:
            open(filepath,'w').write(text)
        #end if
        return text
    #end def write


    def read_text(self,text,check=True):
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
                            grid_data.shape = tuple(grid)
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
                            for bi,bv in bands.iteritems():
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
        if check and not self.is_valid():
            self.error('read failed, not a valid xsf file')
        #end if
    #end def read_text


    def write_text(self):
        if not self.is_valid():
            self.error('write failed, not a valid xsf file')
        #end if
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
                    for v in dg.values.ravel():
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


    def is_valid(self):
        ha = self.has_animation() 
        hb = self.has_bands()
        hs = self.has_structure()
        hd = self.has_data()
        v = ha or hb or hs or hd
        return v
    #end def is_valid


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
            elem.append(ptable.elements[e].atomic_number)
        #end for
        self.filetype    = 'xsf'
        self.periodicity = 'crystal' # assumed
        self.primvec     = s.axes
        self.elem        = array(elem,dtype=int)
        self.pos         = s.pos
    #end def incorporate_structure


    def add_density(self,cell,density,name='density',corner=None,grid=None,centered=False,add_ghost=False,transpose=False):
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
        #end if

        if transpose: # shift from row major to column major
            g = grid
            d = density
            density = zeros((d.size,))
            n = 0
            for k in xrange(g[2]):
                for j in xrange(g[1]):
                    for i in xrange(g[0]):
                        density[n] = d[i,j,k]
                        n+=1
                    #end for
                #end for
            #end for
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


    def change_units(self,in_unit,out_unit):
        fac = 1.0/convert(1.0,in_unit,out_unit)**3
        density = self.get_density()
        density.values *= fac
        if 'values_noghost' in density:
            density.values_noghost *= fac
        #end if
    #end def change_units


    def remove_ghost(self,density=None,transpose=True):
        if density is None:
            density = self.get_density()
        #end if
        if 'values_noghost' in density:
            return density.values_noghost
        #end if
        data = density.values
        if transpose: # switch from column major to row major
            g = data.shape
            d = data.ravel()
            data = zeros(g,dtype=float)
            n = 0
            for k in xrange(g[2]):
                for j in xrange(g[1]):
                    for i in xrange(g[0]):
                        data[i,j,k] = d[n]
                        n+=1
                    #end for
                #end for
            #end for
        #end if
        # remove the ghost cells
        d = data
        g = array(d.shape,dtype=int)-1
        data = zeros(tuple(g),dtype=float)
        data[:,:,:] = d[:g[0],:g[1],:g[2]]
        density.values_noghost = data
        return data
    #end def remove_ghost

    
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
            r = range(0,ndim)
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
        savetxt(filepath,array(zip(r,d)))
    #end def line_plot
#end class XsfFile
