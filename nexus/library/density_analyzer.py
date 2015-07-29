##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  density_analyzer.py                                               #
#    Support for manipulation of gridded density files.              #
#                                                                    #
#  Content summary:                                                  #
#    DensityAnalyzer                                                 #
#      Class to represent a density file.                            #
#      Supports density addition, subtraction, integration, etc.     #
#      Currently reads xsf files.                                    #
#                                                                    #                                        
#====================================================================#


import os
from numpy import array,zeros,ndarray,around,sqrt,arange
from numpy.linalg import norm
from developer import DevBase
from generic import obj
from structure import Structure
from plotting import *


class DensityAnalyzer(DevBase):
    def __init__(self,densfile=None,denserrfile=None,transpose=False):
        density       = None
        density_error = None
        corner        = None
        cell          = None
        grid          = None
        structure     = None
        if densfile!=None:
            self.read(densfile,denserrfile,transpose)
        #end if
    #end def __init__


    def read(self,densfile,denserrfile=None,transpose=False):
        if not os.path.exists(densfile):
            self.error('density file {0} does not exist'.format(densfile))
        #end if
        ext = densfile.rsplit('.',1)[1]
        if denserrfile!=None:
            if not os.path.exists(denserrfile):
                self.error('density error file {0} does not exist'.format(denserrfile))
            elif densfile.rsplit('.',1)[1]!=ext:
                self.error('file extensions must match\n  density file extension: {0}\n  density error file extension: {1}'.format(ext,densfile.rsplit('.',1)[1]))
        #end if
        if ext=='xsf':
            self.read_xsf(densfile,transpose)
            if denserrfile!=None:
                self.read_xsf(denserrfile,transpose,error=True)
            #end if
        else:
            self.error('file format {0} is not yet supported'.format(ext))
        #end if
        if tuple(self.grid)!=self.density.shape:
            self.error('density shape does not match grid\n  density.shape = {0}\n  grid = {1}'.format(seld.density.shape,self.grid))
        #end if
    #end def read


    def read_xsf(self,densfile,transpose=False,error=False):
        primvec = None
        natoms  = None
        elem    = None
        pos     = None
        grid    = None
        corner  = None
        cell    = None
        dens    = None
        lines = open(densfile,'r').read().splitlines()
        i=0
        while(i<len(lines)):
            line = lines[i].strip()
            if line=='PRIMVEC':
                primvec = array((lines[i+1]+' '+
                                 lines[i+2]+' '+
                                 lines[i+3]).split(),dtype=float)
                primvec.shape = 3,3
                i+=3
            elif line=='PRIMCOORD':
                natoms = int(lines[i+1].split()[0])
                elem = []
                pos  = []
                for iat in range(natoms):
                    tokens = lines[i+2+iat].split()
                    elem.append(tokens[0])
                    pos.extend(tokens[1:])
                #end for
                pos = array(pos,dtype=float)
                pos.shape = natoms,3
                elem = elem
                pos  = pos
                i+=natoms+1
            elif line.startswith('DATAGRID_3D'):
                grid   = array(lines[i+1].split(),dtype=int)
                corner = array(lines[i+2].split(),dtype=float)
                cell   = array((lines[i+3]+' '+
                                lines[i+4]+' '+
                                lines[i+5]).split(),dtype=float)
                cell.shape = 3,3
                i+=6
                dstr = ''
                line = lines[i].strip()
                while line!='END_DATAGRID_3D':
                    dstr += line+' '
                    i+=1
                    line = lines[i].strip()
                #end while
                dens = array(dstr.split(),dtype=float)
                if transpose:
                    dtrans = zeros(dens.shape)
                    gdims = grid.copy()
                    gdims[0] = grid[1]*grid[2]
                    gdims[1] = grid[2]
                    gdims[2] = 1
                    nd=0
                    for k in xrange(grid[2]):
                        for j in xrange(grid[1]):
                            for i in xrange(grid[0]):
                                p = i*gdims[0]+j*gdims[1]+k*gdims[2]
                                dtrans[nd] = dens[p]
                            #end for
                        #end for
                    #end for
                    dens = dtrans
                #end if
                dens.shape = tuple(grid)
            #end if
            i+=1
        #end while
        v=obj(primvec=primvec,natoms=natoms,elem=elem,pos=pos,
              grid=grid,corner=corner,cell=cell,dens=dens)
        for name,val in v.iteritems():
            if val is None:
                self.error(name+' not found in xsf file '+densfile)
            #end if
        #end for
        if not error:
            self.set(
                corner  = corner,
                cell    = cell,
                grid    = grid,
                density = dens,
                structure = Structure(
                    axes  = primvec,
                    elem  = elem,
                    pos   = pos,
                    units = 'A'
                    )
                )
        else:
            self.density_error = dens
        #end if
    #end def read_xsf


    def sum(self):
        return self.density.sum()
    #end def sum

    
    def normalize(self,norm):
        f = norm/self.sum()
        self.density*=f
        if self.density_error!=None:
            self.density_error*=f
        #end if
    #end def normalize


    def assert_shape(self,other):
        if self.density.shape!=other.density.shape:
            self.error('two density grids are not the same shape')
        #end if
    #end def assert_shape


    def __add__(self,other):
        self.assert_shape(other)
        s = self.copy()
        s.density += other.density
        if s.density_error!=None and other.density_error!=None:
            s.density_error = sqrt(s.density_error**2+other.density_error**2)
        #end if
        return s
    #end def __add__


    def __sub__(self,other):
        self.assert_shape(other)
        s = self.copy()
        s.density -= other.density
        if s.density_error!=None and other.density_error!=None:
            s.density_error = sqrt(s.density_error**2+other.density_error**2)
        #end if
        return s
    #end def __sub__


    def integrate(self,xr=(0,1),yr=(0,1),zr=(0,1)):
        nx,ny,nz = self.grid
        xr = array(around(array(xr,dtype=float)*nx),dtype=int)
        yr = array(around(array(yr,dtype=float)*ny),dtype=int)
        zr = array(around(array(zr,dtype=float)*nz),dtype=int)
        ds = self.density[xr[0]:xr[1],yr[0]:yr[1],zr[0]:zr[1]].sum()
        if self.density_error is None:
            return ds
        else:
            dse = sqrt((self.density_error[xr[0]:xr[1],yr[0]:yr[1],zr[0]:zr[1]]**2).sum())
            return ds,dse
        #end if
    #end def integrate


    def project_line(self,xr=None,yr=None,zr=None):
        if (xr==None) + (yr==None) + (zr==None) != 1:
            self.error('to project onto a line, two ranges must be given to sum over')
        #end if
        if xr is None:
            d  = 0
            xr = 0,1
        elif yr is None:
            d  = 1
            yr = 0,1
        elif zr is None:
            d  = 2
            zr = 0,1
        #end if
        r = norm(self.cell[d])/self.grid[d]*arange(self.grid[d])
        ds = [0,1,2]
        ds.pop(d)
        ds.reverse()
        ds = tuple(ds)
        nx,ny,nz = self.grid
        xr = array(around(array(xr,dtype=float)*nx),dtype=int)
        yr = array(around(array(yr,dtype=float)*ny),dtype=int)
        zr = array(around(array(zr,dtype=float)*nz),dtype=int)
        lp = self.density[xr[0]:xr[1],yr[0]:yr[1],zr[0]:zr[1]]
        for d in ds:
            lp = lp.sum(axis=d)
        #end for
        if self.density_error is None:
            return r,lp
        else:
            lpe = self.density_error[xr[0]:xr[1],yr[0]:yr[1],zr[0]:zr[1]]**2
            for d in ds:
                lpe = lpe.sum(axis=d)
            #end for
            lpe = sqrt(lpe)
            return r,lp,lpe
        #end if
    #end def project_line


    def plot_line(self,xr=None,yr=None,zr=None,fmt='b.-',label='line dens',disp=False,labels=False):
        labels = labels or disp
        if disp:
            figure()
        #end if
        if self.density_error is None:
            r,lp = self.project_line(xr,yr,zr)
            plot(r,lp,fmt,label=label)
        else:
            r,lp,lpe = self.project_line(xr,yr,zr)
            errorbar(r,lp,lpe,fmt=fmt,label=label)
        #end if
        if labels:
            xlabel('r (Angstrom)')
            ylabel('density')
            title('Line projection of density')
        #end if
        if disp:
            show()
        #end if
    #end def plot_line
#end class DensityAnalyzer

