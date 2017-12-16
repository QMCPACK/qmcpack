##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  structure.py                                                      #
#    Support for atomic structure I/O, generation, and manipulation. #
#                                                                    #
#  Content summary:                                                  #
#    Structure                                                       #
#      Represents a simulation cell containing a set of atoms.       #
#      Many functions for manipulating structures or obtaining       #
#        data regarding local atomic structure.                      #
#                                                                    #
#    generate_cell                                                   #
#      User-facing function to generate an empty simulation cell.    #
#                                                                    #
#    generate_structure                                              #
#      User-facing function to specify arbitrary atomic structures   #
#      or generate structures corresponding to atoms, dimers, or     #
#      crystals.                                                     #
#                                                                    #
#====================================================================#


#! /usr/bin/env python

import os
from copy import deepcopy
from random import randint
from numpy import array,floor,empty,dot,diag,sqrt,pi,mgrid,exp,append,arange,ceil,cross,cos,sin,identity,ndarray,atleast_2d,around,ones,zeros,logical_not,flipud,uint64,sign
from numpy.linalg import inv,det,norm
from types import NoneType
from unit_converter import convert
from numerics import nearest_neighbors,convex_hull,voronoi_neighbors
from periodic_table import pt,is_element
from fileio import XsfFile
from generic import obj
from developer import DevBase,unavailable,error,warn
from debug import ci,ls,gs


try:
    from scipy.special import erfc
except ImportError:
    erfc = unavailable('scipy.special','erfc')
#end try
try:
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import plot,subplot,title,xlabel,ylabel
except (ImportError,RuntimeError):
    plot,subplot,title,xlabel,ylabel,plt = unavailable('matplotlib.pyplot','plot','subplot','title','xlabel','ylabel','plt')
#end try





# installation instructions to enable cif file read
#   
#   cif file support in Nexus currently requires two external libraries
#     PyCifRW  - base interface to read cif files into object format: CifFile
#     cif2cell - translation layer from CifFile object to cell reconstruction: CellData
#     (note: cif2cell installation includes PyCifRW)
#
#  installation of cif2cell
#    go to http://sourceforge.net/projects/cif2cell/
#    click on Download (example: cif2cell-1.2.10.tar.gz)
#    unpack directory (tar -xzf cif2cell-1.2.10.tar.gz)
#    enter directory (cd cif2cell-1.2.10)
#    install cif2cell (python setup.py install)
#    check python installation
#      >python
#      >>>from CifFile import CifFile
#      >>>from uctools import CellData
#   
#   Nexus is currently compatible with
#     cif2cell-1.2.10 and PyCifRW-3.3
#     cif2cell-1.2.7  and PyCifRW-4.1.1 
#     compatibility last tested: 20 Mar 2017
#
try:
    from CifFile import CifFile
except:
    CifFile = unavailable('CifFile','CifFile')
#end try
try:
    from uctools import CellData
except:
    CellData = unavailable('uctools','CellData')
#end try


cif2cell_unit_dict = dict(angstrom='A',bohr='B',nm='nm')



def read_cif_celldata(filepath,block=None,grammar='1.1'):
    # read cif file with PyCifRW
    path,cif_file = os.path.split(filepath)
    cwd = os.getcwd()
    os.chdir(path)
    cf = CifFile(cif_file,grammar=grammar)
    #cf = ReadCif(cif_file,grammar=grammar)
    os.chdir(cwd)
    if block is None:
        block = cf.keys()[0]
    #end if
    cb = cf.get(block)
    if cb is None:
        error('block {0} was not found in cif file {1}'.format(block,filepath),'read_cif_celldata')
    #end if

    # repack H-M symbols as normal strings so CellData.getFromCIF won't choke on unicode
    for k in ['_symmetry_space_group_name_H-M','_space_group_name_H-M_alt','_symmetry_space_group_name_h-m','_space_group_name_h-m_alt']:
        if k in cb.block:
            v = cb.block[k]
            if isinstance(v,(list,tuple)):
                for i in range(len(v)):
                    if isinstance(v[i],unicode):
                        v[i] = str(v[i])
                    #end if
                #end for
            #end if
        #end if
    #end for

    # extract structure from CifFile with uctools CellData class
    cd = CellData()
    cd.getFromCIF(cb)

    return cd
#end def read_cif_celldata



def read_cif_cell(filepath,block=None,grammar='1.1',cell='prim'):
    cd = read_cif_celldata(filepath,block,grammar)

    if cell.startswith('prim'):
        cell = cd.primitive()
    elif cell.startswith('conv'):
        cell = cd.conventional()
    else:
        error('cell argument must be primitive or conventional\nyou provided: {0}'.format(cell),'read_cif_cell')
    #end if

    return cell
#end def read_cif_cell



def read_cif(filepath,block=None,grammar='1.1',cell='prim',args_only=False):
    if isinstance(filepath,str):
        cell = read_cif_cell(filepath,block,grammar,cell)
    else:
        cell = filepath
    #end if

    # create Structure object from cell
    if cell.alloy:
        error('cannot handle alloys','read_cif')
    #end if
    units = cif2cell_unit_dict[cell.unit]
    scale = float(cell.lengthscale)
    scale = convert(scale,units,'A')
    units = 'A'
    axes  = scale*array(cell.latticevectors,dtype=float)
    elem  = []
    pos   = []
    for wyckoff_atoms in cell.atomdata:
        for atom in wyckoff_atoms:
            elem.append(str(atom.species.keys()[0]))
            pos.append(atom.position)
        #end for
    #end for
    pos = dot(array(pos,dtype=float),axes)

    if not args_only:
        s = Structure(
            axes  = axes,
            elem  = elem,
            pos   = pos,
            units = units
            )
        return s
    else:
        return axes,elem,pos,units
    #end if
#end def read_cif





# installation instructions for spglib interface
#
#  this is bootstrapped off of spglib's ASE Python interface
#
#  installation of spglib
#    go to http://sourceforge.net/projects/spglib/files/
#    click on Download spglib-1.8.2.tar.gz (952.6 kB)
#    unpack directory (tar -xzf spglib-1.8.2.tar.gz)
#    enter ase directory (cd spglib-1.8.2/python/ase/)
#    build and install (sudo python setup.py install)
from periodic_table import pt as ptable
try:
    from pyspglib import spglib
except:
    spglib = unavailable('pyspglib','spglib')
#end try





def equate(expr):
    return expr
#end def equate

def negate(expr):
    return not expr
#end def negate



def kmesh(kaxes,dim,shift=None):
    '''
    Create a Monkhorst-Pack k-point mesh 
    '''
    if shift is None:
        shift = (0.,0,0)
    #end if
    ndim = len(dim)
    d = array(dim)
    s = array(shift)
    s.shape = 1,ndim
    d.shape = 1,ndim
    kp = empty((1,ndim),dtype=float)
    kgrid = empty((d.prod(),ndim))
    n=0
    for k in range(dim[2]):
        for j in range(dim[1]):
            for i in range(dim[0]):
                kp[:] = i,j,k
                kp = dot((kp+s)/d,kaxes)
                #kp = (kp+s)/d
                kgrid[n] = kp
                n+=1
            #end for
        #end for
    #end for
    return kgrid
#end def kmesh



def reduce_tilematrix(tiling):
    tiling = array(tiling)
    t = array(tiling,dtype=int)
    if abs(tiling-t).sum()>1e-6:
        Structure.class_error('requested tiling is non-integer\n tiling requested: '+str(tiling))
    #end if

    dim = len(t)
    matrix_tiling = t.shape == (dim,dim)
    if matrix_tiling:
        #find a tiling tuple from the tiling matrix
        # do this by shearing the tiling matrix (or equivalently the tiled cell)
        # until it is orthogonal (in the untiled cell axes)
        # this is just rearranging triangular tiles of space to reshape the cell
        # so that t1*t2*t3 = det(T) = det(A_tiled)/det(A_untiled)
        #this way the atoms in the (perhaps oddly shaped) supercell can be 
        # obtained from simple translations of the untiled cell positions
        T = t  #tiling matrix
        tilematrix = T.copy()
        del t
        Tnew = array(T,dtype=float) #sheared/orthogonal tiling matrix
        tbar = identity(dim) #basis for shearing
        dr = range(dim)
        #dr = [1,0,2]
        other = dim*[0] # other[d] = dimensions other than d
        for d in dr: 
            other[d] = set(dr)-set([d])
        #end for
        #move each axis to be parallel to barred directions
        # these are volume preserving shears of the supercell
        # each shear keeps two cell face planes fixed while moving the others
        for dp in [(0,1,2),(2,0,1),(1,2,0),(2,1,0),(0,2,1),(1,0,2)]:
            success = True
            for d in dr:
                tb = tbar[dp[d]] 
                t  = T[d]
                d2,d3 = other[d]
                n = cross(Tnew[d2],Tnew[d3])  #vector normal to 2 cell faces
                vol   = dot(n,t)
                bcomp = dot(n,tb)
                if abs(bcomp)<1e-6:
                    success = False
                    break
                #end if
                tn = vol*1./bcomp*tb #new axis vector
                Tnew[d] = tn
            #end for
            if success:
                break
            #end if
        #end for
        # apply inverse permutation, if needed
        Tn = Tnew.copy()
        for d in dr:
            d2 = dp[d]
            Tnew[d2] = Tn[d]
        #end for
        #the resulting tiling matrix should be diagonal and integer
        tr = diag(Tnew)
        t  = array(around(tr),dtype=int)
        nondiagonal = abs(Tnew-diag(tr)).sum()>1e-6
        noninteger  = abs(tr-t).sum()>1e-6
        if nondiagonal:
            Structure.class_error('could not find a diagonal tiling matrix for generating tiled coordinates')
        #end if
        #non-integer tile vectors are handled directly by tile_points now
        #if noninteger:
        #    Structure.class_error('calculated diagonal tiling matrix is non-integer\n  tiled coordinates cannot be determined')
        ##end if
        #tilevector = abs(t)
        tilevector = abs(tr)
    else:
        tilevector = t
        tilematrix = diag(t)
    #end if

    return tilematrix,tilevector
#end def reduce_tilematrix


def tile_magnetization(mag,tilevec,mag_order,mag_prim):
    # jtk mark current
    #  implement true magnetic tiling based on the magnetic order
    #  Structure needs a function magnetic_period which takes magnetic order
    #   and translates it into a magnetic period tuple
    #  magnetic_period should divide evenly into the tiling vector
    #  the magnetic unit cell is self.tile(magnetic_period)
    #  re-representing the folded cell as the magnetic primitive cell
    #   requires different axes, often with a non-diagonal tiling matrix
    #  if magnetic primitive is requested, a small tiling vector should first be used
    #   (ie 221,212,122,222 periods all involve a 211 prim tiling w/ a simple reshaping/reassignment of the cell axes
    #  Structure should have a function providing labels to each magnetic species
    mag = array(int(round( tilevec.prod() ))*list(mag),dtype=object)
    return mag
#end def tile_magnetization


def rotate_plane(plane,angle,points,units='degrees'):
    if units=='degrees':
        angle *= pi/180
    elif not units.startswith('rad'):
        error('angular units must be degrees or radians\nyou provided: {0}'.format(angle),'rotate_plane')
    #end if
    c = cos(angle)
    s = sin(angle)
    if plane=='xy':
        R = [[ c,-s, 0],
             [ s, c, 0],
             [ 0, 0, 1]]
    elif plane=='yx':
        R = [[ c, s, 0],
             [-s, c, 0],
             [ 0, 0, 1]]
    elif plane=='yz':
        R = [[ 1, 0, 0],
             [ 0, c,-s],
             [ 0, s, c]]
    elif plane=='zy':
        R = [[ 1, 0, 0],
             [ 0, c, s],
             [ 0,-s, c]]
    elif plane=='zx':
        R = [[ c, 0, s],
             [ 0, 1, 0],
             [-s, 0, c]]
    elif plane=='xz':
        R = [[ c, 0,-s],
             [ 0, 1, 0],
             [ s, 0, c]]
    else:
        error('plane must be xy/yx/yz/zy/zx/xz\nyou provided: {0}'.format(plane),'rotate_plane')
    #end if
    R = array(R,dtype=float)
    return dot(R,points.T).T
#end def rotate_plane



opt_tm_matrices    = obj()
opt_tm_wig_indices = obj()

def trivial_filter(T):
    return True
#end def trival_filter

class MaskFilter(DevBase):
    def set(self,mask,dim=3):
        omask = array(mask)
        mask  = array(mask,dtype=bool)
        if mask.size==dim:
            mvec = mask.ravel()
            mask = empty((dim,dim),dtype=bool)
            i=0
            for mi in mvec:
                j=0
                for mj in mvec:
                    mask[i,j] = mi==mj
                    j+=1
                #end for
                i+=1
            #end for
        elif mask.shape!=(dim,dim):
            error('shape of mask array must be {0},{0}\nshape received: {1},{2}\nmask array received: {3}'.format(dim,mask.shape[0],mask.shape[1],omask),'optimal_tilematrix')
        #end if
        self.mask = mask==False
    #end def set

    def __call__(self,T):
        return (T[self.mask]==0).all()
    #end def __call__
#end class MaskFilter
mask_filter = MaskFilter()
            

def optimal_tilematrix(axes,volfac,dn=1,tol=1e-3,filter=trivial_filter,mask=None,nc=5):
    if mask is not None:
        mask_filter.set(mask)
        filter = mask_filter
    #end if
    dim = 3
    if isinstance(axes,Structure):
        axes = axes.axes
    else:
        axes = array(axes,dtype=float)
    #end if
    if not isinstance(volfac,int):
        volfac = int(around(volfac))
    #end if
    volume = abs(det(axes))*volfac
    axinv  = inv(axes)
    cube   = volume**(1./3)*identity(dim)
    Tref   = array(around(dot(cube,axinv)),dtype=int)
    # calculate and store all tiling matrix variations
    if dn not in opt_tm_matrices:
        mats = []
        rng = tuple(range(-dn,dn+1))
        for n1 in rng:
            for n2 in rng:
                for n3 in rng:
                    for n4 in rng:
                        for n5 in rng:
                            for n6 in rng:
                                for n7 in rng:
                                    for n8 in rng:
                                        for n9 in rng:
                                            mats.append((n1,n2,n3,n4,n5,n6,n7,n8,n9))
                                        #end for
                                    #end for
                                #end for
                            #end for
                        #end for
                    #end for
                #end for
            #end for
        #end for
        mats = array(mats,dtype=int)
        mats.shape = (2*dn+1)**(dim*dim),dim,dim
        opt_tm_matrices[dn] = mats
    else:
        mats = opt_tm_matrices[dn]
    #end if
    # calculate and store all wigner image indices
    if nc not in opt_tm_wig_indices:
        inds = []
        rng = tuple(range(-nc,nc+1))
        for k in rng:
            for j in rng:
                for i in rng:
                    if i!=0 or j!=0 or k!=0:
                        inds.append((i,j,k))
                    #end if
                #end for
            #end for
        #end for
        inds = array(inds,dtype=int)
        opt_tm_wig_indices[nc] = inds
    else:
        inds = opt_tm_wig_indices[nc]
    #end if
    # track counts of tiling matrices
    ntilings        = len(mats)
    nequiv_volume   = 0
    nfilter         = 0
    nequiv_inscribe = 0
    nequiv_wigner   = 0
    nequiv_cubicity = 0
    nequiv_shape    = 0
    # try a faster search for cells w/ target volume
    det_inds_p = [
        [(0,0),(1,1),(2,2)],
        [(0,1),(1,2),(2,0)],
        [(0,2),(1,0),(2,1)]
        ]
    det_inds_m = [
        [(0,0),(1,2),(2,1)],
        [(0,1),(1,0),(2,2)],
        [(0,2),(1,1),(2,0)]
        ]
    volfacs = zeros((len(mats),),dtype=int)
    for (i1,j1),(i2,j2),(i3,j3) in det_inds_p:
        volfacs += (Tref[i1,j1]+mats[:,i1,j1])*(Tref[i2,j2]+mats[:,i2,j2])*(Tref[i3,j3]+mats[:,i3,j3])
    #end for
    for (i1,j1),(i2,j2),(i3,j3) in det_inds_m:
        volfacs -= (Tref[i1,j1]+mats[:,i1,j1])*(Tref[i2,j2]+mats[:,i2,j2])*(Tref[i3,j3]+mats[:,i3,j3])
    #end for
    Tmats = mats[abs(volfacs)==volfac]
    nequiv_volume = len(Tmats)    
    # find the set of cells with maximal inscribing radius
    inscribe_tilings = []
    rmax = -1e99
    for mat in Tmats:
        T = Tref + mat
        if filter(T):
            nfilter+=1
            Taxes = dot(T,axes)
            rc1 = norm(cross(Taxes[0],Taxes[1]))
            rc2 = norm(cross(Taxes[1],Taxes[2]))
            rc3 = norm(cross(Taxes[2],Taxes[0]))
            r   = 0.5*volume/max(rc1,rc2,rc3) # inscribing radius
            if r>rmax or abs(r-rmax)<tol:
                inscribe_tilings.append((r,T,Taxes))
                rmax = r
            #end if
        #end if
    #end for
    # find the set of cells w/ maximal wigner radius out of the inscribing set
    wigner_tilings = []
    rwmax = -1e99
    for r,T,Taxes in inscribe_tilings:
        if abs(r-rmax)<tol:
            nequiv_inscribe+=1
            rw = 1e99
            for ind in inds:
                rw = min(rw,0.5*norm(dot(ind,Taxes)))
            #end for
            if rw>rwmax or abs(rw-rwmax)<tol:
                wigner_tilings.append((rw,T,Taxes))
                rwmax = rw
            #end if
        #end if
    #end for
    # find the set of cells w/ maximal cubicity
    # (minimum cube_deviation)
    cube_tilings = []            
    cmin = 1e99
    for rw,T,Ta in wigner_tilings:
        if abs(rw-rwmax)<tol:
            nequiv_wigner+=1
            dc = volume**(1./3)*sqrt(2.)
            d1 = abs(norm(Ta[0]+Ta[1])-dc)
            d2 = abs(norm(Ta[1]+Ta[2])-dc)
            d3 = abs(norm(Ta[2]+Ta[0])-dc)
            d4 = abs(norm(Ta[0]-Ta[1])-dc)
            d5 = abs(norm(Ta[1]-Ta[2])-dc)
            d6 = abs(norm(Ta[2]-Ta[0])-dc)
            cube_dev = (d1+d2+d3+d4+d5+d6)/(6*dc)
            if cube_dev<cmin or abs(cube_dev-cmin)<tol:
                cube_tilings.append((cube_dev,rw,T,Ta))
                cmin = cube_dev
            #end if
        #end if
    #end for
    # prioritize selection by "shapeliness" of tiling matrix
    #   prioritize positive diagonal elements
    #   penalize off-diagonal elements
    #   penalize negative off-diagonal elements
    shapely_tilings = []
    smax = -1e99
    for cd,rw,T,Taxes in cube_tilings:
        if abs(cd-cmin)<tol:
            nequiv_cubicity+=1
            d = diag(T)
            o = (T-diag(d)).ravel()
            s = sign(d).sum()-(abs(o)>0).sum()-(o<0).sum()
            if s>smax or abs(s-smax)<tol:
                shapely_tilings.append((s,rw,T,Taxes))
                smax = s
            #end if
        #end if
    #end for
    # prioritize selection by symmetry of tiling matrix
    ropt   = -1e99
    Topt   = None
    Taxopt = None
    diagonal      = []
    symmetric     = []
    antisymmetric = []
    other         = []
    for s,rw,T,Taxes in shapely_tilings:
        if abs(s-smax)<tol:
            nequiv_shape+=1
            Td = diag(diag(T))
            if abs(Td-T).sum()==0:
                diagonal.append((rw,T,Taxes))
            elif abs(T.T-T).sum()==0:
                symmetric.append((rw,T,Taxes))
            elif abs(T.T+T-2*Td).sum()==0:
                antisymmetric.append((rw,T,Taxes))
            else:
                other.append((rw,T,Taxes))
            #end if
        #end if
    #end for
    s = 1
    if len(diagonal)>0:
        cells = diagonal
    elif len(symmetric)>0:
        cells = symmetric
    elif len(antisymmetric)>0:
        cells = antisymmetric
        s = -1
    elif len(other)>0:
        cells = other
    #end if
    skew_min = 1e99
    if len(cells)>0:
        for rw,T,Taxes in cells:
            Td = diag(diag(T))
            skew = abs(T.T-s*T-(1-s)*Td).sum()
            if skew<skew_min:
                ropt = rw
                Topt = T
                Taxopt = Taxes
                skew_min = skew
            #end if
        #end for
    #end if
    if Taxopt is None:
        error('optimal tilematrix for volfac={0} not found with tolerance {1}\ndifference range (dn): {2}\ntiling matrices searched: {3}\ncells with target volume: {4}\ncells that passed the filter: {5}\ncells with equivalent inscribing radius: {6}\ncells with equivalent wigner radius: {7}\ncells with equivalent cubicity: {8}\nmatrices with equivalent shapeliness: {9}\nplease try again with dn={10}'.format(volfac,tol,dn,ntilings,nequiv_volume,nfilter,nequiv_inscribe,nequiv_wigner,nequiv_cubicity,nequiv_shape,dn+1))
    #end if
    if det(Taxopt)<0:
        Topt = -Topt
    #end if
    return Topt,ropt
#end def optimal_tilematrix




class Sobj(DevBase):
    None
#end class Sobj



class Structure(Sobj): 

    operations = obj()

    @classmethod
    def set_operations(cls):
        cls.operations.set(
            remove_folded_structure = cls.remove_folded_structure,
            recenter = cls.recenter,
            )
    #end def set_operations


    def __init__(self,axes=None,scale=1.,elem=None,pos=None,mag=None,
                 center=None,kpoints=None,kweights=None,kgrid=None,kshift=None,
                 permute=None,units=None,tiling=None,rescale=True,dim=3,
                 magnetization=None,magnetic_order=None,magnetic_prim=True,
                 operations=None,background_charge=0,frozen=None,bconds=None,
                 posu=None):

        if center is None:
            if axes is not None:
                center = array(axes,dtype=float).sum(0)/2
            else:
                center = dim*[0]
            #end if
        #end if
        if bconds is None or bconds=='periodic':
            bconds = dim*['p']
        #end if
        if axes is None:
            axes   = []
            bconds = []
        #end if
        if elem is None:
            elem = []
        #end if
        if posu!=None:
            pos = posu
        #end if
        if pos is None:
            pos = empty((0,dim))
        #end if
        if kshift is None:
            kshift = 0,0,0
        #end if
        if mag is None:
            mag = len(elem)*[None]
        #end if
        self.scale    = 1.
        self.units    = units
        self.dim      = dim
        self.center   = array(center,dtype=float)
        self.axes     = array(axes,dtype=float)
        self.set_bconds(bconds)
        self.set_elem(elem)
        self.pos      = array(pos,dtype=float)
        self.frozen   = None
        self.mag      = array(mag,dtype=object)
        self.kpoints  = empty((0,dim))            
        self.kweights = empty((0,))         
        self.background_charge = background_charge
        self.remove_folded_structure()
        if len(axes)==0:
            self.kaxes=array([])
        else:
            self.kaxes=2*pi*inv(self.axes).T
        #end if
        if posu!=None:
            self.pos_to_cartesian()
        #end if
        if frozen!=None:
            self.frozen = array(frozen,dtype=bool)
            if self.frozen.shape!=self.pos.shape:
                self.error('frozen directions must have the same shape as positions\n  positions shape: {0}\n  frozen directions shape: {1}'.format(self.pos.shape,self.frozen.shape))
            #end if
        #end if
        self.magnetize(magnetization)
        if tiling is not None:
            self.tile(tiling,
                      in_place = True,
                      magnetic_order = magnetic_order,
                      magnetic_prim  = magnetic_prim
                      )
        #end if
        if kpoints!=None:
            self.add_kpoints(kpoints,kweights)
        #end if
        if kgrid!=None:
            self.add_kmesh(kgrid,kshift)
        #end if        
        if rescale:
            self.rescale(scale)
        else:
            self.scale = scale
        #end if
        if permute!=None:
            self.permute(permute)
        #end if
        if operations!=None:
            self.operate(operations)
        #end if
    #end def __init__


    def check_consistent(self,tol=1e-8,exit=True,message=False):
        msg = ''
        if self.has_axes():
            kaxes = 2*pi*inv(self.axes).T
            abs_diff = abs(self.kaxes-kaxes).sum()
            if abs_diff>tol:
                msg += 'direct and reciprocal space axes are not consistent\naxes present:\n{0}\nkaxes present:\n{1}\nconsistent kaxes:\n{2}\nabsolute difference: {3}\n'.format(self.axes,self.kaxes,kaxes,abs_diff)
            #end if
        #end if
        consistent = len(msg)==0
        if not consistent and exit:
            self.error(msg)
        #end if
        if not message:
            return consistent
        else:
            return consistent,msg
        #end if
    #end def check_consistent


    def set_axes(self,axes):
        self.reset_axes(axes)
    #end def set_axes


    def set_bconds(self,bconds):
        self.bconds = array(tuple(bconds),dtype=str)
    #end def bconds


    def set_elem(self,elem):
        self.elem = array(elem,dtype=object)
    #end def set_elem

    
    def set_pos(self,pos):
        self.pos = array(pos,dtype=float)
    #end def set_pos


    def size(self):
        return len(self.elem)
    #end def size

    
    def has_axes(self):
        return len(self.axes)==self.dim
    #end def has_axes


    def operate(self,operations):
        for op in operations:
            if not op in self.operations:
                self.error('{0} is not a known operation\n  valid options are:\n    {1}'.format(op,list(self.operations.keys())))
            else:
                self.operations[op](self)
            #end if
        #end for
    #end def operate


    def set_folded(self,folded):
        self.set_folded_structure(folded)
    #end def set_folded


    def remove_folded(self):
        self.remove_folded_structure()
    #end def remove_folded

    
    def has_folded(self):
        return self.has_folded_structure()
    #end def has_folded


    def set_folded_structure(self,folded):
        self.folded_structure = folded
        self.tmatrix = self.tilematrix(folded)
    #end def set_folded_structure


    def remove_folded_structure(self):
        self.folded_structure = None
        self.tmatrix = None
    #end def remove_folded_structure

        
    def has_folded_structure(self):
        return self.folded_structure!=None
    #end def has_folded_structure

            
    def group_atoms(self,folded=True):
        if len(self.elem)>0:
            order = self.elem.argsort()
            if (self.elem!=self.elem[order]).any():
                self.elem = self.elem[order]
                self.pos  = self.pos[order]
            #end if
        #end if
        if self.folded_structure!=None and folded:
            self.folded_structure.group_atoms(folded)
        #end if
    #end def group_atoms


    def rename(self,folded=True,**name_pairs):
        elem = self.elem
        for old,new in name_pairs.iteritems():
            for i in xrange(len(self.elem)):
                if old==elem[i]:
                    elem[i] = new
                #end if
            #end for
        #end for
        if self.folded_structure!=None and folded:
            self.folded_structure.rename(folded=folded,**name_pairs)
        #end if
    #end def rename


    def reset_axes(self,axes=None):
        if axes is None:
            axes = self.axes
        else:
            axes = array(axes)
            self.remove_folded_structure()
        #end if
        self.axes  = axes
        self.kaxes = 2*pi*inv(axes).T
        self.center = axes.sum(0)/2
    #end def reset_axes


    def adjust_axes(self,axes):
        self.skew(dot(inv(self.axes),axes))
    #end def adjust_axes
        

    def reshape_axes(self,reshaping):
        R = array(reshaping)
        if abs(abs(det(R))-1)<1e-6:
            self.axes = dot(self.axes,R)
        else:
            R = dot(inv(self.axes),R)
            if abs(abs(det(R))-1)<1e-6:
                self.axes = dot(self.axes,R)
            else:
                self.error('reshaping matrix must not change the volume\n  reshaping matrix:\n  {0}\n  volume change ratio: {1}'.format(R,abs(det(R))))
            #end if
        #end if
    #end def reshape_axes

    
    def corners(self):
        a = self.axes
        c = array([(0,0,0),
                   a[0],
                   a[1],
                   a[2],
                   a[0]+a[1],
                   a[1]+a[2],
                   a[2]+a[0],
                   a[0]+a[1]+a[2],
                   ])
        return c
    #end def corners

    
    def miller_direction(self,h,k,l,normalize=False):
        d = dot((h,k,l),self.axes)
        if normalize:
            d/=norm(d)
        #end if
        return d
    #end def miller_direction

    
    def miller_normal(self,h,k,l,normalize=False):
        d = dot((h,k,l),self.kaxes)
        if normalize:
            d/=norm(d)
        #end if
        return d
    #end def miller_normal


    def project_plane(self,a1,a2,points=None):
        # a1/a2: in plane vectors
        if points is None:
            points = self.pos
        #end if
        a1n = norm(a1)
        a2n = norm(a2)
        a1/=a1n
        a2/=a2n
        n = cross(a1,a2)
        plane_coords = []
        for p in points:
            p -= dot(n,p)*n # project point into plane
            c1 = dot(a1,p)/a1n
            c2 = dot(a2,p)/a2n
            plane_coords.append((c1,c2))
        #end for
        return array(plane_coords,dtype=float)
    #end def project_plane

        
    def bounding_box(self,scale=1.0,box='tight',recenter=False):
        pmin    = self.pos.min(0)
        pmax    = self.pos.max(0)
        pcenter = (pmax+pmin)/2
        prange  = pmax-pmin
        if box=='tight':
            axes = diag(prange)
        elif box=='cubic' or box=='cube':
            prmax = prange.max()
            axes = diag((prmax,prmax,prmax))
        elif isinstance(box,ndarray) or isinstance(box,list):
            box = array(box)
            if box.shape!=(3,3):
                self.error('requested box must be 3-dimensional (3x3 axes)\n  you provided: '+str(box)+'\n shape: '+str(box.shape))
            #end if
            binv = inv(box)
            pu = dot(self.pos,binv)
            pmin    = pu.min(0)
            pmax    = pu.max(0)
            pcenter = (pmax+pmin)/2
            prange  = pmax-pmin
            axes    = dot(diag(prange),box)
        else:
            self.error("invalid request for box\n  valid options are 'tight', 'cubic', or axes array (3x3)\n  you provided: "+str(box))
        #end if
        self.reset_axes(scale*axes)
        self.slide(self.center-pcenter,recenter)
    #end def bounding_box


    def center_molecule(self):
        self.slide(self.center-self.pos.mean(0),recenter=False)
    #end def center_molecule


    def center_solid(self):
        u = self.pos_unit()
        du = (1-u.min(0)-u.max(0))/2
        self.slide(dot(du,self.axes),recenter=False)
    #end def center_solid


    def permute(self,permutation):
        dim = self.dim
        P = empty((dim,dim),dtype=int)
        if len(permutation)!=dim:
            self.error(' permutation vector must have {0} elements\n you provided {1}'.format(dim,permutation))
        #end if
        for i in range(dim):
            p = permutation[i]
            pv = zeros((dim,),dtype=int)
            if p=='x' or p=='0':
                pv[0] = 1
            elif p=='y' or p=='1':
                pv[1] = 1
            elif p=='z' or p=='2':
                pv[2] = 1
            #end if
            P[:,i] = pv[:]
        #end for
        self.center = dot(self.center,P)
        if self.has_axes():
            self.axes = dot(self.axes,P)
        #end if
        if len(self.pos)>0:
            self.pos = dot(self.pos,P)
        #end if
        if len(self.kaxes)>0:
            self.kaxes = dot(self.kaxes,P)
        #end if
        if len(self.kpoints)>0:
            self.kpoints = dot(self.kpoints,P)
        #end if
        if self.folded_structure!=None:
            self.folded_structure.permute(permutation)
        #end if
    #end def permute


    def rotate_plane(self,plane,angle,units='degrees'):
        self.pos = rotate_plane(plane,angle,self.pos,units)
        if self.has_axes():
            axes = rotate_plane(plane,angle,self.axes,units)
            self.reset_axes(axes)
        #end if
    #end def rotate_plane


    def upcast(self,DerivedStructure):
        if not issubclass(DerivedStructure,Structure):
            self.error(DerivedStructure.__name__,'is not derived from Structure')
        #end if
        ds = DerivedStructure()
        for name,value in self.iteritems():
            ds[name] = deepcopy(value)
        #end for
        return ds
    #end def upcast

    
    def incorporate(self,other):
        self.set_elem(list(self.elem)+list(other.elem))
        self.pos=array(list(self.pos)+list(other.pos))
    #end def incorporate


    def add_atoms(self,elem,pos):
        self.set_elem(list(self.elem)+list(elem))
        self.pos=array(list(self.pos)+list(pos))
    #end def add_atoms


    def is_open(self):
        return not self.any_periodic()
    #end def is_open


    def is_periodic(self):
        return self.any_periodic()
    #end def is_periodic


    def any_periodic(self):
        has_cell    = self.has_axes()
        pbc = False
        for bc in self.bconds:
            pbc |= bc=='p'
        #end if
        periodic = has_cell and pbc
        return periodic
    #end def any_periodic

    
    def all_periodic(self):
        has_cell = self.has_axes()
        pbc = True
        for bc in self.bconds:
            pbc &= bc=='p'
        #end if
        periodic = has_cell and pbc
        return periodic
    #end def all_periodic


    def distances(self,pos1=None,pos2=None):
        if isinstance(pos1,Structure):
            pos1 = pos1.pos
        #end if
        if pos2==None:
            if pos1==None:
                return sqrt((self.pos**2).sum(1))
            else:
                pos2 = self.pos
            #end if
        #end if
        if len(pos1)!=len(pos2):
            self.error('positions arrays are not the same length')
        #end if
        return sqrt(((pos1-pos2)**2).sum(1))
    #end def distances

    
    def volume(self):
        if not self.has_axes():
            return None
        else:
            return abs(det(self.axes))
        #end if
    #end def volume


    def rwigner(self,nc=5):
        if self.dim!=3:
            self.error('rwigner is currently only implemented for 3 dimensions')
        #end if
        rmin = 1e90
        n=empty((1,3))
        rng = tuple(range(-nc,nc+1))
        for k in rng:
            for j in rng:
                for i in rng:
                    if i!=0 or j!=0 or k!=0:
                        n[:] = i,j,k
                        rmin = min(rmin,.5*norm(dot(n,self.axes)))
                    #end if
                #end for
            #end for
        #end for
        return rmin
    #end def rwigner


    def rinscribe(self):
        if self.dim!=3:
            self.error('rinscribe is currently only implemented for 3 dimensions')
        #end if
        radius = 1e99
        dim=3
        axes=self.axes
        for d in range(dim):
            i = d
            j = (d+1)%dim
            rc = cross(axes[i,:],axes[j,:])
            radius = min(radius,.5*abs(det(axes))/norm(rc))
        #end for
        return radius
    #end def rinscribe


    def rwigner_cube(self,*args,**kwargs):
        cube = Structure()
        a = self.volume()**(1./3)
        cube.set_axes([[a,0,0],[0,a,0],[0,0,a]])
        return cube.rwigner(*args,**kwargs)
    #end def rwigner_cube


    def rinscribe_cube(self,*args,**kwargs):
        cube = Structure()
        a = self.volume()**(1./3)
        cube.set_axes([[a,0,0],[0,a,0],[0,0,a]])
        return cube.rinscribe(*args,**kwargs)
    #end def rinscribe_cube


    def rmin(self):
        return self.rwigner()
    #end def rmin


    def rcell(self):
        return self.rinscribe()
    #end def rcell

    # scale invariant measure of deviation from cube shape
    #   based on deviation of face diagonals from cube
    def cube_deviation(self):
        a = self.axes
        dc = self.volume()**(1./3)*sqrt(2.)
        d1 = abs(norm(a[0]+a[1])-dc)
        d2 = abs(norm(a[1]+a[2])-dc)
        d3 = abs(norm(a[2]+a[0])-dc)
        d4 = abs(norm(a[0]-a[1])-dc)
        d5 = abs(norm(a[1]-a[2])-dc)
        d6 = abs(norm(a[2]-a[0])-dc)
        return (d1+d2+d3+d4+d5+d6)/(6*dc)
    #end def cube_deviation
    
    # apply volume preserving shear-removing transformations to cell axes
    #   resulting unsheared cell has orthogonal axes
    #    while remaining periodically correct
    #   note that the unshearing procedure is not unique
    #   it depends on the order of unshearing operations
    def unsheared_axes(self,axes=None,distances=False):
        if self.dim!=3:
            self.error('unsheared_axes is currently only implemented for 3 dimensions')
        #end if
        if axes is None:
            axes = self.axes
        #end if
        dim=3
        axbar = identity(dim)
        axnew = array(axes,dtype=float)
        dists = empty((dim,))
        for d in range(dim):
            d2 = (d+1)%dim
            d3 = (d+2)%dim
            n = cross(axnew[d2],axnew[d3])  #vector normal to 2 cell faces
            axdist = dot(n,axes[d])/dot(n,axbar[d])
            axnew[d]  = axdist*axbar[d]
            dists[d] = axdist
        #end for
        if not distances:
            return axnew
        else:
            return axnew,dists
        #end if
    #end def unsheared_axes


    # vectors parallel to cell faces
    #   length of vectors is distance between parallel face planes
    #   note that the product of distances is not the cell volume in general
    #   see "unsheared_axes" function
    #   (e.g. a volume preserving shear may bring two face planes arbitrarily close)
    def face_vectors(self,axes=None,distances=False):
        if axes is None:
            axes = self.axes
        #end if
        fv = inv(axes).T
        for d in range(len(fv)): 
            fv[d] /= norm(fv[d]) # face normals
        #end for
        dv = dot(axes,fv.T) # axis projections onto face normals
        fv = dot(dv,fv)     # face normals lengthened by plane separation
        if not distances:
            return fv
        else:
            return fv,diag(dv)
        #end if
    #end def face_vectors


    def face_distances(self):
        return self.face_vectors(distances=True)[1]
    #end def face_distances

    
    def rescale(self,scale):
        self.scale  *= scale
        self.axes   *= scale
        self.pos    *= scale
        self.center *= scale
        self.kaxes  /= scale
        self.kpoints/= scale
        if self.folded_structure!=None:
            self.folded_structure.rescale(scale)
        #end if
    #end def rescale


    def stretch(self,s1,s2,s3):
        if self.dim!=3:
            self.error('stretch is currently only implemented for 3 dimensions')
        #end if
        d = diag((s1,s2,s3))
        self.skew(d)
    #end def stretch

        
    def skew(self,skew):
        axinv  = inv(self.axes)
        axnew  = dot(self.axes,skew)
        kaxinv = inv(self.kaxes)
        kaxnew = dot(inv(skew),self.kaxes)
        self.pos     = dot(dot(self.pos,axinv),axnew)
        self.center  = dot(dot(self.center,axinv),axnew)
        self.kpoints = dot(dot(self.kpoints,kaxinv),kaxnew)
        self.axes  = axnew
        self.kaxes = kaxnew
        if self.folded_structure!=None:
            self.folded_structure.skew(skew)
        #end if
    #end def skew
        
    
    def change_units(self,units,folded=True):
        if units!=self.units:
            scale = convert(1,self.units,units)
            self.scale  *= scale
            self.axes   *= scale
            self.pos    *= scale
            self.center *= scale
            self.kaxes  /= scale
            self.kpoints/= scale
            self.units  = units
        #end if
        if self.folded_structure!=None and folded:
            self.folded_structure.change_units(units,folded=folded)
        #end if
    #end def change_units
                              
        
    # insert sep space at loc along axis
    #   if sep<0, space is removed instead
    def cleave(self,axis,loc,sep=None,remove=False,tol=1e-6):
        self.remove_folded_structure()
        if isinstance(axis,int):
            if sep is None:
                self.error('separation induced by cleave must be provided')
            #end if
            v = self.face_vectors()[axis]
            if isinstance(loc,float):
                c = loc*v/norm(v)
            #end if
        else:
            v = axis
        #end if
        c = array(c)  # point on cleave plane
        v = array(v)  # normal vector to cleave plane, norm is cleave separation
        if sep!=None:
            v = abs(sep)*v/norm(v)
        #end if
        if norm(v)<tol:
            return
        #end if
        vn = array(v/norm(v))
        if sep!=None and sep<0:
            v = -v # preserve the normal direction for atom identification, but reverse the shift direction
        #end if
        self.recorner()  # want box contents to be static
        if self.has_axes():
            components = 0
            dim = self.dim
            axes = self.axes
            for i in xrange(dim):
                i2 = (i+1)%dim
                i3 = (i+2)%dim
                a2 = axes[i2]/norm(axes[i2])
                a3 = axes[i3]/norm(axes[i3])
                comp = abs(dot(a2,vn))+abs(dot(a3,vn))
                if comp < 1e-6:
                    components+=1
                    iaxis = i
                #end if
            #end for
            commensurate = components==1
            if not commensurate:
                self.error('cannot insert vacuum because cleave is incommensurate with the cell\n  cleave plane must be parallel to a cell face')
            #end if
            a = self.axes[iaxis]
            #self.axes[iaxis] = (1.+dot(v,a)/dot(a,a))*a
            self.axes[iaxis] = (1.+dot(v,v)/dot(v,a))*a
        #end if
        indices = []
        pos = self.pos
        for i in xrange(len(pos)):
            p = pos[i]
            comp = dot(p-c,vn)
            if comp>0 or abs(comp)<tol:
                pos[i] += v
                indices.append(i)
            #end if
        #end for
        if remove:
            self.remove(indices)
        #end if
    #end def cleave


    def translate(self,v):
        v = array(v)
        pos = self.pos
        for i in range(len(pos)):
            pos[i]+=v
        #end for
        self.center+=v
        if self.folded_structure!=None:
            self.folded_structure.translate(v)
        #end if
    #end def translate

                              
    def slide(self,v,recenter=True):
        v = array(v)
        pos = self.pos
        for i in range(len(pos)):
            pos[i]+=v
        #end for
        if recenter:
            self.recenter()
        #end if
        if self.folded_structure!=None:
            self.folded_structure.slide(v,recenter)
        #end if
    #end def slide


    def zero_corner(self):
        corner = self.center-self.axes.sum(0)/2
        self.translate(-corner)
    #end def zero_corner


    def locate_simple(self,pos):
        pos = array(pos)
        if pos.shape==(self.dim,):
            pos = [pos]
        #end if
        nn = nearest_neighbors(1,self.pos,pos)
        return nn.ravel()
    #end def locate_simple

    
    def locate(self,identifiers,radii=None,exterior=False):
        indices = None
        if isinstance(identifiers,Structure):
            cell = identifiers
            indices = cell.inside(self.pos)
        elif isinstance(identifiers,ndarray) and identifiers.dtype==bool:
            indices = arange(len(self.pos))[identifiers]
        elif isinstance(identifiers,int):
            indices = [identifiers]
        elif len(identifiers)>0 and isinstance(identifiers[0],int):
            indices = identifiers
        elif isinstance(identifiers,str):
            atom = identifiers
            indices = []
            for i in xrange(len(self.elem)):
                if self.elem[i]==atom:
                    indices.append(i)
                #end if
            #end for
        elif len(identifiers)>0 and isinstance(identifiers[0],str):
            indices = []
            for atom in identifiers:
                for i in xrange(len(self.elem)):
                    if self.elem[i]==atom:
                        indices.append(i)
                    #end if
                #end for
            #end for
        #end if
        if radii!=None or indices==None:
            if indices is None:
                pos = identifiers
            else:
                pos = self.pos[indices]
            #end if
            if isinstance(radii,float) or isinstance(radii,int):
                radii = len(pos)*[radii]
            elif radii!=None and len(radii)!=len(pos):
                self.error('lengths of input radii and positions do not match\n  len(radii)={0}\n  len(pos)={1}'.format(len(radii),len(pos)))
            #end if
            dtable = self.min_image_distances(pos)
            indices = []
            if radii is None:
                for i in range(len(pos)):
                    indices.append(dtable[i].argmin())
                #end for
            else:
                ipos = arange(len(self.pos))
                for i in range(len(pos)):
                    indices.extend(ipos[dtable[i]<radii[i]])
                #end for
            #end if
        #end if
        if exterior:
            indices = list(set(range(len(self.pos)))-set(indices))
        #end if
        return indices
    #end def locate

    
    def freeze(self,identifiers=None,radii=None,exterior=False,negate=False,directions='xyz'):
        if isinstance(identifiers,ndarray) and identifiers.shape==self.pos.shape and identifiers.dtype==bool:
            if negate:
                self.frozen = ~identifiers
            else:
                self.frozen = identifiers.copy()
            #end if
            return
        #end if
        if identifiers is None:
            indicies = arange(len(self.pos),dtype=int)
        else:
            indices = self.locate(identifiers,radii,exterior)
        #end if
        if len(indices)==0:
            self.error('failed to select any atoms to freeze')
        #end if
        if isinstance(directions,str):
            d = empty((3,),dtype=bool)
            d[0] = 'x' in directions
            d[1] = 'y' in directions
            d[2] = 'z' in directions
            directions = len(indices)*[d]
        else:
            directions = array(directions,dtype=bool)
        #end if
        if self.frozen is None:
            self.frozen = zeros(self.pos.shape,dtype=bool)
        #end if
        frozen = self.frozen
        i=0
        if not negate:
            for index in indices:
                frozen[index] = directions[i]
                i+=1
            #end for
        else:
            for index in indices:
                frozen[index] = directions[i]==False
                i+=1
            #end for
        #end if
    #end def freeze


    def magnetize(self,identifiers=None,magnetization='',**mags):
        magsin = None
        if isinstance(identifiers,obj):
            magsin = identifiers.copy()
        elif isinstance(magnetization,obj):
            magsin = magnetization.copy()
        #endif
        if magsin!=None:
            magsin.transfer_from(mags)
            mags = magsin
            identifiers = None
            magnetization = ''
        #end if
        for e,m in mags.iteritems():
            if not e in self.elem:
                self.error('cannot magnetize non-existent element {0}'.format(e))
            elif not isinstance(m,(NoneType,int)):
                self.error('magnetizations provided must be either None or integer\n  you provided: {0}\n  full magnetization request provided:\n {1}'.format(m,mags))
            #end if
            self.mag[self.elem==e] = m
        #end for
        if identifiers is None and magnetization=='':
            return
        elif magnetization=='':
            magnetization = identifiers
            indices = range(len(self.elem))
        else:
            indices = self.locate(identifiers)
        #end if
        if not isinstance(magnetization,(list,tuple,ndarray)):
            magnetization = [magnetization]
        #end if
        for m in magnetization:
            if not isinstance(m,(NoneType,int)):
                self.error('magnetizations provided must be either None or integer\n  you provided: {0}\n  full magnetization list provided: {1}'.format(m,magnetization))
            #end if
        #end for
        if len(magnetization)==1:
            m = magnetization[0]
            for i in indices:
                self.mag[i] = m
            #end for
        elif len(magnetization)==len(indices):
            for i in range(len(indices)):
                self.mag[indices[i]] = magnetization[i]
            #end for
        else:
            self.error('magnetization list and list selected atoms differ in length\n  length of magnetization list: {0}\n  number of atoms selected: {1}\n  magnetization list: {2}\n  atom indices selected: {3}\n  atoms selected: {4}'.format(len(magnetization),len(indices),magnetization,indices,self.elem[indices]))
        #end if
    #end def magnetize


    def carve(self,identifiers):
        indices = self.locate(identifiers)
        if isinstance(identifiers,Structure):
            sub = identifiers
            sub.elem = self.elem[indices].copy()
            sub.pos  = self.pos[indices].copy()
        else:
            sub = self.copy()
            sub.elem = self.elem[indices]
            sub.pos  = self.pos[indices]
        #end if
        sub.host_indices = array(indices)
        return sub
    #end def carve

        
    def remove(self,identifiers):
        indices = self.locate(identifiers)
        keep = list(set(range(len(self.pos)))-set(indices))
        erem = self.elem[indices]
        prem = self.pos[indices]
        self.elem = self.elem[keep]
        self.pos  = self.pos[keep]
        self.remove_folded_structure()
        return erem,prem
    #end def remove

    
    def replace(self,identifiers,elem=None,pos=None,radii=None,exterior=False):
        indices = self.locate(identifiers,radii,exterior)
        if isinstance(elem,Structure):
            cell = elem
            elem = cell.elem
            pos  = cell.pos
        elif elem==None:
            elem = self.elem
        #end if
        indices=array(indices)
        elem=array(elem,dtype=object)
        pos =array(pos)
        nrem = len(indices)
        nadd = len(pos)
        if nadd<nrem:
            ar = array(range(0,nadd))
            rr = array(range(nadd,nrem))
            self.elem[indices[ar]] = elem[:]
            self.pos[indices[ar]]  = pos[:]
            self.remove(indices[rr])
        elif nadd>nrem:
            ar = array(range(0,nrem))
            er = array(range(nrem,nadd))
            self.elem[indices[ar]] = elem[ar]
            self.pos[indices[ar]]  = pos[ar]
            ii = indices[ar[-1]]
            self.set_elem( list(self.elem[0:ii])+list(elem[er])+list(self.elem[ii:]) )
            self.pos = array( list(self.pos[0:ii])+list(pos[er])+list(self.pos[ii:]) )
        else:
            self.elem[indices] = elem[:]
            self.pos[indices]  = pos[:]
        #end if
        self.remove_folded_structure()
    #end def replace


    def replace_nearest(self,elem,pos=None):
        if isinstance(elem,Structure):
            cell = elem
            elem = cell.elem
            pos  = cell.pos
        #end if
        nn = nearest_neighbors(1,self.pos,pos)
        np = len(pos)
        nps= len(self.pos)
        d = empty((np,))
        ip = array(range(np))
        ips= nn.ravel()
        for i in ip:
            j = ips[i]
            d[i]=sqrt(((pos[i]-self.pos[j])**2).sum())
        #end for
        order = d.argsort()
        ip = ip[order]
        ips=ips[order]
        replacable = empty((nps,))
        replacable[:] = False
        replacable[ips]=True
        insert = []
        last_replaced=nps-1
        for n in range(np):
            i = ip[n]
            j = ips[n]
            if replacable[j]:
                self.pos[j] = pos[i]
                self.elem[j]=elem[i]
                replacable[j]=False
                last_replaced = j
            else:
                insert.append(i)
            #end if
        #end for
        insert=array(insert)
        ii = last_replaced
        if len(insert)>0:
            self.set_elem( list(self.elem[0:ii])+list(elem[insert])+list(self.elem[ii:]) )
            self.pos = array( list(self.pos[0:ii])+list(pos[insert])+list(self.pos[ii:]) )
        #end if
        self.remove_folded_structure()
    #end def replace_nearest


    def point_defect(self,identifiers=None,elem=None,dr=None):
        if isinstance(elem,str):
            elem = [elem]
            if dr!=None:
                dr = [dr]
            #end if
        #end if
        if not 'point_defects' in self:
            self.point_defects = obj()
        #end if
        point_defects = self.point_defects
        ncenters = len(point_defects)
        if identifiers is None:
            index = ncenters
            if index>=len(self.pos):
                self.error('attempted to add a point defect at index {0}, which does not exist\n  for reference there are {1} atoms in the structure'.format(index,len(self.pos)))
            #end if
        else:
            indices = self.locate(identifiers)
            if len(indices)>1:
                self.error('{0} atoms were located by identifiers provided\n  a point defect replaces only a single atom\n  atom indices located: {1}'.format(len(indices),indices))
            #end if
            index = indices[0]
        #end if
        if elem is None:
            self.error('must supply substitutional elements comprising the point defect\n  expected a list or similar for input argument elem')
        elif len(elem)>1 and dr is None:
            self.error('must supply displacements (dr) since many atoms comprise the point defect')
        elif dr!=None and len(elem)!=len(dr):
            self.error('elem and dr must have the same length')
        #end if
        r = self.pos[index]
        e = self.elem[index]
        elem = array(elem)
        pos = zeros((len(elem),len(r)))
        if dr is None:
            rc = r
            for i in range(len(elem)):
                pos[i] = r
            #end for
        else:
            nrc = 0
            rc  = 0*r
            dr = array(dr)
            for i in range(len(elem)):
                pos[i] = r + dr[i]
                if norm(dr[i])>1e-5:
                    rc+=dr[i]
                    nrc+=1
                #end if
            #end for
            if nrc==0:
                rc = r
            else:
                rc = r + rc/nrc
            #end if
        #end if
        point_defect = obj(
            center = rc,
            elem_replaced = e,
            elem = elem,
            pos  = pos
            )
        point_defects.append(point_defect)
        elist = list(self.elem)
        plist = list(self.pos)
        if len(elem)==0 or len(elem)==1 and elem[0]=='':
            elist.pop(index)
            plist.pop(index)
        else:
            elist[index] = elem[0]
            plist[index] = pos[0]
            for i in range(1,len(elem)):
                elist.append(elem[i])
                plist.append(pos[i])
            #end for
        #end if
        self.set_elem(elist)
        self.pos  = array(plist)
        self.remove_folded_structure()
    #end def point_defect


    def species(self,symbol=False):
        if not symbol:
            return set(self.elem)
        else:
            species_labels = set(self.elem)
            species = set()
            for e in species_labels:
                is_elem,symbol = is_element(e,symbol=True)
                species.add(symbol)
            #end for
            return species_labels,species
        #end if
    #end def species

        
    def ordered_species(self,symbol=False):
        speclab_set    = set()
        species_labels = []
        if not symbol:
            for e in self.elem:
                if e not in speclab_set:
                    speclab_set.add(e)
                    species_labels.append(e)
                #end if
            #end for
            return species_labels
        else:
            species  = []
            spec_set = set()
            for e in self.elem:
                is_elem,symbol = is_element(e,symbol=True)
                if e not in speclab_set:
                    speclab_set.add(e)
                    species_labels.append(e)
                #end if
                if symbol not in spec_set:
                    spec_set.add(symbol)
                    species.append(symbol)
                #end if
            #end for
            return species_labels,species
        #end if
    #end def ordered_species


    def order_by_species(self,folded=False):
        species        = []
        species_counts = []
        elem_indices   = []

        spec_set = set()
        for i in xrange(len(self.elem)):
            e = self.elem[i]
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
        self.reorder(elem_order)

        if folded and self.folded_structure!=None:
            self.folded_structure.order_by_species(folded)
        #end if

        return species,species_counts
    #end def order_by_species


    def reorder(self,order):
        order = array(order)
        self.elem = self.elem[order]
        self.pos  = self.pos[order]
    #end def reorder

    
    # find layers parallel to a particular cell face
    #   layers are found by scanning a window of width dtol along the axis and counting
    #     the number of atoms within the window.  window position w/ max number of atoms
    #     defines the layer.  layer distance is the window position.
    #   the resolution of the scan is determined by dbin
    #   (axis length)/dbin is the number of fine bins
    #   dtol/dbin is the number of fine bins in the moving (boxcar) window
    #   plot=True: plot the layer histogram (fine hist and moving average)
    #   composition=True: return the composition of each layer (count of each species)
    # returns an object containing indices of atoms in each layer by distance along axis
    #   example: structure w/ 3 layers of 4 atoms each at distances 3.0, 6.0, and 9.0 Angs.
    #   layers
    #     3.0 = [ 0, 1, 2, 3 ]
    #     6.0 = [ 4, 5, 6, 7 ]
    #     9.0 = [ 8, 9,10,11 ]
    def layers(self,axis=0,dtol=0.03,dbin=0.01,plot=False,composition=False):
        nbox = int(dtol/dbin)
        if nbox%2==0:
            nbox+=1
        #end if
        nwind = (nbox-1)/2
        s = self.copy()
        s.recenter()
        vaxis = s.axes[axis]
        daxis = norm(vaxis)
        naxis = vaxis/daxis
        dbin  = dtol/nbox
        nbins = int(ceil(daxis/dbin))
        dbin  = daxis/nbins
        dbins = daxis*(arange(nbins)+.5)/nbins
        dists = daxis*s.pos_unit()[:,axis]
        hist  = zeros((nbins,),dtype=int)
        boxhist = zeros((nbins,),dtype=int)
        ihist = obj()
        iboxhist = obj()
        index = 0
        for d in dists:
            ibin = int(floor(d/dbin))
            hist[ibin]+=1
            if not ibin in ihist:
                ihist[ibin] = []
            #end if
            ihist[ibin].append(index)
            index+=1
        #end for
        for ib in xrange(nbins):
            for i in xrange(ib-nwind,ib+nwind+1):
                n = hist[i%nbins]
                if n>0:
                    boxhist[ib]+=n
                    if not ib in iboxhist:
                        iboxhist[ib] = []
                    #end if
                    iboxhist[ib].extend(ihist[i%nbins])
                #end if
            #end for
        #end for
        peaks = []
        nlast=0
        for ib in xrange(nbins):
            n = boxhist[ib]
            if nlast==0 and n>0:
                pcur = []
                peaks.append(pcur)
            #end if
            if n>0:
                pcur.append(ib)
            #end if
            nlast = n
        #end for
        if boxhist[0]>0 and boxhist[-1]>0:
            peaks[0].extend(peaks[-1])
            peaks.pop()
        #end if
        layers = obj()
        ip = []
        for peak in peaks:
            ib = peak[boxhist[peak].argmax()]
            ip.append(ib)
            pindices = iboxhist[ib]
            ldist = dbins[ib] # distance is along an axis vector
            faxis = self.face_vectors()[axis]
            ldist = dot(ldist*naxis,faxis/norm(faxis))
            layers[ldist] = array(pindices,dtype=int)
        #end for
        if plot:
            plt.plot(dbins,boxhist,'b.-',label='boxcar histogram')
            plt.plot(dbins,hist,'r.-',label='fine histogram')
            plt.plot(dbins[ip],boxhist[ip],'rv',markersize=20)
            plt.show()
            plt.legend()
        #end if
        if not composition:
            return layers
        else:
            return layers,self.layer_composition(layers)
        #end if
    #end def layers


    def layer_composition(self,layers):
        lcomp = obj()
        for d,ind in layers.iteritems():
            comp = obj()
            elem = self.elem[ind]
            for e in elem:
                if e not in comp:
                    comp[e] = 1
                else:
                    comp[e] += 1
                #end if
            #end for
            lcomp[d]=comp
        #end for
        return lcomp
    #end def layer_composition


    def shells(self,identifiers,radii=None,exterior=False,cumshells=False,distances=False,dtol=1e-6):
        # get indices for 'core' and 'bulk'
        #   core is selected by identifiers, forms core for shells to be built around
        #   bulk is all atoms except for core
        if identifiers=='point_defects':
            if not 'point_defects' in self:
                self.error('requested shells around point defects, but structure has no point defects')
            #end if
            core = []
            for pd in self.point_defects:
                core.append(pd.center)
            #end for
            core = array(core)
            bulk_ind = self.locate(core,radii=dtol,exterior=True)
            core_ind = self.locate(bulk_ind,exterior=True)
            bulk = self.pos[bulk_ind]
        else:
            core_ind = self.locate(identifiers,radii,exterior)
            bulk_ind = self.locate(core_ind,exterior=True)
            core = self.pos[core_ind]
            bulk = self.pos[bulk_ind]
        #end if
        bulk_ind = array(bulk_ind,dtype=int)
        # build distance table between bulk and core
        dtable = self.distance_table(bulk,core)
        # find shortest distance for each bulk atom to any core atom and order by distance
        dist   = dtable.min(1)
        ind    = arange(len(bulk))
        order  = dist.argsort()
        dist   = dist[order]
        ind    = bulk_ind[ind[order]]
        # find shells around the core
        #   the closest atom to the core starts the first shell and defines a shell distance
        #   other atoms are in the shell if within dtol distance of the first atom
        #   otherwise a new shell is started
        ns = 0
        ds = -1
        shells = obj()
        shells[ns] = list(core_ind)  # first shell is all core atoms
        dshells = [0.]
        for n in xrange(len(dist)):
            if abs(dist[n]-ds)>dtol:
                shell = [ind[n]]   # new shell starts with single atom
                ns+=1
                shells[ns] = shell
                ds = dist[n]       # shell distance is distance of this atom from core
                dshells.append(ds)
            else:
                shell.append(ind[n])
            #end if
        #end for
        dshells = array(dshells,dtype=float)
        results = [shells]
        if cumshells:
            # assemble cumulative shells, ie cumshell[ns] = sum(shells[n],n=0 to ns)
            cumshells = obj()
            cumshells[0] = list(shells[0])
            for ns in xrange(1,len(shells)):
                cumshells[ns] = cumshells[ns-1]+shells[ns]
            #end for
            for ns,cshell in cumshells.iteritems():
                cumshells[ns] = array(cshell,dtype=int)
            #end for
            results.append(cumshells)
        #end if
        for ns,shell in shells.iteritems():
            shells[ns] = array(shell,dtype=int)
        if distances:
            results.append(dshells)
        #end if
        if len(results)==1:
            results = results[0]
        #end if
        return results
    #end def shells


    # find connected sets of atoms.
    #   indices is a list of atomic indices to consider (self.pos[indices] are their positions)
    #   atoms are considered connected if they are within rmax of each other
    #   order sets the maximum number of atoms in any connected graph
    #     order = 1 returns single atoms
    #     order = 2 returns dimers + order=1 results
    #     order = 3 returns trimers + order=2 results
    #     ...
    #   degree is explained w/ an example: a triangle of atoms 0,1,2  and a line of atoms 3,4,5 (3 & 5 are not neighbors)
    #     degree = False : returned object (cgraphs) has following structure:
    #       cgraphs[1] = [ (0,), (1,), (2,), (3,), (4,), (5,) ]  # first  order connected graphs (atoms)
    #       cgraphs[2] = [ (0,1), (0,2), (1,2), (3,4), (4,5) ]   # second order connected graphs (dimers)
    #       cgraphs[3] = [ (0,1,2), (3,4,5) ]                    # third  order connected graphs (trimers)
    #     degree = True : returned object (cgraphs) has following structure:
    #       cgraphs
    #         1      # first  order connected graphs (atoms)
    #           0    #   sum of vertex degrees is 0 (a single atom has no neighbors)
    #             (0,) = [ (0,), (1,), (2,), (3,), (4,), (5,) ]   # graphs with vertex degree (0,)
    #         2      # second order connected graphs (dimers)
    #           2    #   sum of vertex degrees is 2 (each atom is connected to 1 neighbor)
    #             (1,1) = [ (0,1), (0,2), (1,2), (3,4), (4,5) ]   # graphs with vertex degree (1,1)
    #         3      # third  order connected graphs (trimers)
    #           4    #   sum of vertex degrees is 4 (2 atoms have 1 neighbor and 1 atom has 2)
    #             (1,1,2) = [ (3,5,4) ]
    #           6    #   sum of vertex degrees is 6 (each atom is connected to 2 others)
    #             (2,2,2) = [ (0,1,2) ]           # graphs with vertex degree (2,2,2)  
    def connected_graphs(self,order,indices=None,rmax=None,nmax=None,voronoi=False,degree=False,site_maps=False,**spec_max):
        if indices is None:
            indices = arange(len(self.pos),dtype=int)
            pos = self.pos
        else:
            pos = self.pos[indices]
        #end if
        np = len(indices)
        neigh_table = []
        actual_indices = None
        if voronoi:
            actual_indices = True
            neighbors = self.voronoi_neighbors(indices,restrict=True,distance_ordered=False)
            for nilist in neighbors:
                neigh_table.append(nilist)
            #end for
        else:
            actual_indices = False
            elem = set(self.elem[indices])
            spec = set(spec_max.keys())
            if spec==elem or rmax!=None:
                None
            elif spec<elem and nmax!=None:
                for e in elem:
                    if e not in spec:
                        spec_max[e] = nmax
                    #end if
                #end for
            #end if
            # get neighbor table for subset of atoms specified by indices
            nt,dt = self.neighbor_table(pos,pos,distances=True)
            # determine how many neighbors to consider based on rmax (all are neighbors if rmax is None)
            nneigh = zeros((np,),dtype=int)
            if len(spec_max)>0:
                for n in xrange(np):
                    nneigh[n] = min(spec_max[self.elem[n]],len(nt[n]))
                #end for
            elif rmax is None:
                nneigh[:] = np
            else:
                nneigh = (dt<rmax).sum(1)                    
            #end if
            for i in xrange(np):
                neigh_table.append(nt[i,1:nneigh[i]])
            #end for
            del nt,dt,nneigh,elem,spec,rmax
        #end if
        neigh_table = array(neigh_table,dtype=int)
        # record which atoms are neighbors to each other
        neigh_pairs = set()
        if actual_indices:
            for i in xrange(np):
                for ni in neigh_table[i]:
                    neigh_pairs.add((i,ni))
                    neigh_pairs.add((ni,i))
                #end for
            #end for
        else:
            for i in xrange(np):
                for ni in neigh_table[i]:
                    ii = indices[i]
                    jj = indices[ni]
                    neigh_pairs.add((ii,jj))
                    neigh_pairs.add((jj,ii))
                #end for
            #end for
        #end if
        # find the connected graphs
        graphs_found = set()  # map to contain tuples of connected atom's indices
        cgraphs = obj()
        for o in range(1,order+1): # organize by order
            cgraphs[o] = []
        #end for
        if order>0:
            cg = cgraphs[1]
            for i in xrange(np):  # list of single atoms
                gi = (i,)              # graph indices
                cg.append(gi)          # add graph to graph list of order 1
                graphs_found.add(gi)   # add graph to set of all graphs
            #end for
            for o in range(2,order+1): # graphs of order o are found by adding all
                cglast = cgraphs[o-1]  # possible single neighbors to each graph of order o-1 
                cg     = cgraphs[o]
                for gilast in cglast:    # all graphs of order o-1
                    for i in gilast:       # all indices in each graph of order o-1
                        for ni in neigh_table[i]: # neighbors of selected atom in o-1 graph
                            gi = tuple(sorted(gilast+(ni,))) # new graph with neighbor added
                            if gi not in graphs_found and len(set(gi))==o: # add it if it is new and really is order o
                                graphs_found.add(gi)  # add graph to set of all graphs
                                cg.append(gi)         # add graph to graph list of order o
                            #end if
                        #end for
                    #end for
                #end for
            #end for
        #end if
        if actual_indices:
            for o,cg in cgraphs.iteritems():
                cgraphs[o] = array(cg,dtype=int)
            #end for
        else:
            # map indices back to actual atomic indices
            for o,cg in cgraphs.iteritems():
                cgmap = []
                for gi in cg:
                    #gi = array(gi)
                    gimap = tuple(sorted(indices[array(gi)]))
                    cgmap.append(gimap)
                #end for
                cgraphs[o] = array(sorted(cgmap),dtype=int)
            #end for
        #end if
        # reorganize the graph listing by cluster and vertex degree, if desired
        if degree:
            #degree_map = obj()
            cgraphs_deg = obj()
            for o,cg in cgraphs.iteritems():
                dgo = obj()
                cgraphs_deg[o] = dgo
                for gi in cg:
                    di = zeros((o,),dtype=int)
                    for m in xrange(o):
                        i = gi[m]
                        for n in xrange(m+1,o):
                            j = gi[n]
                            if (i,j) in neigh_pairs:
                                di[m]+=1
                                di[n]+=1
                            #end if
                        #end for
                    #end for
                    d = int(di.sum())
                    dorder = di.argsort()
                    di = tuple(di[dorder])
                    gi = tuple(array(gi)[dorder])
                    if not d in dgo:
                        dgo[d]=obj()
                    #end if
                    dgd = dgo[d]
                    if not di in dgd:
                        dgd[di] = []
                    #end if
                    dgd[di].append(gi)
                    #degree_map[gi] = d,di
                #end for
                for dgd in dgo:
                    for di,dgi in dgd.iteritems():
                        dgd[di]=array(sorted(dgi),dtype=int)
                    #end for
                #end for
            #end for
            cgraphs = cgraphs_deg
        #end if

        if not site_maps:
            return cgraphs
        else:
            cmaps = obj()
            if not degree:
                for order,og in cgraphs.iteritems():
                    cmap = obj()
                    for slist in og:
                        for s in slist:
                            if not s in cmap:
                                cmap[s] = obj()
                            #end if
                            cmap[s].append(slist)
                        #end for
                    #end for
                    cmaps[order] = cmap
                #end for
            else:
                for order,og in cgraphs.iteritems():
                    for total_degree,tg in og.iteritems():
                        for local_degree,lg in tg.iteritems():
                            cmap = obj()
                            for slist in lg:
                                n=0
                                for s in slist:
                                    d = local_degree[n]
                                    if not s in cmap:
                                        cmap[s] = obj()
                                    #end if
                                    if not d in cmap[s]:
                                        cmap[s][d] = obj()
                                    #end if
                                    cmap[s][d].append(slist)
                                    n+=1
                                #end for
                            #end for
                            cmaps.add_attribute_path((order,total_degree,local_degree),cmap)
                        #end for
                    #end for
                #end for
            #end if
            return cgraphs,cmaps
        #end if
    #end def connected_graphs


    # returns connected graphs that are rings up to the requested order
    #   rings are constructed by pairing lines that share endpoints
    #   all vertices of a ring have degree two
    def ring_graphs(self,order,**kwargs):
        # get all half order connected graphs
        line_order = order/2+order%2+1
        cgraphs = self.connected_graphs(line_order,degree=True,site_maps=False,**kwargs)
        # collect half order graphs that are lines
        lgraphs = obj()
        for o in range(2,line_order+1):
            total_degree  = 2*o-2
            vertex_degree = tuple([1,1]+(o-2)*[2])
            lg = None
            if o in cgraphs:
                cg = cgraphs[o]
                if total_degree in cg:
                    dg = cg[total_degree]
                    if vertex_degree in dg:
                        lg = dg[vertex_degree]
                    #end if
                #end if
            #end if
            if lg!=None:
                lg_end = obj()
                for gi in lg:
                    end_key = tuple(sorted(gi[0:2])) # end points
                    if end_key not in lg_end:
                        lg_end[end_key] = []
                    #end if
                    lg_end[end_key].append(tuple(gi))
                #end for
                lgraphs[o] = lg_end
            #end if
        #end for
        # contruct rings from lines that share endpoints
        rgraphs = obj()
        for o in range(3,order+1):
            o1 = o/2+1    # split half order for odd, same for even, 
            o2 = o1+o%2
            lg1 = lgraphs.get_optional(o1,None) # sets of half order lines
            lg2 = lgraphs.get_optional(o2,None)
            if lg1!=None and lg2!=None:
                rg = []
                rset = set()
                for end_key,llist1 in lg1.iteritems(): # list of lines sharing endpoints
                    if end_key in lg2:
                        llist2 = lg2[end_key]          # second list of lines sharing endpoints
                        for gi1 in llist1:             # combine line pairs into rings
                            for gi2 in llist2:
                                ri = tuple(sorted(set(gi1+gi2[2:]))) # ring indices
                                if ri not in rset and len(ri)==o:    # exclude repeated lines or rings
                                    rg.append(ri)
                                    rset.add(ri)
                                #end if
                            #end for
                        #end for
                    #end if
                #end for
                rgraphs[o] = array(sorted(rg),dtype=int)
            #end if
        #end for
        return rgraphs
    #end def ring_graphs


    # find the centroid of a set of points/atoms in min image convention
    def min_image_centroid(self,points=None,indices=None):
        if indices!=None:
            points = self.pos[indices]
        elif points is None:
            self.error('points or images must be provided to min_image_centroid')
        #end if
        p     = array(points,dtype=float)
        cprev = p[0]+1e99
        c     = p[0]
        while(norm(c-cprev)>1e-8):
            p = self.cell_image(p,center=c)
            cprev = c
            c = p.mean(axis=0)
        #end def min_image_centroid
        return c
    #end def min_image_centroid


    # find min image centroids of multiple sets of points/atoms
    def min_image_centroids(self,points=None,indices=None):
        cents = []
        if points!=None:
            for p in points:
                cents.append(self.min_image_centroid(p))
            #end for
        elif indices!=None:
            for ind in indices:
                cents.append(self.min_image_centroid(indices=ind))
            #end for
        else:
            self.error('points or images must be provided to min_image_centroid')
        #end if
        return array(cents,dtype=float)
    #end def min_image_centroids
    
    
    def min_image_vectors(self,points=None,points2=None,axes=None,pairs=True):
        if points is None:
            points = self.pos
        #end if
        if axes is None:
            axes  = self.axes
        #end if
        axinv = inv(axes)
        points = array(points)
        single = points.shape==(self.dim,)
        if single:
            points = [points]
        #end if
        if points2 is None:
            points2 = self.pos
        elif points2.shape==(self.dim,):
            points2 = [points2]
        #end if
        npoints  = len(points)
        npoints2 = len(points2)
        if pairs:
            vtable = empty((npoints,npoints2,self.dim),dtype=float)
            i=-1
            for p in points:
                i+=1
                j=-1
                for pp in points2:
                    j+=1
                    u = dot(pp-p,axinv)
                    vtable[i,j] = dot(u-floor(u+.5),axes)
                #end for
            #end for
            result = vtable
        else:
            if npoints!=npoints2:
                self.error('cannot create one to one minimum image vectors, point sets differ in length\n  npoints1 = {0}\n  npoints2 = {1}'.format(npoints,npoints2))
            #end if
            vectors = empty((npoints,self.dim),dtype=float)
            n = 0
            for p in points:
                pp = points2[n]
                u = dot(pp-p,axinv)
                vectors[n] = dot(u-floor(u+.5),axes)
                n+=1
            #end for
            result = vectors
        #end if
                
        return result
    #end def min_image_vectors


    def min_image_distances(self,points=None,points2=None,axes=None,vectors=False,pairs=True):
        vtable = self.min_image_vectors(points,points2,axes,pairs=pairs)
        rdim = len(vtable.shape)-1
        dtable = sqrt((vtable**2).sum(rdim))
        if not vectors:
            return dtable
        else:
            return dtable,vtable
        #end if
    #end def min_image_distances


    def distance_table(self,points=None,points2=None,axes=None,vectors=False):
        return self.min_image_distances(points,points2,axes,vectors)
    #end def distance_table


    def vector_table(self,points=None,points2=None,axes=None):
        return self.min_image_vectors(points,points2,axes)
    #end def vector_table

    
    def neighbor_table(self,points=None,points2=None,axes=None,distances=False,vectors=False):
        dtable,vtable = self.min_image_distances(points,points2,axes,vectors=True)
        ntable = empty(dtable.shape,dtype=int)
        for i in range(len(dtable)):
            ntable[i] = dtable[i].argsort()
        #end for
        results = [ntable]
        if distances:
            for i in range(len(dtable)):
                dtable[i] = dtable[i][ntable[i]]
            #end for
            results.append(dtable)
        #end if
        if vectors:
            for i in range(len(vtable)):
                vtable[i] = vtable[i][ntable[i]]
            #end for
            results.append(vtable)
        #end if
        if len(results)==1:
            results = results[0]
        #end if
        return results
    #end def neighbor_table


    def min_image_norms(self,points,norms):
        if isinstance(norms,int) or isinstance(norms,float):
            norms = [norms]
        #end if
        vtable = self.min_image_vectors(points)
        rdim = len(vtable.shape)-1
        nout = []
        for p in norms:
            nout.append( ((abs(vtable)**p).sum(rdim))**(1./p) )
        #end for
        if len(norms)==1:
            nout = nout[0]
        #end if
        return nout
    #end def min_image_norms


    # get all neighbors according to contacting voronoi polyhedra in PBC
    def voronoi_neighbors(self,indices=None,restrict=False,distance_ordered=True):
        if indices is None:
            indices = arange(len(self.pos))
        #end if
        indices = set(indices)
        # make a new version of this (small cell)
        sn = self.copy()
        sn.recenter()
        # tile a large cell periodically
        d = 3
        t = tuple(zeros((d,),dtype=int)+3)
        ss = sn.tile(t)
        ss.recenter(sn.center)
        # get nearest neighbor index pairs in the large cell
        neigh_pairs = voronoi_neighbors(ss.pos)
        # create a mapping from large to small indices
        large_to_small = 3**d*range(len(self.pos))
        # find the neighbor pairs in the small cell
        neighbors = obj()
        small_inds = set(ss.locate(sn.pos))
        for n in xrange(len(neigh_pairs)):
            i,j = neigh_pairs[n,:]
            if i in small_inds or j in small_inds: # pairs w/ at least one in cell image
                i = large_to_small[i]  # mapping to small cell indices
                j = large_to_small[j]
                if not restrict or (i in indices and j in indices): # restrict to orig index set
                    if not i in neighbors:
                        neighbors[i] = [j]
                    else:
                        neighbors[i].append(j)
                    #ned if
                    if not j in neighbors:
                        neighbors[j] = [i]
                    else:
                        neighbors[j].append(i)
                    #end if
                #end if
            #end if
        #end for
        # remove any duplicates and order by distance
        if distance_ordered:
            dt = self.distance_table()
            for i,ni in neighbors.iteritems():
                ni = array(list(set(ni)),dtype=int)
                di = dt[i,ni]
                order = di.argsort()
                neighbors[i] = ni[order]
            #end for
        else:  # just remove duplicates
            for i,ni in neighbors.iteritems():
                neighbors[i] = array(list(set(ni)),dtype=int)
            #end for
        #end if
        return neighbors
    #end def voronoi_neighbors


    # get nearest neighbors according to constrants (voronoi, max distance, coord. number)
    def nearest_neighbors(self,indices=None,rmax=None,nmax=None,restrict=False,voronoi=False,distances=False,**spec_max):
        if indices is None:
            indices = arange(len(self.pos))
        #end if
        elem = set(self.elem[indices])
        spec = set(spec_max.keys())
        if spec==elem or rmax!=None or voronoi:
            None
        elif spec<elem and nmax!=None:
            for e in elem:
                if e not in spec:
                    spec_max[e] = nmax
                #end if
            #end for
        else:
            self.error('must specify nmax for all species\n  species present: {0}\n  you only provided nmax for these species: {1}'.format(sorted(elem),sorted(spec)))
        #end if
        pos = self.pos[indices]
        if not restrict:
            pos2 = self.pos
        else:
            pos2 = pos
        #end if
        if voronoi:
            neighbors = self.voronoi_neighbors(indices=indices,restrict=restrict)
            dt = self.distance_table(pos,pos2)[:,1:]
        else:
            nt,dt = self.neighbor_table(pos,pos2,distances=True)
            dt=dt[:,1:]
            nt=nt[:,1:]
            neighbors = list(nt)
        #end if
        for i in xrange(len(indices)):
            neighbors[i] = indices[neighbors[i]]
        #end for
        dist = list(dt)
        if rmax is None:
            for i in xrange(len(indices)):
                nn = neighbors[i]
                dn = dist[i]
                e = self.elem[indices[i]]
                if e in spec_max:
                    smax = spec_max[e]
                    if len(nn)>smax:
                        neighbors[i] = nn[:smax]
                        dist[i]      = dn[:smax]
                    #end if
                #end if
            #end for
        else:
            for i in xrange(len(indices)):
                neighbors[i] = neighbors[i][dt[i]<rmax]
            #end for
        #end if
        if not distances:
            return neighbors
        else:
            return neighbors,dist
        #end if
    #end def nearest_neighbors


    # determine local chemical coordination limited by constraints
    def chemical_coordination(self,indices=None,nmax=None,rmax=None,restrict=False,voronoi=False,neighbors=False,distances=False,**spec_max):
        if indices is None:
            indices = arange(len(self.pos))
        #end if
        if not distances:
            neigh = self.nearest_neighbors(indices=indices,nmax=nmax,rmax=rmax,restrict=restrict,voronoi=voronoi,**spec_max)
        else:
            neigh,dist = self.nearest_neighbors(indices=indices,nmax=nmax,rmax=rmax,restrict=restrict,voronoi=voronoi,distances=True,**spec_max)
        #end if
        neigh_elem = []
        for i in xrange(len(indices)):
            neigh_elem.extend(self.elem[neigh[i]])
        #end for
        chem_key = tuple(sorted(set(neigh_elem)))
        chem_coord = zeros((len(indices),len(chem_key)),dtype=int)
        for i in xrange(len(indices)):
            counts = zeros((len(chem_key),),dtype=int)
            nn = list(self.elem[neigh[i]])
            for n in xrange(len(counts)):
                chem_coord[i,n] = nn.count(chem_key[n])
            #end for
        #end for
        chem_map = obj()
        i=0
        for coord in chem_coord:
            coord = tuple(coord)
            if not coord in chem_map:
                chem_map[coord] = [indices[i]]
            else:
                chem_map[coord].append(indices[i])
            #end if
            i+=1
        #end for
        for coord,ind in chem_map.iteritems():
            chem_map[coord] = array(ind,dtype=int)
        #end for
        results = [chem_key,chem_coord,chem_map]
        if neighbors:
            results.append(neigh)
        #end if
        if distances:
            results.append(dist)
        #end if
        return results
    #end def chemical_coordination


    def rcore_max(self,units=None):
        nt,dt = self.neighbor_table(self.pos,distances=True)
        d = dt[:,1]
        rcm = d.min()/2
        if units!=None:
            rcm = convert(rcm,self.units,units)
        #end if
        return rcm
    #end def rcore_max


    def cell_image(self,p,center=None):
        pos = array(p,dtype=float)
        if center is None:
            c = self.center.copy()
        else:
            c = array(center,dtype=float)
        #end if
        axes = self.axes
        axinv = inv(axes)
        for i in xrange(len(pos)):
            u = dot(pos[i]-c,axinv)
            pos[i] = dot(u-floor(u+.5),axes)+c
        #end for
        return pos
    #end def cell_image


    def center_distances(self,points,center=None):
        if center is None:
            c = self.center.copy()
        else:
            c = array(center,dtype=float)
        #end if        
        points = self.cell_image(points,center=c)
        for i in xrange(len(points)):
            points[i] -= c
        #end for
        return sqrt((points**2).sum(1))
    #end def center_distances


    def recenter(self,center=None):
        if center!=None:
            self.center=array(center)
        #end if
        pos = self.pos
        c = empty((1,self.dim),dtype=float)
        c[:] = self.center[:]
        axes = self.axes
        axinv = inv(axes)
        for i in xrange(len(pos)):
            u = dot(pos[i]-c,axinv)
            pos[i] = dot(u-floor(u+.5),axes)+c
        #end for
        self.recenter_k()
    #end def recenter


    def recorner(self):
        pos = self.pos
        axes = self.axes
        axinv = inv(axes)
        for i in range(len(pos)):
            u = dot(pos[i],axinv)
            pos[i] = dot(u-floor(u),axes)
        #end for
    #end def recorner

    
    def recenter_k(self,kpoints=None,kaxes=None,kcenter=None,remove_duplicates=False):
        use_self = kpoints==None
        if use_self:
            kpoints=self.kpoints
        #end if
        if kaxes==None:
            kaxes=self.kaxes
        #end if
        if len(kpoints)>0:
            axes = kaxes
            axinv = inv(axes)
            if kcenter is None:
                c = axes.sum(0)/2
            else:
                c = array(kcenter)
            #end if
            for i in range(len(kpoints)):
                u = dot(kpoints[i]-c,axinv)
                u -= floor(u+.5)
                u[abs(u-.5)<1e-12] -= 1.0
                u[abs(u   )<1e-12]  = 0.0
                kpoints[i] = dot(u,axes)+c
            #end for
            if remove_duplicates:
                inside = self.inside(kpoints,axes,c)
                kpoints  = kpoints[inside]
                nkpoints = len(kpoints)
                unique = empty((nkpoints,),dtype=bool)
                unique[:] = True
                nn = nearest_neighbors(1,kpoints)
                if nkpoints>1:
                    nn.shape = nkpoints,
                    dist = self.distances(kpoints,kpoints[nn])
                    tol = 1e-8
                    duplicates = arange(nkpoints)[dist<tol]
                    for i in duplicates:
                        if unique[i]:
                            for j in duplicates:
                                if sqrt(((kpoints[i]-kpoints[j])**2).sum(1))<tol:
                                    unique[j] = False
                                #end if
                            #end for
                        #end if
                    #end for
                #end if
                kpoints = kpoints[unique]
            #end if
        #end if
        if use_self:
            self.kpoints = kpoints
        else:
            return kpoints   
        #end if
    #end def recenter_k


    def inside(self,pos,axes=None,center=None,tol=1e-8,separate=False):
        if axes==None:
            axes=self.axes
        #end if
        if center==None:
            center=self.center
        #end if
        axes = array(axes)
        center = array(center)
        inside = []
        surface = []
        su = []
        axinv = inv(axes)
        for i in xrange(len(pos)):
            u = dot(pos[i]-center,axinv)
            umax = abs(u).max()
            if abs(umax-.5)<tol:
                surface.append(i)
                su.append(u)
            elif umax<.5:
                inside.append(i)
            #end if
        #end for
        npos,dim = pos.shape
        drange = range(dim)
        n = len(surface)
        i=0
        while i<n:
            j=i+1
            while j<n:
                du = abs(su[i]-su[j])
                match = False
                for d in drange:
                    match = match or abs(du[d]-1.)<tol
                #end for
                if match:
                    surface[j]=surface[-1]
                    surface.pop()
                    su[j]=su[-1]
                    su.pop()
                    n-=1
                else:
                    j+=1
                #end if
            #end while
            i+=1
        #end while
        if not separate:
            inside+=surface
            return inside
        else:
            return inside,surface
        #end if
    #end def inside


    def tile(self,*td,**kwargs):
        in_place           = kwargs.pop('in_place',False)
        magnetic_order     = kwargs.pop('magnetic_order',None)
        magnetic_primitive = kwargs.pop('magnetic_primitive',True)

        dim = self.dim
        if len(td)==1:
            if isinstance(td[0],int):
                tiling = dim*[td[0]]
            else:
                tiling = td[0]
            #end if
        else:
            tiling = td
        #end if
        tiling = array(tiling)

        matrix_tiling = tiling.shape == (dim,dim)

        tilematrix,tilevector = reduce_tilematrix(tiling)

        ncells = int(round( abs(det(tilematrix)) ))

        if ncells==1 and abs(tilematrix-identity(self.dim)).sum()<1e-1:
            if in_place:
                return self
            else:
                return self.copy()
            #end if
        #end if

        self.recenter()

        elem = array(ncells*list(self.elem))
        pos  = self.tile_points(self.pos,self.axes,tilematrix,tilevector)
        axes = dot(tilematrix,self.axes)

        center   = axes.sum(0)/2
        mag      = tile_magnetization(self.mag,tilevector,magnetic_order,magnetic_primitive)
        kaxes    = dot(inv(tilematrix.T),self.kaxes)
        kpoints  = array(self.kpoints)
        kweights = array(self.kweights)

        ts = self.copy()
        ts.center  = center
        ts.set_elem(elem)
        ts.axes    = axes
        ts.pos     = pos
        ts.mag     = mag
        ts.kaxes   = kaxes
        ts.kpoints = kpoints
        ts.kweights= kweights
        ts.background_charge = ncells*self.background_charge

        ts.recenter()
        ts.unique_kpoints()
        if self.folded_structure!=None:
            ts.tmatrix = dot(tilematrix,self.tmatrix)
            ts.folded_structure = self.folded_structure.copy()
        else:
            ts.tmatrix = tilematrix
            ts.folded_structure = self.copy()
        #end if

        if in_place:
            self.clear()
            self.transfer_from(ts)
            ts = self
        #end if

        return ts
    #end def tile


    def tile_points(self,points,axes,tilemat,tilevec=None):
        if tilevec is None:
            tilemat,tilevec = reduce_tilematrix(tilemat)
        #end if
        if not isinstance(points,ndarray):
            points = array(points)
        #end if
        if not isinstance(tilevec,ndarray):
            tilevec = array(tilevec)
        #end if
        if not isinstance(axes,ndarray):
            axes = array(axes)
        #end if
        t = tilevec
        ti = array(around(t),dtype=int)
        noninteger = abs(t-ti).sum()>1e-6
        if len(points.shape)==1:
            npoints,dim = len(points),1
        else:
            npoints,dim = points.shape
        #end if
        ntpoints = npoints*int(round( t.prod() ))
        if not noninteger:
            t = ti
            if ntpoints==0:
                tpoints = array([])
            else:
                tpoints  = empty((ntpoints,dim))
                ns=0
                ne=npoints
                for k in range(t[2]):
                    for j in range(t[1]):
                        for i in range(t[0]):
                            v = dot(array([[i,j,k]]),axes)
                            for d in range(dim):
                                tpoints[ns:ne,d] = points[:,d]+v[0,d]
                            #end for
                            ns+=npoints 
                            ne+=npoints
                        #end for
                    #end for
                #end for
            #end if
        else:
            if abs(ntpoints-int(around(ntpoints)))>1e-6:
                self.error('tiling vector does not correspond to an integer volume change\ntiling vector: {0}\nvolume change: {1}  {2}  {3}'.format(tilevec,tilevec.prod(),ntpoints,int(ntpoints)))
            #end if
            ntpoints = int(around(ntpoints))
            # round up to larger tiling
            #  +1 added for greater cell coverage
            #  add more if error below is tripped w/ fewer points than expected
            t = array(ceil(t),dtype=int)+1 
            # get the tiled points
            tpoints = self.tile_points(points,axes,tilemat,t)
            # remove any that are not unique
            taxes = dot(tilemat,axes)
            #tpoints,weights,pmap = self.unique_points(tpoints,taxes)
            tpoints,weights,pmap = self.unique_points_fast(tpoints,taxes)
            if len(tpoints)!=ntpoints:
                self.error('tiling by non-integer tiling vector failed\npoints expected after tiling: {0}\npoints resulted from tiling: {1}'.format(ntpoints,len(tpoints)))
            #end if
        #end if
        return tpoints
    #end def tile_points


    def opt_tilematrix(self,*args,**kwargs):
        return optimal_tilematrix(self.axes,*args,**kwargs)
    #end def opt_tilematrix


    def tile_opt(self,*args,**kwargs):
        Topt,ropt = self.opt_tilematrix(*args,**kwargs)
        return self.tile(Topt)
    #end def tile_opt


    def kfold(self,tiling,kpoints,kweights):
        if isinstance(tiling,int):
            tiling = self.dim*[tiling]
        #end if
        tiling = array(tiling)
        if tiling.shape==(self.dim,self.dim):
            tiling = tiling.T
        #end if
        tilematrix,tilevector = reduce_tilematrix(tiling)
        ncells = int(round( abs(det(tilematrix)) ))
        kp     = self.tile_points(kpoints,self.kaxes,tilematrix,tilevector)
        kw     = array(ncells*list(kweights),dtype=float)/ncells
        return kp,kw
    #end def kfold


    def get_primitive(self):
        if self.folded_structure is None:
            fs = self
        else:
            fs = self.folded_structure
            while fs.folded_structure!=None:
                fs = fs.folded_structure
            #end while
        #end if
        return fs
    #end def get_primitive


    def fold(self,small,*requests):
        self.error('fold needs a developers attention to make it equivalent with tile')
        if self.dim!=3:
            self.error('fold is currently only implemented for 3 dimensions')
        #end if
        self.recenter_k()
        corners = []
        ndim = len(small.axes)
        imin = empty((ndim,),dtype=int)
        imax = empty((ndim,),dtype=int)
        imin[:] =  1000000
        imax[:] = -1000000
        axinv  = inv(self.kaxes)
        center = self.kaxes.sum(0)/2 
        c = empty((1,3))
        for k in -1,2:
            for j in -1,2:
                for i in -1,2:
                    c[:] = i,j,k
                    c = dot(c,small.kaxes)
                    u = dot(c-center,axinv)
                    for d in range(ndim):
                        imin[d] = min(int(floor(u[0,d])),imin[d])
                        imax[d] = max(int(ceil(u[0,d])),imax[d])
                    #end for
                #end for
            #end for
        #end for

        axes = small.kaxes
        axinv = inv(small.kaxes)

        center = small.kaxes.sum(0)/2
        nkpoints = len(self.kpoints)
        kindices = []
        kpoints  = []
        shift = empty((ndim,))
        kr = range(nkpoints)
        for k in range(imin[2],imax[2]+1):
            for j in range(imin[1],imax[1]+1):
                for i in range(imin[0],imax[0]+1):
                    for n in kr:
                        shift[:] = i,j,k
                        shift = dot(shift,self.kaxes)
                        kp = self.kpoints[n]+shift
                        u = dot(kp-center,axinv)
                        if abs(u).max()<.5+1e-10:
                            kindices.append(n)
                            kpoints.append(kp)
                        #end if
                    #end for
                #end for
            #end for
        #end for
        kindices = array(kindices)
        kpoints  = array(kpoints)
        inside = self.inside(kpoints,axes,center)
        kindices = kindices[inside]
        kpoints  = kpoints[inside]

        small.kpoints = kpoints
        small.recenter_k()
        kpoints = array(small.kpoints)
        if len(requests)>0:
            results = []
            for request in requests:
                if request=='kmap':
                    kmap = obj()
                    for k in self.kpoints:
                        kmap[tuple(k)] = []
                    #end for
                    for i in range(len(kpoints)):
                        kp = tuple(self.kpoints[kindices[i]])
                        kmap[kp].append(array(kpoints[i]))
                    #end for
                    for kl,ks in kmap.iteritems():
                        kmap[kl] = array(ks)
                    #end for
                    res = kmap
                elif request=='tilematrix':
                    res = self.tilematrix(small)
                else:
                    self.error(request+' is not a recognized input to fold')
                #end if
                results.append(res)
            #end if
            return results
        #end if
    #end def fold


    def tilematrix(self,small=None,tol=1e-6,status=False):
        if small==None:
            if self.folded_structure!=None:
                small = self.folded_structure
            else:
                return identity(self.dim,dtype=int)
            #end if
        #end if
        tm = dot(self.axes,inv(small.axes))
        tilemat = array(around(tm),dtype=int)
        error = abs(tilemat-tm).sum()
        non_integer_elements = error > tol
        if status:
            return tilemat,not non_integer_elements
        else:
            if non_integer_elements:
                self.error('large cell cannot be constructed as an integer tiling of the small cell\nlarge cell axes:\n'+str(self.axes)+'\nsmall cell axes:  \n'+str(small.axes)+'\nlarge/small:\n'+str(self.axes/small.axes)+'\ntiling matrix:\n'+str(tm)+'\nintegerized tiling matrix:\n'+str(tilemat)+'\nerror: '+str(error)+'\ntolerance: '+str(tol))
            #end if
            return tilemat
        #end if
    #end def tilematrix
            

    def add_kpoints(self,kpoints,kweights=None,unique=False):
        if kweights is None:
            kweights = ones((len(kpoints),))
        #end if
        self.kpoints  = append(self.kpoints,kpoints,axis=0)
        self.kweights = append(self.kweights,kweights)
        if unique:
            self.unique_kpoints()
        #end if
        self.recenter_k() #added because qmcpack cannot handle kpoints outside the box
        if self.folded_structure!=None:
            kp,kw = self.kfold(self.tmatrix,kpoints,kweights)
            self.folded_structure.add_kpoints(kp,kw,unique=unique)
        #end if
    #end def add_kpoints


    def clear_kpoints(self):
        self.kpoints  = empty((0,self.dim))
        self.kweights = empty((0,))
        if self.folded_structure!=None:
            self.folded_structure.clear_kpoints()
        #end if
    #end def clear_kpoints


    def add_kmesh(self,kgrid,kshift=None,unique=False):
        self.add_kpoints(kmesh(self.kaxes,kgrid,kshift),unique=unique)
    #end def add_kmesh

    
    def kpoints_unit(self,kpoints=None):
        if kpoints is None:
            kpoints = self.kpoints
        #end if
        return dot(kpoints,inv(self.kaxes))
    #end def kpoints_unit


    def kpoints_reduced(self):
        return self.kpoints*self.scale/(2*pi)
    #end def kpoints_reduced


    def inversion_symmetrize_kpoints(self,tol=1e-10,folded=False):
        kp    = self.kpoints
        kaxes = self.kaxes
        ntable,dtable = self.neighbor_table(kp,-kp,kaxes,distances=True)
        pairs = set()
        keep = empty((len(kp),),dtype=bool)
        keep[:] = True
        for i in range(len(dtable)):
            if keep[i] and dtable[i,0]<tol:
                j = ntable[i,0]
                if j!=i and keep[j]:
                    keep[j] = False
                    self.kweights[i] += self.kweights[j]
                #end if
            #end if
        #end for
        self.kpoints  = self.kpoints[keep]
        self.kweights = self.kweights[keep]
        if folded and self.folded_structure!=None:
            self.folded_structure.inversion_symmetrize_kpoints(tol)
        #end if
    #end def inversion_symmetrize_kpoints


    def unique_points(self,points,axes,weights=None,tol=1e-10):
        pmap = obj()
        npoints = len(points)
        if npoints>0:
            if weights is None:
                weights = ones((npoints,),dtype=int)
            #end if
            ntable,dtable = self.neighbor_table(points,points,axes,distances=True)
            keep = empty((npoints,),dtype=bool)
            keep[:] = True
            pmo = obj()
            for i in xrange(npoints):
                if keep[i]:
                    pm = []
                    jn=0
                    while jn<npoints and dtable[i,jn]<tol:
                        j = ntable[i,jn]
                        pm.append(j)
                        if j!=i and keep[j]:
                            keep[j] = False
                            weights[i] += weights[j]
                        #end if
                        jn+=1
                    #end while
                    pmo[i] = set(pm)
                #end if
            #end for
            points  = points[keep]
            weights = weights[keep]
            j=0
            for i in xrange(len(keep)):
                if keep[i]:
                    pmap[j] = pmo[i]
                    j+=1
                #end if
            #end for
        #end if
        return points,weights,pmap
    #end def unique_points


    def unique_points_fast(self,points,axes,weights=None,tol=1e-10):
        # use an O(N) cell table instead of an O(N^2) neighbor table
        pmap = obj()
        points = array(points)
        axes   = array(axes)
        npoints = len(points)
        if npoints>0:
            if weights is None:
                weights = ones((npoints,),dtype=int)
            else:
                weights = array(weights)
            #end if
            keep = ones((npoints,),dtype=bool)
            # place all the points in the box, converted to unit coords
            upoints = array(points)
            axinv = inv(axes)
            for i in range(len(points)):
                u = dot(points[i],axinv)
                upoints[i] = u-floor(u)
            #end for
            # create an integer array of cell indices
            axmax = -1.0
            for a in axes:
                axmax = max(axmax,norm(a))
            #end for
            #   make an integer space corresponding to 1e-7 self.units spatial resolution
            cmax = uint64(1e7)*uint64(ceil(axmax)) 
            ipoints = array(around(cmax*upoints),dtype=uint64)
            ipoints[ipoints==cmax] = 0 # make the outer boundary the same as the inner boundary
            # load the cell table with point indices
            #   points in the same cell are identical
            ctable = obj()
            i=0
            for ip in ipoints:
                ip = tuple(ip)
                if ip not in ctable:
                    ctable[ip] = i
                    pmap[i] = [i]
                else:
                    j = ctable[ip]
                    keep[i] = False
                    weights[j] += weights[i]
                    pmap[j].append(i)
                #end if
                i+=1
            #end for
            points  = points[keep]
            weights = weights[keep]
        #end if
        return points,weights,pmap
    #end def unique_points_fast


    def unique_positions(self,tol=1e-10,folded=False):
        pos,weights,pmap = self.unique_points(self.pos,self.axes)
        if len(pos)!=len(self.pos):
            self.pos = pos
        #end if
        if folded and self.folded_structure!=None:
            self.folded_structure.unique_positions(tol)
        #end if
        return pmap
    #end def unique_positions

        
    def unique_kpoints(self,tol=1e-10,folded=False):
        kmap = obj()
        kp   = self.kpoints
        if len(kp)>0:
            kaxes = self.kaxes
            ntable,dtable = self.neighbor_table(kp,kp,kaxes,distances=True)
            npoints = len(kp)
            keep = empty((len(kp),),dtype=bool)
            keep[:] = True
            kmo = obj()
            for i in xrange(npoints):
                if keep[i]:
                    km = []
                    jn=0
                    while jn<npoints and dtable[i,jn]<tol:
                        j = ntable[i,jn]
                        km.append(j)
                        if j!=i and keep[j]:
                            keep[j] = False
                            self.kweights[i] += self.kweights[j]
                        #end if
                        jn+=1
                    #end while
                    kmo[i] = set(km)
                #end if
            #end for
            self.kpoints  = self.kpoints[keep]
            self.kweights = self.kweights[keep]
            j=0
            for i in xrange(len(keep)):
                if keep[i]:
                    kmap[j] = kmo[i]
                    j+=1
                #end if
            #end for
        #end if
        if folded and self.folded_structure!=None:
            self.folded_structure.unique_kpoints(tol)
        #end if
        return kmap
    #end def unique_kpoints


    def kmap(self):
        kmap = None
        if self.folded_structure!=None:
            fs = self.folded_structure
            self.kpoints  = array(fs.kpoints)
            self.kweights = array(fs.kweights)
            kmap = self.unique_kpoints()
        #end if
        return kmap
    #end def kmap


    def select_twist(self,selector='smallest',tol=1e-6):
        index = None
        invalid_selector = False
        if isinstance(selector,str):
            if selector=='smallest':
                index = (self.kpoints**2).sum(1).argmin()
            elif selector=='random':
                index = randint(0,len(self.kpoints)-1)
            else:
                invalid_selector = True
            #end if
        elif isinstance(selector,(tuple,list,ndarray)):
            ku_sel = array(selector,dtype=float)
            n = 0
            for ku in self.kpoints_unit():
                if norm(ku-ku_sel)<tol:
                    index = n
                    break
                #end if
                n+=1
            #end for
            if index is None:
                self.error('cannot identify twist number\ntwist requested: {0}\ntwists present: {1}'.format(ku_sel,sorted([tuple(k) for k in self.kpoints_unit()])))
            #end if
        else:
            invalid_selector = True
        #end if
        if invalid_selector:
            self.error('cannot identify twist number\ninvalid selector provided: {0}\nvalid string inputs for selector: smallest, random\nselector can also be a length 3 tuple, list or array (a twist vector)'.format(selector))
        #end if
        return index
    #end def select_twist


    def fold_pos(self,large,tol=0.001):
        vratio = large.volume()/self.volume()
        if abs(vratio-int(around(vratio)))>1e-6:
            self.error('cannot fold positions from large cell into current one\nlarge cell volume is not an integer multiple of the current one\nlarge cell volume: {0}\ncurrent cell volume: {1}\nvolume ratio: {2}'.format(large.volume(),self.volume(),vratio))
        T,success = large.tilematrix(self,status=True)
        if not success:
            self.error('cannot fold positions from large cell into current one\ncells are related by non-integer tilematrix')
        #end if
        nnearest = int(around(vratio))
        self.elem = large.elem.copy() 
        self.pos  = large.pos.copy()
        self.recenter()
        nt,dt = self.neighbor_table(distances=True)
        nt = nt[:,:nnearest]
        dt = dt[:,:nnearest]
        if dt.ravel().max()>tol:
            self.error('cannot fold positions from large cell into current one\npositions of equivalent atoms are further apart than the tolerance\nmax distance encountered: {0}\ntolerance: {1}'.format(dt.ravel().max(),tol))
        #end if
        counts = zeros((len(self.pos),),dtype=int)
        for n in nt.ravel():
            counts[n] += 1
        #end for
        if (counts!=nnearest).any():
            self.error('cannot fold positions from large cell into current one\neach atom must have {0} equivalent positions\nsome atoms found with the following equivalent position counts: {1}'.format(nnearest,counts[counts!=nnearest]))
        #end if
        ind_visited = set()
        neigh_map = obj()
        keep = []
        n=0
        for nset in nt:
            if n not in ind_visited:
                neigh_map[n] = nset
                keep.append(n)
                for ind in nset:
                    ind_visited.add(ind)
                #end for
            #end if
            n+=1
        #end for
        if len(ind_visited)!=len(self.pos):
            self.error('cannot fold positions from large cell into current one\nsome equivalent atoms could not be identified')
        #end if
        new_elem = []
        new_pos  = []
        for n in keep:
            nset = neigh_map[n]
            elist = list(set(self.elem[nset]))
            if len(elist)!=1:
                self.error('cannot fold positions from large cell into current one\nspecies of some equivalent atoms do not match')
            #end if
            new_elem.append(elist[0])
            new_pos.append(self.pos[nset].mean(0))
        #end for
        self.set_elem(new_elem)
        self.set_pos(new_pos)
    #end def fold_pos


    def pos_unit(self,pos=None):
        if pos is None:
            pos = self.pos
        #end if
        return dot(pos,inv(self.axes))
    #end def pos_unit


    def pos_to_cartesian(self):
        self.pos = dot(self.pos,self.axes)
    #end def pos_to_cartesian


    def at_Gpoint(self):
        kpu = self.kpoints_unit()
        kg = array([0,0,0])
        return len(kpu)==1 and norm(kg-kpu[0])<1e-6
    #end def at_Gpoint


    def at_Lpoint(self):
        kpu = self.kpoints_unit()
        kg = array([.5,.5,.5])
        return len(kpu)==1 and norm(kg-kpu[0])<1e-6
    #end def at_Lpoint


    def at_real_kpoint(self):
        kpu = 2*self.kpoints_unit()
        return len(kpu)==1 and abs(kpu-around(kpu)).sum()<1e-6
    #end def at_real_kpoint


    def bonds(self,neighbors,vectors=False):
        if self.dim!=3:
            self.error('bonds is currently only implemented for 3 dimensions')
        #end if
        natoms,dim = self.pos.shape
        centers = empty((natoms,neighbors,dim))
        distances = empty((natoms,neighbors))
        vect      = empty((natoms,neighbors,dim))
        t = self.tile((3,3,3))
        t.recenter(self.center)
        nn = nearest_neighbors(neighbors+1,t.pos,self.pos)
        for i in range(natoms):
            ii = nn[i,0]
            n=0
            for jj in nn[i,1:]:
                p1 = t.pos[ii]
                p2 = t.pos[jj]
                centers[i,n,:] = (p1+p2)/2
                distances[i,n]= sqrt(((p1-p2)**2).sum())
                vect[i,n,:] = p2-p1
                n+=1
            #end for
        #end for
        sn = self.copy()
        nnr = nn[:,1:].ravel()
        sn.set_elem(t.elem[nnr])
        sn.pos  = t.pos[nnr]
        sn.recenter()
        indices = self.locate(sn.pos)
        indices = indices.reshape(natoms,neighbors)
        if not vectors:
            return indices,centers,distances
        else:
            return indices,centers,distances,vect
        #end if
    #end def bonds

        
    def displacement(self,reference,map=False):
        if self.dim!=3:
            self.error('displacement is currently only implemented for 3 dimensions')
        #end if
        ref = reference.tile((3,3,3))
        ref.recenter(reference.center)
        rmap = array(3**3*range(len(reference.pos)),dtype=int)
        nn = nearest_neighbors(1,ref.pos,self.pos).ravel()
        displacement = self.pos - ref.pos[nn]
        if not map:
            return displacement
        else:
            return displacement,rmap[nn]
        #end if
    #end def displacement


    def scalar_displacement(self,reference):
        return sqrt((self.displacement(reference)**2).sum(1))
    #end def scalar_displacement

    
    def distortion(self,reference,neighbors):
        if self.dim!=3:
            self.error('distortion is currently only implemented for 3 dimensions')
        #end if
        if reference.volume()/self.volume() < 1.1:
            ref = reference.tile((3,3,3))
            ref.recenter(reference.center)
        else:
            ref = reference
        #end if
        rbi,rbc,rbl,rbv =  ref.bonds(neighbors,vectors=True)
        sbi,sbc,sbl,sbv = self.bonds(neighbors,vectors=True)
        nn = nearest_neighbors(1,reference.pos,self.pos).ravel()
        distortion = empty(sbv.shape)
        magnitude  = empty((len(self.pos),))
        for i in range(len(self.pos)):
            ir = nn[i]
            bonds  = sbv[i]
            rbonds = rbv[ir]
            ib  = empty((neighbors,),dtype=int)
            ibr = empty((neighbors,),dtype=int)
            r  = range(neighbors)
            rr = range(neighbors)
            for n in range(neighbors):
                mindist = 1e99
                ibmin  = -1
                ibrmin = -1
                for nb in r:
                    for nbr in rr:
                        d = norm(bonds[nb]-rbonds[nbr])
                        if d<mindist:
                            mindist=d
                            ibmin=nb
                            ibrmin=nbr
                        #end if
                    #end for
                #end for
                ib[n]=ibmin
                ibr[n]=ibrmin
                r.remove(ibmin)
                rr.remove(ibrmin)
                #end for
            #end for
            d = bonds[ib]-rbonds[ibr]
            distortion[i] = d
            magnitude[i] = (sqrt((d**2).sum(axis=1))).sum()
        #end for
        return distortion,magnitude
    #end def distortion


    def bond_compression(self,reference,neighbors):
        ref = reference
        rbi,rbc,rbl =  ref.bonds(neighbors)
        sbi,sbc,sbl = self.bonds(neighbors)
        bondlen = rbl.mean()
        return abs(1.-sbl/bondlen).max(axis=1)
    #end def bond_compression


    def boundary(self,dims=(0,1,2),dtol=1e-6):
        dim_eff = len(dims)
        natoms,dim = self.pos.shape
        bdims = array(dim*[False])
        for d in dims:
            bdims[d] = True
        #end for
        p = self.pos[:,bdims]
        indices = convex_hull(p,dim_eff,dtol)
        return indices
    #end def boundary


    def embed(self,small,dims=(0,1,2),dtol=1e-6,utol=1e-6):
        small = small.copy()
        small.recenter()
        center = array(self.center)
        self.recenter(small.center)
        bind = small.boundary(dims,dtol)
        bpos = small.pos[bind]
        belem= small.elem[bind]
        nn = nearest_neighbors(1,self.pos,bpos).ravel()
        mpos = self.pos[nn]
        dr = (mpos-bpos).mean(0)
        for i in xrange(len(bpos)):
            bpos[i]+=dr
        #end for
        dmax = sqrt(((mpos-bpos)**2).sum(1)).max()
        for i in xrange(len(small.pos)):
            small.pos[i]+=dr
        #end for
        ins,surface = small.inside(self.pos,tol=utol,separate=True)
        replaced = empty((len(self.pos),),dtype=bool)
        replaced[:] = False
        inside = replaced.copy()
        inside[ins] = True
        nn = nearest_neighbors(1,self.pos,small.pos).ravel()
        elist = list(self.elem)
        plist = list(self.pos)
        pos  = small.pos
        elem = small.elem
        for i in xrange(len(pos)):
            n = nn[i]
            if not replaced[n]:
                elist[n] = elem[i]
                plist[n] = pos[i]
                replaced[n] = True
            else:
                elist.append(elem[i])
                plist.append(pos[i])
            #end if
        #end for
        remove = arange(len(self.pos))[inside & logical_not(replaced)]
        remove.sort()
        remove = flipud(remove)
        for i in remove:
            elist.pop(i)
            plist.pop(i)
        #end for
        self.set_elem(elist)
        self.pos  = array(plist)
        self.recenter(center)
        return dmax
    #end def embed


    def shell(self,cell,neighbors,direction='in'):
        if self.dim!=3:
            self.error('shell is currently only implemented for 3 dimensions')
        #end if
        dd = {'in':equate,'out':negate}
        dir = dd[direction]
        natoms,dim=self.pos.shape
        ncells=3**3
        ntile = ncells*natoms
        pos = empty((ntile,dim))
        ind = empty((ntile,),dtype=int)
        oind = range(natoms)
        for nt in range(ncells):
            n=nt*natoms
            ind[n:n+natoms]=oind[:]
            pos[n:n+natoms]=self.pos[:]
        #end for
        nt=0
        for k in -1,0,1:
            for j in -1,0,1:
                for i in -1,0,1:
                    iv = array([[i,j,k]])
                    v = dot(iv,self.axes)
                    for d in range(dim):
                        ns = nt*natoms
                        ne = ns+natoms
                        pos[ns:ne,d] += v[0,d]
                    #end for
                    nt+=1
                #end for
            #end for
        #end for
        
        inside = empty(ntile,)
        inside[:]=False
        ins = cell.inside(pos)
        inside[ins]=True

        iishell = set()
        nn = nearest_neighbors(neighbors,pos)
        for ii in range(len(nn)):
            for jj in nn[ii]:
                in1 = inside[ii]
                in2 = inside[jj]
                if dir(in1 and not in2):
                    iishell.add(ii)
                #end if
                if dir(in2 and not in1):
                    iishell.add(jj)
                #end if
            #end if
        #end if
        ishell = ind[list(iishell)]
        return ishell
    #end def shell


    def interpolate(self,other,images,min_image=True,recenter=True,match_com=False,chained=False):
        s1 = self.copy()
        s2 = other.copy()
        s1.remove_folded()
        s2.remove_folded()
        if s2.units!=s1.units:
            s2.change_units(s1.units)
        #end if
        if (s1.elem!=s2.elem).any():
            self.error('cannot interpolate structures, atoms do not match\n  atoms1: {0}\n  atoms2: {1}'.format(s1.elem,s2.elem))
        #end if
        structures = []
        npath = images+2
        c1   = s1.center
        c2   = s2.center
        ax1  = s1.axes
        ax2  = s2.axes
        pos1 = s1.pos
        pos2 = s2.pos
        min_image &= abs(ax1-ax2).max()<1e-6
        if min_image:
            dp = self.min_image_vectors(pos1,pos2,ax1,pairs=False)
            pos2 = pos1 + dp
        #end if
        if match_com:
            com1 = pos1.mean(axis=0)
            com2 = pos2.mean(axis=1)
            dcom = com1-com2
            for n in xrange(len(pos2)):
                pos2[n] += dcom
            #end for
            if chained:
                other.pos = pos2
            #end if
        #end if
        for n in xrange(npath):
            f1 = 1.-float(n)/(npath-1)
            f2 = 1.-f1
            center = f1*c1   + f2*c2
            axes   = f1*ax1  + f2*ax2
            pos    = f1*pos1 + f2*pos2
            s = s1.copy()
            s.reset_axes(axes)
            s.center = center
            s.pos    = pos
            if recenter:
                s.recenter()
            #end if
            structures.append(s)
        #end for
        return structures
    #end def interpolate


    # returns madelung potential constant v_M
    #   see equation 7 in PRB 78 125106 (2008) 
    def madelung(self,axes=None,tol=1e-10):
        if self.dim!=3:
            self.error('madelung is currently only implemented for 3 dimensions')
        #end if
        if axes is None:
            a = self.axes.T.copy()
        else:
            a = axes.T.copy()
        #end if
        if self.units!='B':
            a = convert(a,self.units,'B')
        #end if
        volume = abs(det(a))
        b = 2*pi*inv(a).T
        rconv = 8*(3.*volume/(4*pi))**(1./3)
        kconv = 2*pi/rconv
        gconst = -1./(4*kconv**2)
        vmc = -pi/(kconv**2*volume)-2*kconv/sqrt(pi)
        
        nshells = 20
        vshell = [0.]
        p = Sobj()
        m = Sobj()
        for n in range(1,nshells+1):
            i = mgrid[-n:n+1,-n:n+1,-n:n+1]
            i = i.reshape(3,(2*n+1)**3)
            R = sqrt((dot(a,i)**2).sum(0))
            G2 = (dot(b,i)**2).sum(0)
            
            izero = n + n*(2*n+1) + n*(2*n+1)**2

            p.R  = R[0:izero]
            p.G2 = G2[0:izero]
            m.R  = R[izero+1:]
            m.G2 = G2[izero+1:]
            domains = [p,m]
            vshell.append(0.)
            for d in domains:
                vshell[n] += (erfc(kconv*d.R)/d.R).sum() + 4*pi/volume*(exp(gconst*d.G2)/d.G2).sum()
            #end for
            if abs(vshell[n]-vshell[n-1])<tol:
                break
            #end if
        #end for
        vm = vmc + vshell[-1]

        if axes is None:
            self.Vmadelung = vm
        #end if
        return vm
    #end def madelung


    def makov_payne(self,q=1,eps=1.0,units='Ha',order=1):
        if order!=1:
            self.error('Only first order Makov-Payne correction is currently supported.')
        #end if
        if 'Vmadelung' not in self:
            vm = self.madelung()
        else:
            vm = self.Vmadelung
        #end if
        mp = -0.5*q**2*vm/eps
        if units!='Ha':
            mp = convert(mp,'Ha',units)
        #end if
        return mp
    #end def makov_payne


    def read(self,filepath,format=None,elem=None,block=None,grammar='1.1',cell='prim',contents=False):
        if os.path.exists(filepath):
            path,file = os.path.split(filepath)
            if format is None:
                if '.' in file:
                    name,format = file.rsplit('.',1)
                elif file.lower().endswith('poscar'):
                    format = 'poscar'
                else:
                    self.error('file format could not be determined\nunrecognized file: {0}'.format(filepath))
                #end if
            #end if
        elif not contents:
            self.error('file does not exist: {0}'.format(filepath))
        #end if
        if format is None:
            self.error('file format must be provided')
        #end if
        format = format.lower()
        if format=='xyz':
            self.read_xyz(filepath)
        elif format=='xsf':
            self.read_xsf(filepath)
        elif format=='poscar':
            self.read_poscar(filepath,elem=elem)
        elif format=='cif':
            self.read_cif(filepath,block=block,grammar=grammar,cell=cell)
        elif format=='fhi-aims':
            self.read_fhi_aims(filepath)
        else:
            self.error('cannot read structure from file\nunsupported file format: {0}'.format(format))
        #end if
        if self.has_axes():
            self.set_bconds('ppp')
        #end if
    #end def read


    def read_xyz(self,filepath):
        elem = []
        pos  = []
        if os.path.exists(filepath):
            lines = open(filepath,'r').read().splitlines()
        else:
            lines = filepath.splitlines() # "filepath" is file contents
        #end if
        ntot = 1000000
        natoms = 0
        for l in lines:
            ls = l.strip()
            if ls.isdigit():
                ntot = int(ls)
            #end if
            tokens = ls.split()
            if len(tokens)==4:
                elem.append(tokens[0])
                pos.append(array(tokens[1:],float))
                natoms+=1
                if natoms==ntot:
                    break
                #end if
            #end if
        #end for
        self.dim   = 3
        self.set_elem(elem)
        self.pos   = array(pos)
        self.units = 'A'
    #end def read_xyz


    def read_xsf(self,filepath):
        if os.path.exists(filepath):
            f = XsfFile(filepath)
        else:
            f = XsfFile()
            f.read_text(filepath) # "filepath" is file contents
        #end if
        elem = []
        for n in f.elem:
            elem.append(pt.simple_elements[n].symbol)
        #end for
        self.dim   = 3
        self.units = 'A'
        self.reset_axes(f.primvec)
        self.set_elem(elem)
        self.pos = f.pos
    #end def read_xsf


    def read_poscar(self,filepath,elem=None):
        if os.path.exists(filepath):
            lines = open(filepath,'r').read().splitlines()
        else:
            lines = filepath.splitlines()  # "filepath" is file contents
        #end if
        nlines = len(lines)
        min_lines = 8
        if nlines<min_lines:
            self.error('POSCAR file must have at least {0} lines\n  only {1} lines found'.format(min_lines,nlines))
        #end if
        dim = 3
        scale = float(lines[1].strip())
        axes = empty((dim,dim))
        axes[0] = array(lines[2].split(),dtype=float)
        axes[1] = array(lines[3].split(),dtype=float)
        axes[2] = array(lines[4].split(),dtype=float)
        if scale<0.0:
            scale = abs(scale)/det(axes)
        #end if
        axes = scale*axes
        tokens = lines[5].split()
        if tokens[0].isdigit():
            counts = array(tokens,dtype=int)
            if elem is None:
                self.error('variable elem must be provided to read_poscar() to assign atomic species to positions for POSCAR format')
            elif len(elem)!=len(counts):
                self.error('one elem must be given for each element count in the POSCAR file\n  number of elem counts: {0}\n  number of elem given: {1}'.format(len(counts),len(elem)))
            #end if
            lcur = 6
        else:
            elem   = tokens
            counts = array(lines[6].split(),dtype=int)
            lcur = 7
        #end if
        species = elem
        # relabel species that have multiple occurances
        sset = set(species)
        for spec in sset:
            if species.count(spec)>1:
                cnt=0
                for n in range(len(species)):
                    specn = species[n]
                    if specn==spec:
                        cnt+=1
                        species[n] = specn+str(cnt)
                    #end if
                #end for
            #end if
        #end for
        elem = []
        for i in range(len(counts)):
            elem.extend(counts[i]*[species[i]])
        #end for
        self.dim = dim
        self.units = 'A'
        self.reset_axes(axes)

        if lcur<len(lines) and len(lines[lcur])>0:
            c = lines[lcur].lower().strip()[0]
            lcur+=1
        else:
            return
        #end if
        selective_dynamics = c=='s'
        if selective_dynamics: # Selective dynamics
            if lcur<len(lines) and len(lines[lcur])>0:
                c = lines[lcur].lower().strip()[0]
                lcur+=1
            else:
                return
            #end if
        #end if
        cartesian = c=='c' or c=='k'
        npos = counts.sum()
        if lcur+npos>len(lines):
            return
        #end if
        spos = []
        for i in range(npos):
            spos.append(lines[lcur+i].split())
        #end for
        spos = array(spos)
        pos  = array(spos[:,0:3],dtype=float)
        if cartesian:
            pos = scale*pos
        else:
            pos = dot(pos,axes)
        #end if
        self.set_elem(elem)
        self.pos = pos
        if selective_dynamics or spos.shape[1]>3:
            move = array(spos[:,3:6],dtype=str)
            self.freeze(range(self.size()),directions=move=='F')
        #end if
    #end def read_poscar


    def read_cif(self,filepath,block=None,grammar='1.1',cell='prim'):
        axes,elem,pos,units = read_cif(filepath,block,grammar,cell,args_only=True)
        self.dim = 3
        self.set_axes(axes)
        self.set_elem(elem)
        self.pos = pos
        self.units = units
    #end def read_cif


    def read_fhi_aims(self,filepath):
        if os.path.exists(filepath):
            lines = open(filepath,'r').read().splitlines()
        else:
            lines = filepath.splitlines() # "filepath" is contents
        #end if
        axes = []
        pos  = []
        elem = []
        unit_pos = False
        for line in lines:
            ls = line.strip()
            if len(ls)>0 and ls[0]!='#':
                tokens = ls.split()
                t0 = tokens[0]
                if t0=='lattice_vector':
                    axes.append(tokens[1:])
                elif t0=='atom_frac':
                    pos.append(tokens[1:4])
                    elem.append(tokens[4])
                    unit_pos = True
                elif t0=='atom':
                    pos.append(tokens[1:4])
                    elem.append(tokens[4])
                elif t0.startswith('initial'):
                    None
                else:
                    #None
                    self.error('unrecogonized or not yet supported token in fhi-aims geometry file: {0}'.format(t0))
                #end if
            #end if
        #end for
        axes = array(axes,dtype=float)
        pos  = array(pos,dtype=float)
        if unit_pos:
            pos  = dot(pos,axes)
        #end if
        self.dim = 3
        self.set_axes(axes)
        self.set_elem(elem)
        self.pos   = pos
        self.units = 'A'
    #end def read_fhi_aims


    def write(self,filepath=None,format=None):
        if filepath is None and format is None:
            self.error('please specify either the filepath or format arguments to write()')
        elif format is None:
            if '.' in filepath:
                format = filepath.split('.')[-1]
            else:
                self.error('file format could not be determined\neither request the format directly with the format keyword or add a file format extension to the file name')
            #end if
        #end if
        format = format.lower()
        if format=='xyz':
            c = self.write_xyz(filepath)
        elif format=='xsf':
            c = self.write_xsf(filepath)
        elif format=='fhi-aims':
            c = self.write_fhi_aims(filepath)
        else:
            self.error('file format {0} is unrecognized'.format(format))
        #end if
        return c
    #end def write


    def write_xyz(self,filepath=None,header=True,units='A'):
        if self.dim!=3:
            self.error('write_xyz is currently only implemented for 3 dimensions')
        #end if
        s = self.copy()
        s.change_units(units)
        c=''
        if header:
            c += str(len(s.elem))+'\n\n'
        #end if
        for i in range(len(s.elem)):
            e = s.elem[i]
            p = s.pos[i]
            c+=' {0:2} {1:12.8f} {2:12.8f} {3:12.8f}\n'.format(e,p[0],p[1],p[2])
        #end for
        if filepath!=None:
            open(filepath,'w').write(c)
        #end if
        return c
    #end def write_xyz


    def write_xsf(self,filepath=None):
        if self.dim!=3:
            self.error('write_xsf is currently only implemented for 3 dimensions')
        #end if
        s = self.copy()
        s.change_units('A')
        c  = ' CRYSTAL\n'
        c += ' PRIMVEC\n'
        for a in s.axes:
            c += '   {0:12.8f}  {1:12.8f}  {2:12.8f}\n'.format(*a)
        #end for
        c += ' PRIMCOORD\n'
        c += '   {0} 1\n'.format(len(s.elem))
        for i in range(len(s.elem)):
            e = s.elem[i]
            identified = e in pt.elements
            if not identified:
                if len(e)>2:
                    e = e[0:2]
                elif len(e)==2:
                    e = e[0:1]
                #end if
                identified = e in pt.elements
            #end if
            if not identified:
                self.error('{0} is not an element\nxsf file cannot be written'.format(e))
            #end if
            enum = pt.elements[e].atomic_number
            r = s.pos[i]
            c += '   {0:>3} {1:12.8f}  {2:12.8f}  {3:12.8f}\n'.format(enum,r[0],r[1],r[2])
        #end for
        if filepath!=None:
            open(filepath,'w').write(c)
        #end if
        return c
    #end def write_xsf


    def write_fhi_aims(self,filepath=None):
        s = self.copy()
        s.change_units('A')
        c = ''
        c+='\n'
        for a in s.axes:
            c += 'lattice_vector   {0: 12.8f}  {1: 12.8f}  {2: 12.8f}\n'.format(*a)
        #end for
        c+='\n'
        for p,e in zip(self.pos,self.elem):
            c += 'atom_frac   {0: 12.8f}  {1: 12.8f}  {2: 12.8f}  {3}\n'.format(p[0],p[1],p[2],e)
        #end for
        if filepath!=None:
            open(filepath,'w').write(c)
        #end if
        return c
    #end def write_fhi_aims


    def plot2d_ax(self,ix,iy,*args,**kwargs):
        if self.dim!=3:
            self.error('plot2d_ax is currently only implemented for 3 dimensions')
        #end if
        iz = list(set([0,1,2])-set([ix,iy]))[0]
        ax = self.axes.copy()
        a  = self.axes[iz]
        dc = self.center-ax.sum(0)/2
        pp = array([0*a,ax[ix],ax[ix]+ax[iy],ax[iy],0*a])
        for i in range(len(pp)):
            pp[i]+=dc
            pp[i]-=dot(a,pp[i])/dot(a,a)*a
        #end for
        plot(pp[:,ix],pp[:,iy],*args,**kwargs)
    #end def plot2d_ax


    def plot2d_pos(self,ix,iy,*args,**kwargs):
        if self.dim!=3:
            self.error('plot2d_pos is currently only implemented for 3 dimensions')
        #end if
        iz = list(set([0,1,2])-set([ix,iy]))[0]
        pp = self.pos.copy()
        a = self.axes[iz]
        for i in range(len(pp)):
            pp[i] -= dot(a,pp[i])/dot(a,a)*a
        #end for
        plot(pp[:,ix],pp[:,iy],*args,**kwargs)
    #end def plot2d


    def plot2d(self,pos_style='b.',ax_style='k-'):
        if self.dim!=3:
            self.error('plot2d is currently only implemented for 3 dimensions')
        #end if
        subplot(1,3,1)
        self.plot2d_ax(0,1,ax_style,lw=2)
        self.plot2d_pos(0,1,pos_style)
        title('a1,a2')
        subplot(1,3,2)
        self.plot2d_ax(1,2,ax_style,lw=2)
        self.plot2d_pos(1,2,pos_style)
        title('a2,a3')
        subplot(1,3,3)
        self.plot2d_ax(2,0,ax_style,lw=2)
        self.plot2d_pos(2,0,pos_style)
        title('a3,a1')
    #end def plot2d

    def show(self,viewer='vmd',filepath='/tmp/tmp.xyz'):
        if self.dim!=3:
            self.error('show is currently only implemented for 3 dimensions')
        #end if
        self.write_xyz(filepath)
        os.system(viewer+' '+filepath)
    #end def show


    # minimal ASE Atoms-like interface to Structure objects for spglib
    def get_cell(self):
        return self.axes
    #end def get_cell

    def get_scaled_positions(self):
        return self.pos_unit()
    #end def get_scaled_positions

    def get_number_of_atoms(self):
        return len(self.elem)
    #end def get_number_of_atoms

    def get_atomic_numbers(self):
        an = []
        for e in self.elem:
            an.append(ptable[e].atomic_number)
        #end for
        return array(an,dtype='intc')
    #end def get_atomic_numbers

    def get_magnetic_moments(self):
        self.error('structure objects do not currently support magnetic moments')
    #end def get_magnetic_moments

#end class Structure
Structure.set_operations()


def interpolate_structures(struct1,struct2=None,images=None,min_image=True,recenter=True,match_com=False,repackage=False,chained=False):
    if images is None:
        Structure.class_error('images must be provided','interpolate_structures')
    #end if

    # if a list of structures is provided,
    # interpolate between pairs in the chain of structures
    if isinstance(struct1,(list,tuple)): 
        structures_in = struct1
        structures = []
        for n in xrange(len(structures_in)-1):
            struct1 = structures_in[n]
            struct2 = structures_in[n+1]
            structs = interpolate_structures(struct1,struct2,images,min_image,recenter,match_com,repackage,chained=True)
            if n==0:
                structures.append(structs[0])
            #end if
            structures.extend(structs[1:-1])
            if n==len(structures_in)-2:
                structures.append(structs[-1])
            #end if
        #end for
        return structures
    #end if

    # handle PhysicalSystem objects indirectly
    system1 = None
    system2 = None
    if not isinstance(struct1,Structure):
        system1 = struct1.copy()
        system1.remove_folded()
        struct1 = system1.structure
    #end if
    if not isinstance(struct2,Structure):
        system2 = struct2.copy()
        system2.remove_folded()
        struct2 = system2.structure
    #end if

    # perform the interpolation
    structures = struct1.interpolate(struct2,images,min_image,recenter,match_com)

    # repackage into physical system objects if requested
    if repackage:
        if system1!=None:
            system = system1
        elif system2!=None:
            system = system2
        else:
            Structure.class_error('cannot repackage into physical systems since no system object was provided in place of a structure','interpolate_structures')
        #end if
        systems = []
        for s in structures:
            ps = system.copy()
            ps.structure = s
            systems.append(ps)
        #end for
        result = systems
    else:
        result = structures
    #end if

    return result
#end def interpolate_structures


def structure_animation(filepath,structures,tiling=None):
    path,file = os.path.split(filepath)
    if not file.endswith('xyz'):
        Structure.class_error('only xyz files are supported for now','structure_animation')
    #end if
    anim = ''
    for s in structures:
        if tiling is None:
            anim += s.write_xyz()
        else:
            anim += s.tile(tiling).write_xyz()
        #end if
    #end for
    open(filepath,'w').write(anim)
#end def structure_animation



class DefectStructure(Structure):
    def __init__(self,*args,**kwargs):
        if len(args)>0 and isinstance(args[0],Structure):
            self.transfer_from(args[0],copy=True)
        else:
            Structure.__init__(self,*args,**kwargs)
        #end if
    #end def __init__


    def defect_from_bond_compression(self,compression_cutoff,bond_eq,neighbors):
        bind,bcent,blens = self.bonds(neighbors)
        ind = bind[ abs(blens/bond_eq - 1.) > compression_cutoff ]
        idefect = array(list(set(ind.ravel())))
        defect = self.carve(idefect)
        return defect
    #end def defect_from_bond_compression

    
    def defect_from_displacement(self,displacement_cutoff,reference):
        displacement = self.scalar_displacement(reference)
        idefect = displacement > displacement_cutoff
        defect = self.carve(idefect)
        return defect
    #end def defect_from_displacement


    def compare(self,dist_cutoff,d1,d2=None):
        if d2==None:
            d2 = d1
            d1 = self
        #end if
        res = Sobj()
        natoms1 = len(d1.pos)
        natoms2 = len(d2.pos)
        if natoms1<natoms2:
            dsmall,dlarge = d1,d2
        else:
            dsmall,dlarge = d2,d1
        #end if
        nn = nearest_neighbors(1,dlarge,dsmall)
        dist = dsmall.distances(dlarge[nn.ravel()])
        dmatch = dist<dist_cutoff
        ismall = array(range(len(dsmall.pos)))
        ismall = ismall[dmatch]
        ilarge = nn[ismall]
        if natoms1<natoms2:
            i1,i2 = ismall,ilarge
        else:
            i2,i1 = ismall,ilarge
        #end if
        natoms_match = dmatch.sum()
        res.all_match = natoms1==natoms2 and natoms1==natoms_match
        res.natoms_match = natoms_match
        res.imatch1 = i1
        res.imatch2 = i2
        return res
    #end def compare
#end class DefectStructure



class Crystal(Structure):
    lattice_constants = obj(
        triclinic    = ['a','b','c','alpha','beta','gamma'],
        monoclinic   = ['a','b','c','beta'],
        orthorhombic = ['a','b','c'],
        tetragonal   = ['a','c'],
        hexagonal    = ['a','c'],
        cubic        = ['a'],
        rhombohedral = ['a','alpha']
        )

    lattices = list(lattice_constants.keys())

    centering_types = obj(
        primitive             = 'P',
        base_centered         = ('A','B','C'),
        face_centered         = 'F',
        body_centered         = 'I',
        rhombohedral_centered = 'R'        
        )

    lattice_centerings = obj(
        triclinic = ['P'],
        monoclinic = ['P','A','B','C'],
        orthorhombic = ['P','C','I','F'],
        tetragonal   = ['P','I'],
        hexagonal    = ['P','R'],
        cubic        = ['P','I','F'],
        rhombohedral = ['P']
        )

    centerings = obj(
        P = [],
        A = [[0,.5,.5]],
        B = [[.5,0,.5]],
        C = [[.5,.5,0]],
        F = [[0,.5,.5],[.5,0,.5],[.5,.5,0]],
        I = [[.5,.5,.5]],
        R = [[2./3, 1./3, 1./3],[1./3, 2./3, 2./3]]
        )

    cell_types = set(['primitive','conventional'])

    cell_aliases = obj(
        prim = 'primitive',
        conv = 'conventional'
        )
    cell_classes = obj(
        sc  = 'cubic',
        bcc = 'cubic',
        fcc = 'cubic',
        hex = 'hexagonal'
        )
    for lattice in lattices:
        cell_classes[lattice]=lattice
    #end for

    
    #helpful websites for structures
    #  wikipedia.org
    #  webelements.com
    #  webmineral.com
    #  springermaterials.com


    known_crystals = {
        ('diamond','fcc'):obj(
            lattice   = 'cubic',
            cell      = 'primitive',
            centering = 'F',
            constants = 3.57,
            units     = 'A',
            atoms     = 'C',
            basis     = [[0,0,0],[.25,.25,.25]]
            ),
        ('diamond','sc'):obj(
            lattice   = 'cubic',
            cell      = 'conventional',
            centering = 'F',
            constants = 3.57,
            units     = 'A',
            atoms     = 'C',
            basis     = [[0,0,0],[.25,.25,.25]]
            ),
        ('diamond','prim'):obj(
            lattice   = 'cubic',
            cell      = 'primitive',
            centering = 'F',
            constants = 3.57,
            units     = 'A',
            atoms     = 'C',
            basis     = [[0,0,0],[.25,.25,.25]]
            ),
        ('diamond','conv'):obj(
            lattice   = 'cubic',
            cell      = 'conventional',
            centering = 'F',
            constants = 3.57,
            units     = 'A',
            atoms     = 'C',
            basis     = [[0,0,0],[.25,.25,.25]]
            ),
        ('wurtzite','prim'):obj(
            lattice   = 'hexagonal',
            cell      = 'primitive',
            centering = 'P',
            constants = (3.35,5.22),
            units     = 'A',
            #atoms     = ('Zn','O'),
            #basis     = [[1./3, 2./3, 3./8],[1./3, 2./3, 0]]
            atoms     = ('Zn','O','Zn','O'),
            basis     = [[0,0,5./8],[0,0,0],[2./3,1./3,1./8],[2./3,1./3,1./2]]
            ),
        ('ZnO','prim'):obj(
            lattice   = 'wurtzite',
            cell      = 'prim',
            constants = (3.35,5.22),
            units     = 'A',
            atoms     = ('Zn','O','Zn','O')
            ),
        ('NaCl','prim'):obj(
            lattice   = 'cubic',
            cell      = 'primitive',
            centering = 'F',
            constants = 5.64,
            units     = 'A',
            atoms     = ('Na','Cl'),
            basis     = [[0,0,0],[.5,0,0]],
            basis_vectors = 'conventional'
            ),
        ('rocksalt','prim'):obj(
            lattice   = 'cubic',
            cell      = 'primitive',
            centering = 'F',
            constants = 5.64,
            units     = 'A',
            atoms     = ('Na','Cl'),
            basis     = [[0,0,0],[.5,0,0]],
            basis_vectors = 'conventional'
            ),
        ('copper','prim'):obj(
            lattice   = 'cubic',
            cell      = 'primitive',
            centering = 'F',
            constants = 3.615,
            units     = 'A',
            atoms     = 'Cu'
            ),
        ('calcium','prim'):obj(
            lattice   = 'cubic',
            cell      = 'primitive',
            centering = 'F',
            constants = 5.588,
            units     = 'A',
            atoms     = 'Ca'
            ),
        # http://www.webelements.com/oxygen/crystal_structure.html
        #   Phys Rev 160 694
        ('oxygen','prim'):obj(
            lattice   = 'monoclinic',
            cell      = 'primitive',
            centering = 'C',
            constants = (5.403,3.429,5.086,132.53),
            units     = 'A',
            angular_units = 'degrees',
            atoms     = ('O','O'),
            basis     = [[0,0,1.15/2],[0,0,-1.15/2]],
            basis_vectors = identity(3)
            ),
        # http://en.wikipedia.org/wiki/Calcium_oxide
        # http://www.springermaterials.com/docs/info/10681719_224.html
        ('CaO','prim'):obj(
            lattice   = 'NaCl',
            cell      = 'prim',
            constants = 4.81,
            atoms     = ('Ca','O')
            ),
        ('CaO','conv'):obj(
            lattice   = 'NaCl',
            cell      = 'conv',
            constants = 4.81,
            atoms     = ('Ca','O')
            ),
        # http://en.wikipedia.org/wiki/Copper%28II%29_oxide
        #   http://iopscience.iop.org/0953-8984/3/28/001/
        # http://www.webelements.com/compounds/copper/copper_oxide.html
        # http://www.webmineral.com/data/Tenorite.shtml
        ('CuO','prim'):obj(
            lattice   = 'monoclinic',
            cell      = 'primitive',
            centering = 'C',
            constants = (4.683,3.422,5.128,99.54),
            units     = 'A',
            angular_units = 'degrees',
            atoms     = ('Cu','O','Cu','O'),
            basis     = [[.25,.25,0],[0,.418,.25],
                         [.25,.75,.5],[.5,.5-.418,.75]],
            basis_vectors = 'conventional'
            ),
        ('Ca2CuO3','prim'):obj(# kateryna foyevtsova
            lattice   = 'orthorhombic',
            cell      = 'primitive',
            centering = 'I',
            constants = (3.77,3.25,12.23),
            units     = 'A',
            atoms     = ('Cu','O','O','O','Ca','Ca'),
            basis     = [[   0,   0,   0 ],
                         [ .50,   0,   0 ],
                         [   0,   0, .16026165],
                         [   0,   0, .83973835],
                         [   0,   0, .35077678],
                         [   0,   0, .64922322]],
            basis_vectors = 'conventional'
            ),
        ('La2CuO4','prim'):obj( #tetragonal structure
            lattice   = 'tetragonal',
            cell      = 'primitive',
            centering = 'I',
            constants = (3.809,13.169),
            units     = 'A',
            atoms     = ('Cu','O','O','O','O','La','La'),
            basis     = [[  0,    0,    0],
                         [ .5,    0,    0],
                         [  0,   .5,    0],
                         [  0,    0,  .182],
                         [  0,    0, -.182],
                         [  0,    0,  .362],
                         [  0,    0, -.362]]
            ),
        ('Cl2Ca2CuO2','prim'):obj(
            lattice   = 'tetragonal',
            cell      = 'primitive',
            centering = 'I',
            constants = (3.869,15.05),
            units     = 'A',
            atoms     = ('Cu','O','O','Ca','Ca','Cl','Cl'),
            basis     = [[   0,   0,    0 ],
                         [  .5,   0,    0 ],
                         [   0,  .5,    0 ],
                         [  .5,  .5,  .104],
                         [   0,   0,  .396],
                         [   0,   0,  .183],
                         [  .5,  .5,  .317]],
            basis_vectors = 'conventional'
            ),
        ('Cl2Ca2CuO2','afm'):obj(
            lattice   = 'tetragonal',
            cell      = 'conventional',
            centering = 'P',
            axes      = [[.5,-.5,0],[.5,.5,0],[0,0,1]],
            constants = (2*3.869,15.05),
            units     = 'A',
            atoms     = 4*['Cu','O','O','Ca','Ca','Cl','Cl'],
            basis     = [[   0,   0,    0 ], #Cu
                         [  .25,  0,    0 ],
                         [   0,  .25,   0 ],
                         [  .25, .25, .104],
                         [   0,   0,  .396],
                         [   0,   0,  .183],
                         [  .25, .25, .317],
                         [  .25, .25, .5  ], #Cu
                         [  .5,  .25, .5  ],
                         [  .25, .5,  .5  ],
                         [  .5,  .5,  .604],
                         [  .25, .25, .896],
                         [  .25, .25, .683],
                         [  .5,  .5,  .817],
                         [  .5,   0,    0 ], #Cu2
                         [  .75,  0,    0 ],
                         [  .5,  .25,   0 ],
                         [  .75, .25, .104],
                         [  .5,   0,  .396],
                         [  .5,   0,  .183],
                         [  .75, .25, .317],
                         [  .75, .25, .5  ], #Cu2
                         [   0,  .25, .5  ],
                         [  .75, .5,  .5  ],
                         [   0,  .5,  .604],
                         [  .75, .25, .896],
                         [  .75, .25, .683],
                         [   0,  .5,  .817]],
            basis_vectors = 'conventional'
            ),
        ('CuO2_plane','prim'):obj( 
            lattice   = 'tetragonal',
            cell      = 'primitive',
            centering = 'P',
            constants = (3.809,13.169),
            units     = 'A',
            atoms     = ('Cu','O','O'),
            basis     = [[  0,    0,    0],
                         [ .5,    0,    0],
                         [  0,   .5,    0]]
            ),
        ('graphite_aa','hex'):obj(
            axes      = [[1./2,-sqrt(3.)/2,0],[1./2,sqrt(3.)/2,0],[0,0,1]],
            constants = (2.462,3.525),
            units     = 'A',
            atoms     = ('C','C'),
            basis     = [[0,0,0],[2./3,1./3,0]]
            ),
        ('graphite_ab','hex'):obj(
            axes      = [[1./2,-sqrt(3.)/2,0],[1./2,sqrt(3.)/2,0],[0,0,1]],
            constants = (2.462,3.525),
            units     = 'A',
            cscale    = (1,2),
            atoms     = ('C','C','C','C'),
            basis     = [[0,0,0],[2./3,1./3,0],[0,0,1./2],[1./3,2./3,1./2]]
            ),
        ('graphene','prim'):obj(
            lattice   = 'hexagonal',
            cell      = 'primitive',
            centering = 'P',
            constants = (2.462,15.0),
            units     = 'A',
            atoms     = ('C','C'),
            basis     = [[0,0,0],[2./3,1./3,0]]
            ),
        ('graphene','rect'):obj(
            lattice   = 'orthorhombic',
            cell      = 'conventional',
            centering = 'C',
            constants = (2.462,sqrt(3.)*2.462,15.0),
            units     = 'A',
            atoms     = ('C','C'),
            basis     = [[0,0,0],[1./2,1./6,0]]
            )
        }

    kc_keys = list(known_crystals.keys())
    for (name,cell) in kc_keys:
        desc = known_crystals[name,cell]
        if cell=='prim' and not (name,'conv') in known_crystals:
            cdesc = desc.copy()
            if cdesc.cell=='primitive':
                cdesc.cell = 'conventional'
                known_crystals[name,'conv'] = cdesc
            elif cdesc.cell=='prim':
                cdesc.cell = 'conv'
                known_crystals[name,'conv'] = cdesc
            #end if
        #end if
    #end if
    del kc_keys


    def __init__(self,
                 lattice        = None,
                 cell           = None,
                 centering      = None,
                 constants      = None, 
                 atoms          = None,
                 basis          = None,
                 basis_vectors  = None,
                 tiling         = None,
                 cscale         = None,
                 axes           = None,
                 units          = None,
                 angular_units  = 'degrees',
                 kpoints        = None,
                 kgrid          = None,
                 frozen         = None,
                 magnetization  = None,
                 magnetic_order = None,
                 magnetic_prim  = True,
                 kshift         = (0,0,0),
                 permute        = None,
                 operations     = None,
                 elem           = None, 
                 pos            = None):

        if lattice is None and cell is None and atoms is None and units is None:
            return
        #end if

        gi = obj(
            lattice        = lattice       ,  
            cell           = cell          ,
            centering      = centering     ,
            constants      = constants     ,   
            atoms          = atoms         ,
            basis          = basis         ,
            basis_vectors  = basis_vectors ,
            tiling         = tiling        ,
            cscale         = cscale        ,
            axes           = axes          ,
            units          = units         ,
            angular_units  = angular_units ,
            frozen         = frozen        ,
            magnetization  = magnetization ,
            magnetic_order = magnetic_order,
            magnetic_prim  = magnetic_prim ,
            kpoints        = kpoints       ,
            kgrid          = kgrid         ,
            kshift         = kshift        ,
            permute        = permute       ,
            operations     = operations    ,
            )
        generation_info = gi.copy()

        lattice_in = lattice
        if isinstance(lattice,str):
            lattice=lattice.lower()
        #end if
        if isinstance(cell,str):
            cell=cell.lower()
        #end if

        known_crystal = False
        if (lattice_in,cell) in self.known_crystals:
            known_crystal = True
            lattice_info = self.known_crystals[lattice_in,cell].copy()
        elif (lattice,cell) in self.known_crystals:
            known_crystal = True
            lattice_info = self.known_crystals[lattice,cell].copy()
        #end if
        
        if known_crystal:
            while 'lattice' in lattice_info and 'cell' in lattice_info and (lattice_info.lattice,lattice_info.cell) in self.known_crystals:
                li_old = lattice_info
                lattice_info = self.known_crystals[li_old.lattice,li_old.cell].copy()
                del li_old.lattice
                del li_old.cell
                lattice_info.transfer_from(li_old,copy=False)
            #end while
            if 'cell' in lattice_info:
                cell = lattice_info.cell
            elif cell in self.cell_aliases:
                cell = self.cell_aliases[cell]
            elif cell in self.cell_classes:
                lattice = self.cell_classes[cell]
            else:
                self.error('cell shape '+cell+' is not recognized\n  the variable cell_classes or cell_aliases must be updated to include '+cell)
            #end if
            if 'lattice' in lattice_info:
                lattice = lattice_info.lattice
            #end if
            if 'angular_units' in lattice_info:
                angular_units = lattice_info.angular_units
            #end if
            inputs = obj(
                centering     = centering,
                constants     = constants,
                atoms         = atoms,
                basis         = basis,
                basis_vectors = basis_vectors,
                tiling        = tiling,
                cscale        = cscale,
                axes          = axes,
                units         = units
                )
            for var,val in inputs.iteritems():
                if val is None and var in lattice_info:
                    inputs[var] = lattice_info[var]
                #end if
            #end for
            centering,constants,atoms,basis,basis_vectors,tiling,cscale,axes,units=inputs.list('centering','constants','atoms','basis','basis_vectors','tiling','cscale','axes','units')
        #end if

        if constants is None:
            self.error('the variable constants must be provided')
        #end if
        if atoms is None:
            self.error('the variable atoms must be provided')
        #end if

        if lattice not in self.lattices:
            self.error('lattice type '+str(lattice)+' is not recognized\n  valid lattice types are: '+str(list(self.lattices)))
        #end if
        if cell=='conventional':
            if centering is None:
                self.error('centering must be provided for a conventional cell\n  options for a '+lattice+' lattice are: '+str(self.lattice_centerings[lattice]))
            elif centering not in self.centerings:
                self.error('centering type '+str(centering)+' is not recognized\n  options for a '+lattice+' lattice are: '+str(self.lattice_centerings[lattice]))
            #end if
        #end if
        if isinstance(constants,int) or isinstance(constants,float):
            constants=[constants]
        #end if
        if len(constants)!=len(self.lattice_constants[lattice]):
            self.error('the '+lattice+' lattice depends on the constants '+str(self.lattice_constants[lattice])+'\n you provided '+str(len(constants))+': '+str(constants))
        #end if
        if isinstance(atoms,str):
            if basis!=None:
                atoms = len(basis)*[atoms]
            else:
                atoms=[atoms]
            #end if
        #end if
        if basis is None:
            if len(atoms)==1:
                basis = [(0,0,0)]
            else:
                self.error('must provide as many basis coordinates as basis atoms\n  atoms provided: '+str(atoms)+'\n  basis provided: '+str(basis))
            #end if
        #end if
        if basis_vectors is not None and not isinstance(basis_vectors,str) and len(basis_vectors)!=3:
            self.error('3 basis vectors must be given, you provided '+str(len(basis))+':\n  '+str(basis_vectors))
        #end if

        if tiling is None:
            tiling = (1,1,1)
        #end if
        if cscale is None:
            cscale = len(constants)*[1]
        #end if
        if len(cscale)!=len(constants):
            self.error('cscale and constants must be the same length')
        #end if
        basis  = array(basis)
        tiling = array(tiling,dtype=int)
        cscale = array(cscale)
        constants = cscale*array(constants)

        a,b,c,alpha,beta,gamma = None,None,None,None,None,None
        if angular_units=='radians':
            pi_1o2 = pi/2
            pi_2o3 = 2*pi/3
        elif angular_units=='degrees':
            pi_1o2 = 90.
            pi_2o3 = 120.
        else:
            self.error('angular units must be radians or degrees\n  you provided '+str(angular_units))
        #end if
        if lattice=='triclinic':
            a,b,c,alpha,beta,gamma = constants
        elif lattice=='monoclinic':
            a,b,c,beta = constants
            alpha = gamma = pi_1o2
        elif lattice=='orthorhombic':
            a,b,c = constants
            alpha=beta=gamma=pi_1o2
        elif lattice=='tetragonal':
            a,c = constants
            b=a
            alpha=beta=gamma=pi_1o2
        elif lattice=='hexagonal':
            a,c = constants
            b=a
            alpha=beta=pi_1o2
            gamma=pi_2o3
        elif lattice=='cubic':
            a=constants[0]
            b=c=a
            alpha=beta=gamma=pi_1o2
        elif lattice=='rhombohedral':
            a,alpha = constants
            b=c=a
            beta=gamma=alpha
        #end if
        if angular_units=='degrees':
            alpha *= pi/180
            beta  *= pi/180
            gamma *= pi/180
        #end if

        points = [[0,0,0]]
        #get the conventional axes
        sa,ca = sin(alpha),cos(alpha)
        sb,cb = sin(beta) ,cos(beta)
        sg,cg = sin(gamma),cos(gamma)
        y     = (ca-cg*cb)/sg
        a1c = a*array([1,0,0])
        a2c = b*array([cg,sg,0])
        a3c = c*array([cb,y,sqrt(sb**2-y**2)])
        #a1c = array([a,0,0])
        #a2c = array([b*cos(gamma),b*sin(gamma),0])
        #a3c = array([c*cos(beta),c*cos(alpha)*sin(beta),c*sin(alpha)*sin(beta)])
        axes_conv = array([a1c,a2c,a3c]).copy()


        #from numpy import dot,arccos
        #from numpy.linalg import norm
        #print a,b,c
        #print alpha,beta,gamma
        #print alpha,arccos(dot(a2c,a3c)/(norm(a2c)*norm(a3c)))
        #print beta, arccos(dot(a3c,a1c)/(norm(a3c)*norm(a1c)))
        #print gamma,arccos(dot(a1c,a2c)/(norm(a1c)*norm(a2c)))
        #exit()


        if axes is None:
            if cell not in self.cell_types:
                self.error('cell must be primitive or conventional\n  You provided: '+str(cell))
            #end if
            if cell=='primitive' and centering=='P':
                cell='conventional'
            #end if
            #get the primitive axes
            if centering=='P':
                a1 = a1c
                a2 = a2c
                a3 = a3c
            elif centering=='A':
                a1 = a1c
                a2 = (a2c+a3c)/2
                a3 = (-a2c+a3c)/2
            elif centering=='B':
                a1 = (a1c+a3c)/2
                a2 = a2c
                a3 = (-a1c+a3c)/2
            elif centering=='C':
                a1 = (a1c-a2c)/2
                a2 = (a1c+a2c)/2
                a3 = a3c
            elif centering=='I':
                a1=[ a/2, b/2,-c/2]
                a2=[-a/2, b/2, c/2]
                a3=[ a/2,-b/2, c/2]
            elif centering=='F':
                a1=[a/2, b/2,   0]
                a2=[  0, b/2, c/2]
                a3=[a/2,   0, c/2]
            elif centering=='R':
                a1=[   a,              0,   0]
                a2=[ a/2,   a*sqrt(3.)/2,   0]
                a3=[-a/6, a/(2*sqrt(3.)), c/3]
            else:
                self.error('the variable centering must be specified\n  valid options are: P,A,B,C,I,F,R')
            #end if
            axes_prim = array([a1,a2,a3])
            if cell=='primitive':
                axes = axes_prim
            elif cell=='conventional':            
                axes = axes_conv
                points.extend(self.centerings[centering])
            #end if            
        elif known_crystal:
            axes = dot(diag([a,b,c]),array(axes))
        #end if
        points = array(points,dtype=float)

        elem = []
        pos  = []
        if basis_vectors is None:
            basis_vectors = axes
        elif basis_vectors is 'primitive':
            basis_vectors = axes_prim
        elif basis_vectors is 'conventional':
            basis_vectors = axes_conv
        #end if
        nbasis = len(atoms)
        for point in points:
            for i in range(nbasis):
                atom   = atoms[i]
                bpoint = basis[i]
                p = dot(point,axes) + dot(bpoint,basis_vectors)
                elem.append(atom)
                pos.append(p)
            #end for
        #end for
        pos = array(pos)

        self.set(
            constants = array([a,b,c]),
            angles    = array([alpha,beta,gamma]),
            generation_info = generation_info
            )

        Structure.__init__(
            self,
            axes           = axes,
            scale          = a,
            elem           = elem,
            pos            = pos,
            center         = axes.sum(0)/2,
            units          = units,
            frozen         = frozen,
            magnetization  = magnetization,
            magnetic_order = magnetic_order,
            magnetic_prim  = magnetic_prim,
            tiling         = tiling,
            kpoints        = kpoints,
            kgrid          = kgrid,
            kshift         = kshift,
            permute        = permute,
            rescale        = False,
            operations     = operations)

    #end def __init__
#end class Crystal


class Jellium(Structure):
    prefactors = obj()
    prefactors.transfer_from({1:2*pi,2:4*pi,3:4./3*pi})

    def __init__(self,charge=None,background_charge=None,cell=None,volume=None,density=None,rs=None,dim=3,
                 axes=None,kpoints=None,kweights=None,kgrid=None,kshift=None,units=None,tiling=None):
        del tiling
        if rs!=None:
            if not dim in self.prefactors:
                self.error('only 1,2, or 3 dimensional jellium is currently supported\n  you requested one with dimension {0}'.format(dim))
            #end if
            density = 1.0/(self.prefactors[dim]*rs**dim)
        #end if
        if axes!=None:
            cell = axes
        #end if
        if background_charge!=None:
            charge = background_charge
        #end if
        if cell!=None:
            cell   = array(cell)
            dim    = len(cell)
            volume = det(cell)
        elif volume!=None:
            volume = float(volume)
            cell   = volume**(1./dim)*identity(dim)
        #end if
        if density!=None:
            density = float(density)
            if charge is None and volume!=None:
                charge = density*volume
            elif volume is None and charge!=None:
                volume = charge/density
                cell   = volume**(1./dim)*identity(dim)
            #end if
        #end if
        if charge is None or cell is None:
            self.error('not enough information to form jellium structure\n  information provided:\n  charge: {0}\n  cell: {1}\n  volume: {2}\n  density: {3}\n  rs: {4}\n  dim: {5}'.format(charge,cell,volume,density,rs,dim))
        #end if
        Structure.__init__(self,background_charge=charge,axes=cell,dim=dim,kpoints=kpoints,kweights=kweights,kgrid=kgrid,kshift=kshift,units=units)
    #end def __init__
    
    def density(self):
        return self.background_charge/self.volume()
    #end def density

    def rs(self):
        return 1.0/(self.density()*self.prefactors[self.dim])**(1./self.dim)
    #end def rs

    def tile(self):
        self.not_implemented()
    #end def tile
#end class Jellium


    
    




def generate_cell(shape,tiling=None,scale=1.,units=None,struct_type=Structure):
    if tiling is None:
        tiling = (1,1,1)
    #end if
    axes = Sobj()
    axes.sc  =  1.*array([[ 1,0,0],[0, 1,0],[0,0, 1]])
    axes.bcc = .5*array([[-1,1,1],[1,-1,1],[1,1,-1]])
    axes.fcc = .5*array([[ 1,1,0],[1, 0,1],[0,1, 1]])
    ax     = dot(diag(tiling),axes[shape])
    center = ax.sum(0)/2
    c = Structure(axes=ax,scale=scale,center=center,units=units)
    if struct_type!=Structure:
        c=c.upcast(struct_type)
    #end if
    return c
#end def generate_cell



def generate_structure(type='crystal',*args,**kwargs):
    if type=='crystal':
        s = generate_crystal_structure(*args,**kwargs)
    elif type=='defect':
        s = generate_defect_structure(*args,**kwargs)
    elif type=='atom':
        s = generate_atom_structure(*args,**kwargs)
    elif type=='dimer':
        s = generate_dimer_structure(*args,**kwargs)
    elif type=='trimer':
        s = generate_trimer_structure(*args,**kwargs)
    elif type=='jellium':
        s = generate_jellium_structure(*args,**kwargs)
    elif type=='empty':
        s = Structure()
    elif type=='basic':
        s = Structure(*args,**kwargs)
    else:
        Structure.class_error(str(type)+' is not a valid structure type\noptions are crystal, defect, atom, dimer, trimer, jellium, empty, or basic')
    #end if
    return s
#end def generate_structure




def generate_atom_structure(atom=None,units='A',Lbox=None,skew=0,axes=None,kgrid=(1,1,1),kshift=(0,0,0),bconds=tuple('nnn'),struct_type=Structure):
    if atom is None:
        Structure.class_error('atom must be provided','generate_atom_structure')
    #end if
    if Lbox!=None:
        axes = [[Lbox*(1-skew),0,0],[0,Lbox,0],[0,0,Lbox*(1+skew)]]
    #end if
    if axes is None:
        s = Structure(elem=[atom],pos=[[0,0,0]],units=units,bconds=bconds)
    else:
        s = Structure(elem=[atom],pos=[[0,0,0]],axes=axes,kgrid=kgrid,kshift=kshift,bconds=bconds,units=units)
        s.center_molecule()
    #end if

    return s
#end def generate_atom_structure


def generate_dimer_structure(dimer=None,units='A',separation=None,Lbox=None,skew=0,axes=None,kgrid=(1,1,1),kshift=(0,0,0),bconds=tuple('nnn'),struct_type=Structure,axis='x'):
    if dimer is None:
        Structure.class_error('dimer atoms must be provided to construct dimer','generate_dimer_structure')
    #end if
    if separation is None:
        Structure.class_error('separation must be provided to construct dimer','generate_dimer_structure')
    #end if
    if Lbox!=None:
        axes = [[Lbox*(1-skew),0,0],[0,Lbox,0],[0,0,Lbox*(1+skew)]]
    #end if
    if axis=='x':
        p2 = [separation,0,0]
    elif axis=='y':
        p2 = [0,separation,0]
    elif axis=='z':
        p2 = [0,0,separation]
    else:
        Structure.class_error('dimer orientation axis must be x,y,z\n  you provided: {0}'.format(axis),'generate_dimer_structure')
    #end if
    if axes is None:
        s = Structure(elem=dimer,pos=[[0,0,0],p2],units=units,bconds=bconds)
    else:
        s = Structure(elem=dimer,pos=[[0,0,0],p2],axes=axes,kgrid=kgrid,kshift=kshift,units=units,bconds=bconds)
        s.center_molecule()
    #end if
    return s
#end def generate_dimer_structure


def generate_trimer_structure(trimer=None,units='A',separation=None,angle=None,Lbox=None,skew=0,axes=None,kgrid=(1,1,1),kshift=(0,0,0),struct_type=Structure,axis='x',axis2='y',angular_units='degrees',plane_rot=None):
    if trimer is None:
        Structure.class_error('trimer atoms must be provided to construct trimer','generate_trimer_structure')
    #end if
    if separation is None:
        Structure.class_error('separation must be provided to construct trimer','generate_trimer_structure')
    #end if
    if len(separation)!=2:
        Structure.class_error('two separation distances (atom1-atom2,atom1-atom3) must be provided to construct trimer\nyou provided {0} separation distances'.format(len(separation)),'generate_trimer_structure')
    #end if
    if angle is None:
        Structure.class_error('angle must be provided to construct trimer','generate_trimer_structure')
    #end if
    if angular_units=='degrees':
        angle *= pi/180
    elif not angular_units.startswith('rad'):
        Structure.class_error('angular units must be degrees or radians\nyou provided: {0}'.format(angular_units),'generate_trimer_structure')
    #end if
    if axis==axis2:
        Structure.class_error('axis and axis2 must be different to define the trimer plane\nyou provided {0} for both'.format(axis),'generate_trimer_structure')
    #end if
    if Lbox!=None:
        axes = [[Lbox*(1-skew),0,0],[0,Lbox,0],[0,0,Lbox*(1+skew)]]
    #end if
    p1 = [0,0,0]
    if axis=='x':
        p2 = [separation[0],0,0]
    elif axis=='y':
        p2 = [0,separation[0],0]
    elif axis=='z':
        p2 = [0,0,separation[0]]
    else:
        Structure.class_error('trimer bond1 (atom2-atom1) orientation axis must be x,y,z\n  you provided: {0}'.format(axis),'generate_trimer_structure')
    #end if
    r = separation[1]
    c = cos(angle)
    s = sin(angle)
    axpair = axis+axis2
    if axpair=='xy':
        p3 = [r*c,r*s,0]
    elif axpair=='yx':
        p3 = [r*s,r*c,0]
    elif axpair=='yz':
        p3 = [0,r*c,r*s]
    elif axpair=='zy':
        p3 = [0,r*s,r*c]
    elif axpair=='zx':
        p3 = [r*s,0,r*c]
    elif axpair=='xz':
        p3 = [r*c,0,r*s]
    else:
        Structure.class_error('trimer bond2 (atom3-atom1) orientation axis must be x,y,z\n  you provided: {0}'.format(axis2),'generate_trimer_structure')
    #end if
    if axes is None:
        s = Structure(elem=trimer,pos=[p1,p2,p3],units=units)
    else:
        s = Structure(elem=trimer,pos=[p1,p2,p3],axes=axes,kgrid=kgrid,kshift=kshift,units=units)
        s.center_molecule()
    #end if
    if plane_rot!=None:
        s.rotate_plane(axpair,plane_rot,angular_units)
    #end if
    return s
#end def generate_trimer_structure


def generate_jellium_structure(*args,**kwargs):
    return Jellium(*args,**kwargs)
#end def generate_jellium_structure




def generate_crystal_structure(lattice=None,cell=None,centering=None,
                               constants=None,atoms=None,basis=None,
                               basis_vectors=None,tiling=None,cscale=None,
                               axes=None,units=None,angular_units='degrees',
                               magnetization=None,magnetic_order=None,magnetic_prim=True,
                               kpoints=None,kweights=None,kgrid=None,kshift=(0,0,0),permute=None,
                               structure=None,shape=None,element=None,scale=None, #legacy inputs
                               operations=None,
                               struct_type=Crystal,elem=None,pos=None,frozen=None,
                               posu=None):    

    if structure!=None:
        lattice = structure
    #end if
    if shape!=None:
        cell = shape
    #end if
    if element!=None:
        atoms = element
    #end if
    if scale!=None:
        constants = scale
    #end if

    #interface for total manual specification
    # this is only here because 'crystal' is default and must handle other cases
    if elem is not None and (pos is not None or posu is not None):  
        return Structure(
            axes           = axes,
            elem           = elem,
            pos            = pos,
            units          = units,
            frozen         = frozen,
            magnetization  = magnetization,
            magnetic_order = magnetic_order,
            magnetic_prim  = magnetic_prim,
            tiling         = tiling,
            kpoints        = kpoints,
            kgrid          = kgrid,
            kshift         = kshift,
            permute        = permute,
            rescale        = False,
            operations     = operations,
            posu           = posu)
    elif isinstance(structure,Structure):
        if tiling!=None:
            structure = structure.tile(tiling)
        #end if
        if kpoints!=None:
            structure.add_kpoints(kpoints,kweights)
        #end if
        if kgrid!=None:
            structure.add_kmesh(kgrid,kshift)
        #end if        
        return structure
    #end if


    s=Crystal(
        lattice        = lattice       ,  
        cell           = cell          ,
        centering      = centering     ,
        constants      = constants     ,   
        atoms          = atoms         ,
        basis          = basis         ,
        basis_vectors  = basis_vectors ,
        tiling         = tiling        ,
        cscale         = cscale        ,
        axes           = axes          ,
        units          = units         ,
        angular_units  = angular_units ,
        frozen         = frozen        ,
        magnetization  = magnetization ,
        magnetic_order = magnetic_order,
        magnetic_prim  = magnetic_prim ,
        kpoints        = kpoints       ,
        kgrid          = kgrid         ,
        kshift         = kshift        ,
        permute        = permute       ,
        operations     = operations    ,
        elem           = elem          ,
        pos            = pos
        )

    if struct_type!=Crystal:
        s=s.upcast(struct_type)
    #end if

    return s
#end def generate_crystal_structure



defects = obj(
    diamond = obj(
        H = obj(
            pristine = [[0,0,0]],
            defect   = [[0,0,0],[.625,.375,.375]]
            ),
        T = obj(
            pristine = [[0,0,0]],
            defect   = [[0,0,0],[.5,.5,.5]]
            ),
        X = obj(
            pristine = [[.25,.25,.25]],
            defect   = [[.39,.11,.15],[.11,.39,.15]]
            ),
        FFC = obj(
            #pristine = [[   0,   0,   0],[.25,.25,.25]],
            #defect   = [[.151,.151,-.08],[.10,.10,.33]]
            pristine = [[   0,   0,    0],[.25 ,.25 ,.25 ],[.5  ,.5  ,    0],[.75 ,.75 ,.25 ],[1.5  ,1.5  ,0   ],[1.75 ,1.75 ,.25 ]],
            defect   = [[.151,.151,-.081],[.099,.099,.331],[.473,.473,-.059],[.722,.722,.230],[1.528,1.528,.020],[1.777,1.777,.309]]
            )
        )
    )


def generate_defect_structure(defect,structure,shape=None,element=None,
                              tiling=None,scale=1.,kgrid=None,kshift=(0,0,0),
                              units=None,struct_type=DefectStructure):
    if structure in defects:
        dstruct = defects[structure]
    else:
        DefectStructure.class_error('defects for '+structure+' structure have not yet been implemented')
    #end if
    if defect in dstruct:
        drep = dstruct[defect]
    else:
        DefectStructure.class_error(defect+' defect not found for '+structure+' structure')
    #end if

    ds = generate_crystal_structure(
        structure = structure,
        shape     = shape,
        element   = element,
        tiling    = tiling,
        scale     = 1.0,
        kgrid     = kgrid,
        kshift    = kshift,
        units     = units,
        struct_type = struct_type
        )

    ds.replace(drep.pristine,pos=drep.defect)

    ds.rescale(scale)
    
    return ds
#end def generate_defect_structure


def read_structure(filepath,elem=None,format=None):
    s = generate_structure('empty')
    s.read(filepath,elem=elem,format=format)
    return s
#end def read_structure




#if __name__=='__main__':
#    from numpy.random import rand
#    from matplotlib.pyplot import figure,plot,show
#
#    ax = array([[1.0,.3,.1],[.2,1.2,-.1],[.2,.1,1.]])
#    #ax = array([[1.0,0,0],[0,1.,0],[0,0,1.]])
#    pos = 4*(rand(50,3)-.5)
#    c = (ax[0]+ax[1])/2
#    elem = []
#    for i in range(len(pos)):
#        elem.append('Ge')
#    #end for
#    s = Structure(axes=ax,pos=pos,elem=elem)
#
#    #figure()
#    #plot(s.pos[:,0],s.pos[:,1],'bo')
#    #plot([0,x1,x1+x2,x2,0],[0,y1,y1+y2,y2,0],'k-',lw=2)
#    #s.recenter(c)
#    #plot(s.pos[:,0],s.pos[:,1],'r.')
#
#    #figure()
#    #s.plot2d('bo')
#    #s.recenter(c)
#    #s.plot2d('r.')
#    #show()
#
#    figure()
#    s.recenter(c)
#    s.plot2d('bo')
#    cs=s.carve(s.axes/2,s.center)
#    cs.plot2d('r.')
#    show()
##end if



if __name__=='__main__':
    large = generate_structure(
        type      = 'crystal',
        structure = 'diamond',
        cell      = 'fcc',
        atoms     = 'Ge',
        constants = 5.639,
        units     = 'A',
        tiling    = (2,2,2),
        kgrid     = (1,1,1),
        kshift    = (0,0,0),
        )
    
    small = large.folded_structure

    print small.kpoints_unit()
#end if
