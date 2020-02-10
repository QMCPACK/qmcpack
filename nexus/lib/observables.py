

# Python standard library imports
import os
import inspect
from time import clock
from copy import deepcopy

# Non-standard Python imports
import numpy as np
import matplotlib.pyplot as plt
import h5py

# Nexus imports
import memory
from unit_converter import convert
from generic import obj
from developer import DevBase,log,error,ci
from numerics import simstats
from grid_functions import grid_function,read_grid,StructuredGrid
from grid_functions import SpheroidGrid
from structure import Structure,read_structure
from fileio import XsfFile




class VLog(DevBase):

    verbosity_levels = obj(
        none = 0,
        low  = 1,
        high = 2,
        )

    def __init__(self):
        self.tstart    = clock()
        self.tlast     = self.tstart
        self.mstart    = memory.resident(children=True)
        self.mlast     = self.mstart
        self.verbosity = self.verbosity_levels.low
        self.indent    = 0
    #end def __init__


    def __call__(self,msg,level='low',n=0,time=False,mem=False,width=75):
        if self.verbosity==self.verbosity_levels.none:
            return
        elif self.verbosity >= self.verbosity_levels[level]:
            if mem or time:
                npad = max(0,width-2*(n+self.indent)-len(msg)-36)
                if npad>0:
                    msg += npad*' '
                #end if
                if mem:
                    dm = 1e6 # MB
                    mnow = memory.resident(children=True)
                    msg += '  (mem add {:6.2f}, tot {:6.2f})'.format((mnow-self.mlast)/dm,(mnow-self.mstart)/dm)
                    self.mlast = mnow
                #end if
                if time:
                    tnow = clock()
                    msg += '  (t elap {:7.3f}, tot {:7.3f})'.format(tnow-self.tlast,tnow-self.tstart)
                    self.tlast = tnow
                #end if
            #end if
            log(msg,n=n+self.indent)
        #end if
    #end def __init__

    def increment(self,n=1):
        self.indent += n
    #end def increment

    def decrement(self,n=1):
        self.indent -= n
    #end def decrement

    def set_none(self):
        self.verbosity = self.verbosity_levels.none
    #end def set_none

    def set_low(self):
        self.verbosity = self.verbosity_levels.low
    #end def set_low

    def set_high(self):
        self.verbosity = self.verbosity_levels.high
    #end def set_high

    def set_verbosity(level):
        if level not in self.verbosity_levels:
            vlinv = self.verbosity_levels.inverse()
            error('Cannot set verbosity level to "{}".\nValid options are: {}'.format(level,[vlinv[i] for i in sorted(vlinv.keys())]))
        #end if
        self.verbosity = self.verbosity_levels[level]
    #end def set_verbosity
#end class VLog
vlog = VLog()



def set_verbosity(level):
    vlog.set_verbosity(level)
#end def set_verbosity


class AttributeProperties(DevBase):
    def __init__(self,**kwargs):
        self.assigned   = set(kwargs.keys())
        self.name       = kwargs.pop('name'      , None )
        self.dest       = kwargs.pop('dest'      , None )
        self.type       = kwargs.pop('type'      , None )
        self.default    = kwargs.pop('default'   , None )
        self.no_default = kwargs.pop('no_default', False)
        self.deepcopy   = kwargs.pop('deepcopy'  , False)
        self.required   = kwargs.pop('required'  , False)
        if len(kwargs)>0:
            self.error('Invalid init variable attributes received.\nInvalid attributes:\n{}\nThis is a developer error.'.format(obj(kwargs)))
        #end if
    #end def __init__
#end class AttributeProperties



class DefinedAttributeBase(DevBase):

    @classmethod
    def set_unassigned_default(cls,default):
        cls.unassigned_default = default
    #end def set_unassigned_default

    @classmethod
    def define_attributes(cls,*other_cls,**attribute_properties):
        if len(other_cls)==1 and issubclass(other_cls[0],DefinedAttributeBase):
            cls.obtain_attributes(other_cls[0])
        #end if
        if cls.class_has('attribute_definitions'):
            attr_defs = cls.attribute_definitions
        else:
            attr_defs = obj()
            cls.class_set(
                attribute_definitions = attr_defs
                )
        #end if
        for name,attr_props in attribute_properties.items():
            attr_props = AttributeProperties(**attr_props)
            attr_props.name = name
            if name not in attr_defs:
                attr_defs[name] = attr_props
            else:
                p = attr_defs[name]
                for n in attr_props.assigned:
                    p[n] = attr_props[n]
                #end for
            #end if
        #end for
        if cls.class_has('unassigned_default'):
            for p in attr_defs:
                if 'default' not in p.assigned:
                    p.default = cls.unassigned_default
                #end if
            #end for
        #end if
        required_attributes = set()
        deepcopy_attributes = set()
        typed_attributes    = set()
        toplevel_attributes = set()
        sublevel_attributes = set()
        for name,props in attr_defs.items():
            if props.required:
                required_attributes.add(name)
            #end if
            if props.deepcopy:
                deepcopy_attributes.add(name)
            #end if
            if props.type is not None:
                typed_attributes.add(name)
            #end if
            if props.dest is None:
                toplevel_attributes.add(name)
            else:
                sublevel_attributes.add(name)
            #end if
        #end for
        cls.class_set(
            required_attributes = required_attributes,
            deepcopy_attributes = deepcopy_attributes,
            typed_attributes    = typed_attributes,
            toplevel_attributes = toplevel_attributes,
            sublevel_attributes = sublevel_attributes,
            )
    #end def define_attributes


    @classmethod
    def obtain_attributes(cls,super_cls):
        cls.class_set(
            attribute_definitions = super_cls.attribute_definitions.copy()
            )
    #end def obtain_attributes


    def __init__(self,**values):
        if len(values)>0:
            self.set_default_attributes()
            self.set_attributes(**values)
        #end if
    #end def __init__


    def initialize(self,**values):
        self.set_default_attributes()
        if len(values)>0:
            self.set_attributes(**values)
        #end if
    #end def initialize


    def set_default_attributes(self):
        cls = self.__class__
        props = cls.attribute_definitions
        for name in cls.toplevel_attributes:
            self._set_default_attribute(name,props[name])
        #end for
        for name in cls.sublevel_attributes:
            self._set_default_attribute(name,props[name])
        #end for
    #end def set_default_attributes


    def set_attributes(self,**values):
        cls = self.__class__
        value_names = set(values.keys())
        attr_names  = set(cls.attribute_definitions.keys())
        invalid     = value_names - attr_names
        if len(invalid)>0:
            v = obj()
            v.transfer_from(values,invalid)
            self.error('Attempted to set unrecognized attributes\nUnrecognized attributes:\n{}'.format(v))
        #end if
        missing = set(cls.required_attributes) - value_names
        if len(missing)>0:
            msg = ''
            for n in sorted(missing):
                msg += '\n  '+n
            #end for
            self.error('Required attributes are missing.\nPlease provide the following attributes during initialization:{}'.format(msg))
        #end if
        props = cls.attribute_definitions
        toplevel_names = value_names & cls.toplevel_attributes
        for name in toplevel_names:
            self._set_attribute(self,name,values[name],props[name])
        #end for
        sublevel_names = value_names - toplevel_names
        for name in sublevel_names:
            p = props[name]
            if p.dest not in self:
                self.error('Attribute destination "{}" does not exist at the top level.\nThis is a developer error.'.format(p.dest))
            #end if
            self._set_attribute(self[p.dest],name,values[name],p)
        #end for
    #end def set_attributes


    def check_attributes(self,exit=False):
        msg = ''
        cls = self.__class__
        a = obj()
        for name in cls.toplevel_attributes:
            if name in self:
                a[name] = self[name]
            #end if
        #end for
        props = cls.attribute_definitions
        for name in cls.sublevel_attributes:
            p = props[name]
            if p.dest in self:
                sub = self[p.dest]
                if name in sub:
                    a[name] = sub[name]
                #end if
            #end if
        #end for
        present = set(a.keys())
        missing = cls.required_attributes - present
        if len(missing)>0:
            m = ''
            for n in sorted(missing):
                m += '\n  '+n
            #end for
            msg += 'Required attributes are missing.\nPlease provide the following attributes during initialization:{}\n'.format(m)
        #end if
        for name in cls.typed_attributes:
            if name in a:
                p = props[name]
                v = a[name]
                if not isinstance(v,p.type):
                    msg += 'Attribute "{}" has invalid type.\n  Type expected: {}\n  Type present: {}\n'.format(name,p.type.__name__,v.__class__.__name__)
                #end if
            #end if
        #end for
        valid = len(msg)==0
        if not valid and exit:
            self.error(msg)
        #end if
        return valid
    #end def check_attributes


    def check_unassigned(self,value):
        cls = self.__class__
        unassigned = cls.class_has('unassigned_default') and value==cls.unassigned_default
        return unassigned
    #end def check_unassigned


    def set_attribute(self,name,value):
        cls = self.__class__
        props = cls.attribute_definitions
        if name not in props:
            self.error('Cannot set unrecognized attribute "{}".\nValid options are: {}'.format(name,sorted(props.keys())))
        #end if
        p = props[name]
        if p.type is not None and not isinstance(value,p.type):
            self.error('Cannot set attribute "{}".\nExpected value with type: {}\nReceived value with type: {}'.format(name,p.type.__name__,value.__class__.__name__))
        #end if
        if p.deepcopy:
            value = deepcopy(value)
        #end if
        if p.dest is None:
            self[name] = value
        elif p.dest not in self:
            self.error('Cannot set attribute "{}".\nAttribute destination "{}" does not exist.'.format(name,p.dest))
        else:
            self[p.dest][name] = value
        #end if
    #end def set_attribute


    def get_attribute(self,name,present=False,assigned=False):
        require_present  = present
        require_assigned = assigned
        require_present  |= require_assigned
        cls = self.__class__
        props = cls.attribute_definitions
        if name not in props:
            self.error('Cannot get unrecognized attribute "{}".\nValid options are: {}'.format(name,sorted(props.keys())))
        #end if
        p = props[name]
        present = False
        value   = None
        if p.dest is None:
            if name in self:
                present = True
                value   = self[name]
            #end if
        elif p.dest in self and name in self[p.dest]:
            present = True
            value   = self[p.dest][name]
        #end if
        if require_present:
            unassigned = self.check_unassigned(value)
            if not present or (unassigned and require_assigned):
                extra = ''
                if p.dest is not None:
                    extra = ' at location "{}"'.format(p.dest)
                #end if
                if not present:
                    msg = 'Cannot get attribute "{}"{}.\nAttribute does not exist.'.format(name,extra)
                else:
                    msg = 'Cannot get attribute "{}"{}.\nAttribute has not been assigned.'.format(name,extra)
                #end if
                self.error(msg)
            #end if
        #end if
        return value
    #end def get_attribute


    def has_attribute(self,name):
        return not (name not in self or self.check_unassigned(self[name]))
    #end def has_attribute


    def _set_default_attribute(self,name,props):
        p = props
        if p.no_default:
            return
        #end if
        value = p.default
        if inspect.isclass(value) or inspect.isfunction(value):
            value = value()
        #end if
        if p.dest is None:
            self[name] = value
        elif p.dest not in self:
            self.error('Attribute destination "{}" does not exist at the top level.\nThis is a developer error.'.format(p.dest))
        else:
            self[p.dest][name] = value
        #end if
    #end def _set_default_attribute


    def _set_attribute(self,container,name,value,props):
        p = props
        if p.type is not None and not isinstance(value,p.type):
            self.error('Cannot set attribute "{}".\nExpected value with type: {}\nReceived value with type: {}'.format(name,p.type.__name__,value.__class__.__name__))
        #end if
        if p.deepcopy:
            value = deepcopy(value)
        #end if
        container[name] = value
    #end def _set_attribute

#end class DefinedAttributeBase



class Observable(DefinedAttributeBase):
    def __init__(self,**values):
        self.initialize(**values)
    #end def __init__

    def initialize(self,**values):
        DefinedAttributeBase.initialize(self,**values)
        if len(values)>0:
            self.info.initialized = True
        #end if
    #end def initialize
#end class Observable

Observable.set_unassigned_default(None)

Observable.define_attributes(
    info = obj( 
        type    = obj, 
        default = obj,
        ),
    initialized = obj(
        dest    = 'info',
        type    = bool,
        default = False,
        ),
    structure = obj(
        dest     = 'info',
        type     = Structure,
        default  = None,
        deepcopy = True,
        ),
    )



class ObservableWithComponents(Observable):

    component_names   = None
    default_component = None


    def process_component_name(self,name,allow_all=False):
        is_all = name=='all'
        if name is None:
            name = self.default_component
        #end if
        if not allow_all and is_all:
            self.error('A specific component must be requested.\nReceived request for "all".\nValid options are: {}'.format(self.compnent_names))
        #end if
        if not allow_all:
            return name
        else:
            return name,is_all
        #end if
    #end def process_component_name


    def component(self,name=None):
        name = self.process_component_name(name)
        comp = self.get_attribute(name,assigned=True)
        return comp
    #end def component


    def components(self):
        comps = obj()
        for c in self.component_names:
            if c in self:
                comps[c] = self[c]
            #end if
        #end for
        return comps
    #end def components


    def operate_over_components(func,*args,**kwargs):
        retvals = obj()
        for c in self.component_names:
            if c in self:
                kwargs['component'] = c
                retvals[c] = func(self,*args,**kwargs)
            #end if
        #end for
        return retvals
    #end def operate_over_components

#end class ObservableWithComponents




def read_eshdf_nofk_data(filename,Ef):
    from numpy import array,pi,dot,sqrt,abs,zeros
    from numpy.linalg import inv,det
    from hdfreader import read_hdf

    def h5int(i):
        return array(i,dtype=int)[0]
    #end def h5int

    # Use slightly shifted Fermi energy
    E_fermi  = Ef + 1e-8

    # Open the HDF file w/o loading the arrays into memory (view mode)
    vlog('Reading '+filename)
    h        = read_hdf(filename,view=True)

    # Get the G-vectors in cell coordinates
    gvu      = array(h.electrons.kpoint_0.gvectors)

    # Get the untiled cell axes
    axes     = array(h.supercell.primitive_vectors)

    # Compute the k-space cell axes
    kaxes    = 2*pi*inv(axes).T

    # Convert G-vectors from cell coordinates to atomic units 
    gv       = dot(gvu,kaxes)

    # Get number of kpoints/twists, spins, and G-vectors
    nkpoints = h5int(h.electrons.number_of_kpoints)
    nspins   = h5int(h.electrons.number_of_spins)
    ngvecs   = len(gv)

    # Process the orbital data
    data     = obj()
    for k in range(nkpoints):
        vlog('Processing k-point {:>3}'.format(k),n=1,time=True)
        kin_k   = obj()
        eig_k   = obj()
        k_k     = obj()
        nk_k    = obj()
        nelec_k = zeros((nspins,),dtype=float)
        kp      = h.electrons['kpoint_'+str(k)]
        gvs     = dot(array(kp.reduced_k),kaxes)
        gvk     = gv.copy()
        for d in range(3):
            gvk[:,d] += gvs[d]
        #end for
        kinetic=(gvk**2).sum(1)/2 # Hartree units
        for s in range(nspins):
            kin_s   = []
            eig_s   = []
            k_s     = gvk
            nk_s    = zeros((ngvecs,),dtype=float)
            nelec_s = 0
            path    = 'electrons/kpoint_{0}/spin_{1}'.format(k,s)
            spin    = h.get_path(path)
            eigs    = convert(array(spin.eigenvalues),'Ha','eV')
            nstates = h5int(spin.number_of_states)
            for st in range(nstates):
                eig = eigs[st]
                if eig<E_fermi:
                    stpath   = path+'/state_{0}/psi_g'.format(st)
                    psi      = array(h.get_path(stpath))
                    nk_orb   = (psi**2).sum(1)
                    kin_orb  = (kinetic*nk_orb).sum()
                    nelec_s += nk_orb.sum()
                    nk_s    += nk_orb
                    kin_s.append(kin_orb)
                    eig_s.append(eig)
                #end if
            #end for
            data[k,s] = obj(
                kpoint = array(kp.reduced_k),
                kin    = array(kin_s),
                eig    = array(eig_s),
                k      = k_s,
                nk     = nk_s,
                ne     = nelec_s,
                )
        #end for
    #end for
    res = obj(
        orbfile  = filename,
        E_fermi  = E_fermi,
        axes     = axes,
        kaxes    = kaxes,
        nkpoints = nkpoints,
        nspins   = nspins,
        data     = data,
        )

    return res
#end def read_eshdf_nofk_data



class MomentumDistribution(ObservableWithComponents):
    component_names = ('tot','pol','u','d')
    
    default_component = 'tot'

    def get_raw_data(self):
        data = self.get_attribute('raw',assigned=True)
        if len(data)==0:
            self.error('Raw n(k) data is not present.')
        #end if
        return data
    #end def get_raw_data


    def filter_raw_data(self,filter_tol=1e-5,store=True):
        vlog('Filtering raw n(k) data with tolerance {:6.4e}'.format(filter_tol))
        prior_tol = self.get_attribute('raw_filter_tol',present=True)
        data  = self.get_raw_data()
        if prior_tol is not None and prior_tol<=filter_tol:
            vlog('Filtering applied previously with tolerance {:6.4e}, skipping.'.format(prior_tol))
            return data
        #end if
        k     = data.first().k
        km    = np.linalg.norm(k,axis=1)
        kmax  = 0.
        order = km.argsort()
        for s,sdata in data.items():
            vlog('Finding kmax for {} data'.format(s),n=1,time=True)
            nk = sdata.nk
            for n in reversed(order):
                if nk[n]>filter_tol:
                    break
                #end if
            #end for
            kmax = max(km[n],kmax)
        #end for
        vlog('Original kmax: {:8.4f}'.format(km.max()),n=2)
        vlog('Filtered kmax: {:8.4f}'.format(kmax),n=2)
        vlog('Applying kmax filter to data',n=1,time=True)
        keep = km<kmax
        k = k[keep]
        vlog('size before filter: {}'.format(len(keep)),n=2)
        vlog('size  after filter: {}'.format(len(k)),n=2)
        vlog('fraction: {:6.4e}'.format(len(k)/len(keep)),n=2)
        if store:
            new_data = data
            self.set_attribute('raw_filter_tol',filter_tol)
        else:
            new_data = obj()
        #end if
        for s in data.keys():
            if s not in new_data:
                new_data[s] = obj()
            #end if
            sdata    = new_data[s]
            sdata.k  = k
            sdata.nk = data[s].nk[keep]
        #end for
        if store:
            vlog('Overwriting original raw n(k) with filtered data',n=1)
        #end if
        vlog('Filtering complete',n=1,time=True)
        return new_data
    #end def filter_raw_data


    def map_raw_data_onto_grid(self,unfold=False,filter_tol=1e-5):
        vlog('\nMapping raw n(k) data onto regular grid')
        data = self.get_raw_data()
        structure = self.get_attribute('structure',assigned=unfold)
        if structure is not None:
            kaxes = structure.kaxes
        else:
            kaxes = self.get_attribute('kaxes',required=True)
        #end if
        if filter_tol is not None:
            vlog.increment()
            data = self.filter_raw_data(filter_tol,store=False)
            vlog.decrement()
        #end if
        if not unfold:
            for s,sdata in data.items():
                vlog('Mapping {} data onto grid'.format(s),n=1,time=True)
                self[s] = grid_function(
                    points = k,
                    values = nk,
                    axes   = kaxes,
                    )
            #end for
        else:
            rotations = structure.point_group_operations()
            for s,sdata in data.items():
                if s=='d' and 'u' in data and id(sdata)==id(data.u):
                    continue
                #end if
                vlog('Unfolding {} data'.format(s),n=1,time=True)
                k   = []
                nk  = []
                ks  = sdata.k
                nks = sdata.nk
                for n,R in enumerate(rotations):
                    vlog('Processing rotation {:<3}'.format(n),n=2,mem=True)
                    k.extend(np.dot(ks,R))
                    nk.extend(nks)
                #end for
                k  = np.array(k ,dtype=float)
                nk = np.array(nk,dtype=float)
                vlog('Unfolding finished',n=2,time=True)

                vlog('Mapping {} data onto grid'.format(s),n=1,time=True)
                vlog.increment(2)
                self[s] = grid_function(
                    points  = k,
                    values  = nk,
                    axes    = kaxes,
                    average = True,
                    )
                vlog.decrement(2)
            #end for
        #end if
        if 'd' not in self and 'u' in self:
            self.d = self.u
        #end if
        vlog('Mapping complete',n=1,time=True)
        vlog('Current memory: ',n=1,mem=True)
    #end def map_raw_data_onto_grid


    def plot_radial_raw(self,quants='all',kmax=None,fmt='b.',fig=True,show=True):
        data = self.get_raw_data()
        if quants=='all':
            quants = list(data.keys())
        #end if
        for q in quants:
            d = data[q]
            k  = np.linalg.norm(d.k,axis=1)
            nk = d.nk
            has_error = 'nk_err' in d
            if has_error:
                nke = d.nk_err
            #end if
            if kmax is not None:
                rng = k<kmax
                k   = k[rng]
                nk  = nk[rng]
                if has_error:
                    nke = nke[rng]
                #end if
            #end if
            if fig:
                plt.figure()
            #end if
            if not has_error:
                plt.plot(k,nk,fmt)
            else:
                plt.errorbar(k,nk,nke,fmt=fmt)
            #end if
            plt.xlabel('k (a.u.)')
            plt.ylabel('n(k) {}'.format(q))
        #end for
        if show:
            plt.show()
        #end if
    #end def plot_radial_raw


    def plot_directional_raw(self,kdir,quants='all',kmax=None,fmt='b.',fig=True,show=True,reflect=False):
        data = self.get_raw_data()
        kdir = np.array(kdir,dtype=float)
        kdir /= np.linalg.norm(kdir)
        if quants=='all':
            quants = list(data.keys())
        #end if
        for q in quants:
            d = data[q]
            k  = d.k
            nk = d.nk
            has_error = 'nk_err' in d
            if has_error:
                nke = d.nk_err
            #end if
            km = np.linalg.norm(d.k,axis=1)
            if kmax is not None:
                rng = km<kmax
                km  = km[rng]
                k   = k[rng]
                nk  = nk[rng]
                if has_error:
                    nke = nke[rng]
                #end if
            #end if
            kd = np.dot(k,kdir)
            along_dir = (np.abs(km-np.abs(kd)) < 1e-8*km) | (km<1e-8)
            kd = kd[along_dir]
            nk = nk[along_dir]
            if has_error:
                nke = nke[along_dir]
            #end if
            if fig:
                plt.figure()
            #end if
            if not has_error:
                plt.plot(kd,nk,fmt)
                if reflect:
                    plt.plot(-kd,nk,fmt)
                #end if
            else:
                plt.errorbar(kd,nk,nke,fmt=fmt)
                if reflect:
                    plt.errorbar(-kd,nk,nke,fmt=fmt)
                #end if
            #end if
            plt.xlabel('k (a.u.)')
            plt.ylabel('directional n(k) {}'.format(q))
        #end for
    #end def plot_directional_raw
#end class MomentumDistribution

MomentumDistribution.define_attributes(
    Observable,
    raw = obj(
        type       = obj,
        no_default = True,
        ),
    u = obj(
        type       = obj,
        no_default = True,
        ),
    d = obj(
        type       = obj,
        no_default = True,
        ),
    tot = obj(
        type       = obj,
        no_default = True,
        ),
    pol = obj(
        type       = obj,
        no_default = True,
        ),
    kaxes = obj(
        dest       = 'info',
        type       = np.ndarray,
        no_default = True,
        ),
    raw_filter_tol = obj(
        dest       = 'info',
        type       = float,
        default    = None,
        ),
    )



class MomentumDistributionDFT(MomentumDistribution):

    def read_eshdf(self,filepath,E_fermi=None,savefile=None,unfold=False,grid=True):

        save = False
        if savefile is not None:
            if os.path.exists(savefile):
                vlog('\nLoading from save file {}'.format(savefile))
                self.load(savefile)
                vlog('Done',n=1,time=True)
                return
            else:
                save = True
            #end if            
        #end if
                
        vlog('\nExtracting n(k) data from {}'.format(filepath))

        if E_fermi is None:
            E_fermi = self.info.E_fermi
        else:
            self.info.E_fermi = E_fermi
        #end if
        if E_fermi is None:
            self.error('Cannot read n(k) from ESHDF file.  Fermi energy (eV) is required to populate n(k) from ESHDF data.\nFile being read: {}'.format(filepath))
        #end if

        vlog.increment()
        d = read_eshdf_nofk_data(filepath,E_fermi)
        vlog.decrement()

        spins = {0:'u',1:'d'}

        spin_data = obj()
        for (ki,si) in sorted(d.data.keys()):
            vlog('Appending data for k-point {:>3} and spin {}'.format(ki,si),n=1,time=True)
            data = d.data[ki,si]
            s = spins[si]
            if s not in spin_data:
                spin_data[s] = obj(k=[],nk=[])
            #end if
            sdata = spin_data[s]
            sdata.k.extend(data.k)
            sdata.nk.extend(data.nk)
        #end for
        for sdata in spin_data:
            sdata.k  = np.array(sdata.k)
            sdata.nk = np.array(sdata.nk)
        #end for
        if 'd' not in spin_data:
            spin_data.d = spin_data.u
        #end if
        spin_data.tot = obj(
            k  = spin_data.u.k,
            nk = spin_data.u.nk + spin_data.d.nk,
            )

        self.set_attribute('raw'  ,spin_data)
        self.set_attribute('kaxes',d.kaxes  )

        if grid:
            self.map_raw_data_onto_grid(unfold=unfold)
        #end if

        if save:
            vlog('Saving to file {}'.format(savefile),n=1)
            self.save(savefile)
        #end if

        vlog('n(k) data extraction complete',n=1,time=True)

    #end def read_eshdf
#end class MomentumDistributionDFT

MomentumDistributionDFT.define_attributes(
    MomentumDistribution,
    E_fermi = obj(
        dest    = 'info',
        type    = float,
        default = None,
        )
    )


class MomentumDistributionQMC(MomentumDistribution):
    def read_stat_h5(self,*files,equil=0,savefile=None):

        save = False
        if savefile is not None:
            if os.path.exists(savefile):
                vlog('\nLoading from save file {}'.format(savefile))
                self.load(savefile)
                vlog('Done',n=1,time=True)
                return
            else:
                save = True
            #end if            
        #end if

        vlog('\nReading n(k) data from stat.h5 files',time=True)
        k   = []
        nk  = []
        nke = []
        if len(files)==1 and isinstance(files[0],(list,tuple)):
            files = files[0]
        #end if
        for file in files:
            if isinstance(file,StatFile):
                stat = file
            else:
                vlog('Reading stat.h5 file',n=1,time=True)
                stat = StatFile(file,observables=['momentum_distribution'])
            #end if
            vlog('Processing n(k) data from stat.h5 file',n=1,time=True)
            vlog('filename = {}'.format(stat.filepath),n=2)
            group = stat.observable_groups(self,single=True)

            kpoints = np.array(group['kpoints'])
            nofk    = np.array(group['value'])

            nk_mean,nk_var,nk_error,nk_kappa = simstats(nofk[equil:],dim=0)

            k.extend(kpoints)
            nk.extend(nk_mean)
            nke.extend(nk_error)
        #end for
        vlog('Converting concatenated lists to arrays',n=1,time=True)
        data = obj(
            tot = obj(
                k      = np.array(k),
                nk     = np.array(nk),
                nk_err = np.array(nke),
                )
            )
        self.set_attribute('raw',data)

        if save:
            vlog('Saving to file {}'.format(savefile),n=1)
            self.save(savefile)
        #end if

        vlog('stat.h5 file read finished',n=1,time=True)
    #end def read_stat_h5
#end class MomentumDistributionQMC



class Density(ObservableWithComponents):
    component_names = ('tot','pol','u','d')

    default_component = 'tot'


    def read_xsf(self,filepath,component=None):
        if isinstance(filepath,XsfFile):
            xsf = filepath
            copy_values = True
        else:
            xsf = XsfFile(filepath)
            copy_values = False
        #end if

        # read structure
        if not self.has_attribute('structure'):
            s = Structure()
            s.read_xsf(xsf)
            self.set_attribute('structure',s)
        #end if

        # read grid
        if not self.has_attribute('grid'):
            g = read_grid(xsf)
            self.set_attribute('grid',g)
            self.set_attribute('distance_units','B')
        #end if

        # read values
        xsf.remove_ghost()
        d = xsf.get_density()
        values = d.values_noghost.ravel()
        if copy_values:
            values = values.copy()
        #end if

        # create grid function for component
        f = grid_function(
            type   = 'parallelotope',
            grid   = self.grid,
            values = values,
            copy   = False,
            )

        component = self.process_component_name(component)
        self.set_attribute(component,f)
    #end def read_xsf

    
    def volume_normalize(self):
        g = self.get_attribute('grid',assigned=True)
        dV = g.volume()/g.ncells
        for c in self.components():
            c.values /= dV
        #end for
    #end def volume_normalize


    def change_distance_units(self,units):
        units_old = self.get_attribute('distance_units',assigned=True)
        rscale    = 1.0/convert(1.0,units_old,units)
        dscale    = 1./rscale**3
        grid      = self.get_attribute('grid',assigned=True)
        grid.points *= rscale
        for c in self.components():
            c.values *= dscale
        #end for
    #end def change_distance_units


    def change_density_units(self,units):
        units_old = self.get_attribute('density_units',assigned=True)
        dscale    = 1.0/convert(1.0,units_old,units)
        for c in self.components():
            c.values *= dscale
        #end for
    #end def change_density_units


    def radial_density(self,dr=0.01,ntheta=100,rmax=None,single=False,component=None,precompute=None):
        component,is_all = self.process_component_name(component,allow_all=True)
        pc = precompute
        s = self.get_attribute('structure',assigned=True)
        struct = s
        if rmax is None and pc is None:
            rmax = s.voronoi_species_radii()
        #end if
        if pc is None:
            equiv_atoms = s.equivalent_atoms()
            species = None
            species_rmax = obj()
            if isinstance(rmax,float):
                species = list(equiv_atoms.keys())
                for s in species:
                    species_rmax[s] = rmax
                #end for
            else:
                species = list(rmax.keys())
                species_rmax.transfer_from(rmax)
            #end if
            species_grids = obj()
            for s in species:
                nr = int(np.ceil(species_rmax[s]/dr))
                species_grids[s] = SpheroidGrid(
                    axes     = species_rmax[s]*np.eye(3),
                    cells    = (nr,ntheta,2*ntheta),
                    centered = True,
                    )
            #end for
            pc = obj(
                single        = single,
                equiv_atoms   = equiv_atoms,
                species_grids = species_grids,
                rmax          = species_rmax,
                )
                         
        #end if
        crdf = None
        if is_all:
            rdf = self.operate_over_components(self.cumulative_radial_density,precompute=pc)
        else:
            rdf = obj()
            d = self.component(component)
            for s,sgrid in pc.species_grids.items():
                rmax = pc.rmax[s]
                rrad    = sgrid.radii()
                rsphere = sgrid.r
                drad    = np.zeros(rrad.shape,dtype=d.dtype)
                if single:
                    atom_indices = [pc.equiv_atoms[s][0]]
                else:
                    atom_indices = pc.equiv_atoms[s]
                #end if
                rcenter = np.zeros((3,),dtype=float)
                for i in atom_indices:
                    new_center = struct.pos[i]
                    dr         = new_center-rcenter
                    rsphere   += dr
                    rcenter    = new_center
                    dsphere    = d.interpolate(rsphere)
                    dsphere.shape = sgrid.shape
                    dsphere.shape = len(dsphere),dsphere.size//len(dsphere)
                    drad += dsphere.mean(axis=1)*4*np.pi*rrad**2
                #end for
                drad /= len(atom_indices)
                rdf[s] = obj(
                    radius  = rrad,
                    density = drad,
                    )
            #end for
            d.clear_ghost()
        #end if
        return rdf
    #end def radial_density


    def plot_radial_density(self,component=None,show=True,**kwargs):
        component,is_all = self.process_component_name(component,allow_all=True)
        rdf = self.radial_density(component=component,**kwargs)
        if not is_all:
            rdf_all = obj()
            rdf_all[component] = rdf
        else:
            rdf_all = rdf
        #end if
        rdf = rdf_all.first()
        species = list(rdf.keys())

        for c in self.component_names:
            if c in rdf_all:
                rdf = rdf_all[c]
                for s in sorted(rdf.keys()):
                    srdf = rdf[s]
                    plt.figure()
                    plt.plot(srdf.radius,srdf.density,'b.-')
                    plt.xlabel('radius')
                    plt.ylabel('density')
                    plt.title('{} {} density'.format(s,c))
                #end for
            #end if
        #end for
        if show:
            plt.show()
        #end if
    #end def plot_radial_density
#end class Density

Density.define_attributes(
    Observable,
    raw = obj(
        type       = obj,
        no_default = True,
        ),
    u = obj(
        type       = obj,
        no_default = True,
        ),
    d = obj(
        type       = obj,
        no_default = True,
        ),
    tot = obj(
        type       = obj,
        no_default = True,
        ),
    pol = obj(
        type       = obj,
        no_default = True,
        ),
    grid = obj(
        type       = StructuredGrid,
        no_default = True,
        ),
    distance_units = obj(
        dest = 'info',
        type = str,
        ),
    density_units = obj(
        dest = 'info',
        type = str,
        )
    )



class ChargeDensity(Density):
    None
#end class ChargeDensity


class EnergyDensity(Density):
    None
#end class EnergyDensity




class StatFile(DevBase):

    scalars = set('''
        LocalEnergy   
        LocalEnergy_sq
        Kinetic       
        LocalPotential
        ElecElec      
        IonIon        
        LocalECP      
        NonLocalECP   
        KEcorr        
        MPC           
        '''.split())

    observable_aliases = obj(
        momentum_distribution = ['nofk'],
        )
    for observable in list(observable_aliases.keys()):
        for alias in observable_aliases[observable]:
            observable_aliases[alias] = observable
        #end for
        observable_aliases[observable] = observable
    #end for

    observable_classes = obj(
        momentum_distribution = MomentumDistributionQMC,
        )

    observable_class_to_stat_group = obj()
    for name,cls in observable_classes.items():
        observable_class_to_stat_group[cls.__name__] = name
    #end for


    def __init__(self,filepath=None,**read_kwargs):
        self.filepath = None

        if filepath is not None:
            self.filepath = filepath
            self.read(filepath,**read_kwargs)
        #end if
    #end def __init__

            
    def read(self,filepath,observables='all'):
        if not os.path.exists(filepath):
            self.error('Cannot read file.\nFile path does not exist: {}'.format(filepath))
        #end if
        h5 = h5py.File(filepath,'r')
        observable_groups = obj()
        for name,group in h5.items():
            # Skip scalar quantities
            if name in self.scalars:
                continue
            #end if
            # Identify observable type by name, for now
            for alias,observable in self.observable_aliases.items():
                cond_name  = self.condense_name(name)
                cond_alias = self.condense_name(alias)
                if cond_name.startswith(cond_alias):
                    if observable not in observable_groups:
                        observable_groups[observable] = obj()
                    #end if
                    observable_groups[observable][name] = group
                #end if
            #end for
        #end for
        if isinstance(observables,str):
            if observables=='all':
                self.transfer_from(observable_groups)
            #end if
        else:
            for obs in observables:
                if obs in observable_groups:
                    self[obs] = observable_groups[obs]
                #end if
            #end for
        #end if
    #end def read


    def condense_name(self,name):
        return name.lower().replace('_','')
    #end def condenst_name


    def observable_groups(self,observable,single=False):
        if inspect.isclass(observable):
            observable = observable.__name__
        elif isinstance(observable,Observable):
            observable = observable.__class__.__name__
        #end if
        groups = None
        if observable in self:
            groups = self[observable]
        elif observable in self.observable_class_to_stat_group:
            observable = self.observable_class_to_stat_group[observable]
            if observable in self:
                groups = self[observable]
            #end if
        #end if
        if single and groups is not None:
            if len(groups)==1:
                return groups.first()
            else:
                self.error('Single stat.h5 observable group requested, but multiple are present.\nGroups present: {}'.format(sorted(groups.keys())))
            #end if
        else:
            return groups
        #end if
    #end def observable_groups
#end class StatFile
