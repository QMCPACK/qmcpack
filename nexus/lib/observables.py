

# Python standard library imports
import os
import inspect
from time import clock

# Non-standard Python imports
import numpy as np
import matplotlib.pyplot as plt
import h5py

# Nexus imports
from unit_converter import convert
from generic import obj
from developer import DevBase,log,error,ci
from numerics import simstats




class VLog(DevBase):

    verbosity_levels = obj(
        none = 0,
        low  = 1,
        high = 2,
        )

    def __init__(self):
        self.tstart    = clock()
        self.tlast     = self.tstart
        self.verbosity = self.verbosity_levels.low
        self.indent    = 0
    #end def __init__

    def __call__(self,msg,level='low',n=0,time=False):
        if self.verbosity==self.verbosity_levels.none:
            return
        elif self.verbosity >= self.verbosity_levels[level]:
            if time:
                tnow = clock()
                msg += '  (elapsed {:8.4f}, total {:8.4f})'.format(tnow-self.tlast,tnow-self.tstart)
                self.tlast = tnow
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



class Observable(DevBase):

    defaults = obj(
        info = obj,
        )

    info_defaults = obj()

    def __init__(self,**values):
        self.set_defaults()
        self.set_values(**values)
    #end def __init__


    def set_defaults(self):
        cls = self.__class__
        for name,value in cls.defaults.items():
            if inspect.isclass(value):
                value = value()
            #end if
            self[name] = value
        #end for
        info = self.info
        for name,value in cls.info_defaults.items():
            if inspect.isclass(value):
                value = value()
            #end if
            info[name] = value
        #end for
    #end def set_defaults


    def set_values(self,**values):
        cls = self.__class__
        value_names = set(values.keys())
        names       = set(cls.defaults.keys()) & value_names
        info_names  = set(cls.info_defaults.keys()) & value_names
        invalid     = value_names - names - info_names
        if len(invalid)>0:
            v = obj()
            v.transfer_from(values,invalid)
            self.error('Attempted to set unrecognized values\nUnrecognized values:\n{}'.format(v))
        #end if
        for name in names:
            self[name] = values[name]
        #end for
        info = self.info
        for name in info_names:
            info[name] = values[name]
        #end for
    #end def set_values
#end class Observable





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



class MomentumDistribution(Observable):
    spins = ('u','d')

    def plot_radial_raw(self,quants='all',kmax=None,fmt='b.',fig=True,show=True):
        data = self.raw
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
        data = self.raw
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



class MomentumDistributionDFT(MomentumDistribution):

    info_defaults = obj(
        E_fermi = None,
        **MomentumDistribution.defaults
        )


    def read_eshdf(self,filepath,E_fermi=None,savefile=None,grid=True):
        from grid_functions import grid_function

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
        self.raw = spin_data

        if grid:
            for s,sdata in spin_data.items():
                vlog('Mapping spin {} data onto grid'.format(si),n=1,time=True)
                self[s] = grid_function(
                    points = sdata.k,
                    values = sdata.nk,
                    axes   = d.kaxes,
                    )
            #end for
            if 'd' not in self:
                self.d = self.u
            #end if
        #end if

        if save:
            vlog('Saving to file {}'.format(savefile),n=1)
            self.save(savefile)
        #end if

        vlog('n(k) data extraction complete',n=1,time=True)

    #end def read_eshdf
#end class MomentumDistributionDFT



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
        self.raw = obj(
            tot = obj(
                k      = np.array(k),
                nk     = np.array(nk),
                nk_err = np.array(nke),
                )
            )

        if save:
            vlog('Saving to file {}'.format(savefile),n=1)
            self.save(savefile)
        #end if

        vlog('stat.h5 file read finished',n=1,time=True)
    #end def read_stat_h5
#end class MomentumDistributionQMC



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
