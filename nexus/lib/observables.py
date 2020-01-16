

import numpy as np

from unit_converter import convert
from generic import obj
from developer import DevBase,log,error,ci



class Observable(DevBase):
    None
#end class Observable





def read_eshdf_nofk_data(filename,Ef):
    from numpy import array,pi,dot,sqrt,abs,zeros
    from numpy.linalg import inv,det
    from hdfreader import read_hdf

    def h5int(i):
        return array(i,dtype=int)[0]
    #end def h5int

    E_fermi   = Ef + 1e-8
    h        = read_hdf(filename,view=True)
    gvu      = array(h.electrons.kpoint_0.gvectors)
    axes     = array(h.supercell.primitive_vectors)
    kaxes    = 2*pi*inv(axes).T
    gv       = dot(gvu,kaxes)
    Ngv      = len(gv[:,0])
    kmag     = sqrt((gv**2).sum(1))
    nk       = h5int(h.electrons.number_of_kpoints)
    ns       = h5int(h.electrons.number_of_spins)
    occpaths = obj()
    data     = obj()
    for k in range(nk):
        kin_k   = obj()
        eig_k   = obj()
        k_k     = obj()
        nk_k    = obj()
        nelec_k = zeros((ns,),dtype=float)
        kp = h.electrons['kpoint_'+str(k)]
        gvs = dot(array(kp.reduced_k),kaxes)
        gvk = gv.copy()
        for d in range(3):
            gvk[:,d] += gvs[d]
        #end for
        kinetic=(gvk**2).sum(1)/2 # Hartree units
        for s in range(ns):
            #print ' ',(k,s),(nk,ns)
            kin_s = []
            eig_s = []
            k_s   = gvk
            nk_s  = 0*kmag
            nelec_s = 0
            path = 'electrons/kpoint_{0}/spin_{1}'.format(k,s)
            spin = h.get_path(path)
            eig = convert(array(spin.eigenvalues),'Ha','eV')
            nst = h5int(spin.number_of_states)
            for st in range(nst):
                e = eig[st]
                if e<E_fermi:
                    stpath = path+'/state_{0}/psi_g'.format(st)
                    occpaths.append(stpath)
                    psi = array(h.get_path(stpath))
                    nk_orb = (psi**2).sum(1)
                    kin_orb = (kinetic*nk_orb).sum()
                    nelec_s += nk_orb.sum()
                    nk_s += nk_orb
                    kin_s.append(kin_orb)
                    eig_s.append(e)
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
        nkpoints = nk,
        nspin    = ns,
        data     = data,
        )
    return res
#end def read_eshdf_nofk_data



class MomentumDistribution(Observable):
    None
#end class MomentumDistribution



class MomentumDistributionDFT(MomentumDistribution):
    def __init__(self,E_fermi=None):
        self.info = obj(
            E_fermi = E_fermi,
            )
    #end def __init__


    def read_eshdf(self,filepath,E_fermi=None):
        from grid_functions import grid_function
        if E_fermi is None:
            E_fermi = self.info.E_fermi
        #end if
        if E_fermi is None:
            self.error('Cannot read n(k) from ESHDF file.  Fermi energy (eV) is required to populate n(k) from ESHDF data.\nFile being read: {}'.format(filepath))
        #end if
        d = read_eshdf_nofk_data(filepath,E_fermi)

        f = grid_function(
            points = d.data[1,0].k,
            values = d.data[1,0].nk,
            axes   = d.kaxes,
            )

        #jtk mark current
        #  move bulk of "read_from_points" up into ParallelotopeGrid
        #  next merge all up/down nk/k data into single arrays
        #  make grid function for up/down nk

        ci()
    #end def read_eshdf
#end class MomentumDistributionDFT
