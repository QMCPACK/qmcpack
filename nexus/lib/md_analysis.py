def md_statistics(self,equil=None,autocorr=None):
    import numpy as np
    from numerics import simstats,simplestats
    mds = obj()
    for q,v in self.md_data.items():
        if equil is not None:
            v = v[equil:]
        else:
            equil=0
        if autocorr is None:
            mean,var,error,kappa = simstats(v)
            mds[q] = mean,error,kappa,equil
        else:
            nv = len(v)
            nb = int(np.floor(float(nv)/autocorr))
            nexclude = nv-nb*autocorr
            v = v[nexclude:]
            v.shape = nb,autocorr
            mean,error = simplestats(v.mean(axis=1))
            mds[q] = mean,error,autocorr,nexclude
    return mds


def md_plots(self,filename,show=True):
    md = self.md_data
    import numpy as np
    import matplotlib.pyplot as plt
    plt.figure()
    plt.subplots_adjust(left=0.2)
    plt.subplots_adjust(bottom=0.15)
    plt.suptitle(filename)
    plt.axis([None, None, -.1, .1])
    plt.subplot(3,1,1)
    plt.ylim(-1.0,1.0)
    plt.plot(md.time,md.total_energy-np.mean(md.total_energy),label='Etot')
    plt.plot(md.time,md.kinetic_energy-np.mean(md.kinetic_energy),label='Ekin')
    plt.plot(md.time,md.potential_energy-np.mean(md.potential_energy),label='Epot')
    plt.ylabel('E (Ryd)')
    plt.legend()
    plt.axis([None, None, None,None])
    plt.subplot(3,1,2)
    plt.plot(md.time,md.temperature)
    plt.ylabel('T (K)')
    plt.subplot(3,1,3)
    plt.plot(md.time,md.pressure)
    plt.ylabel('P (kbar)')
    plt.xlabel('time (ps)')
    if show:
        plt.show()
    return fig
