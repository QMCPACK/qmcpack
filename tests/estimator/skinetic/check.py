#!/usr/bin/env python
import numpy as np
import h5py
import xml.etree.ElementTree as ET

def grab_stat_entries(stat_file_name,name):
    """ stat_file_name: .stat.h5 file name, name: ex. gofr """
    stat_file = h5py.File(stat_file_name)

    data = []
    # find entry and extract data
    for estimator in stat_file.keys():
        if estimator != name:
            continue
        # end if
        entry = {"name":estimator}
        for key,value in stat_file[estimator].iteritems():
            entry[key] = value[:]
            if entry[key].shape == (1,):
                entry[key] = entry[key][0]
            # end if
        # end for
        data.append(entry)
    # end for
    if len(data) == 0:
        raise RuntimeError('no entry named %s was found in %s'%(name,stat_file_name))
    # end if
    return data
# end def

if __name__ == '__main__':

    est_name = 'skinetic'
    # read scalar.dat file with numpy.loadtxt (easier with pandas)
    fxml = 'vmc.xml'
    fscalar = 'bcc.s000.scalar.dat'
    fstat = 'bcc.s000.stat.h5'
    with open(fscalar,'r') as f:
        header = f.readline()
    # end with
    names = header.strip('#').split()
    data = np.loadtxt(fscalar)
    # data should be of shape (nstep,ncol), wheren col=len(names)

    skinetic_idx = [names.index(name) for name in names if name.startswith(est_name)]
    ktot = data[:,skinetic_idx].sum(axis=1) # total kinetic
    
    ref_idx = names.index('Kinetic')
    ktot_ref = data[:,ref_idx]
    # got ktot & ktot_ref

    # read input xml for VMC timestep (not absolutely necessary)
    sim_steps = np.arange(data.shape[0])
    inp = ET.parse(fxml)

    ts = float(inp.find('.//parameter[@name="timestep"]').text)
    sim_time = ts*sim_steps

    # read stat.h5 for skinetic entry
    stat_entry = grab_stat_entries(fstat,est_name)[0] # only 1 skinetic
    h5_data = stat_entry['value']
    
    if not np.allclose(h5_data,data[:,skinetic_idx]):
        print "species kinetic energy estimator h5 test failed"
        import matplotlib.pyplot as plt
        plt.plot(data[:,skinetic_idx],lw=1,c='k',label='scalar')
        plt.plot(h5_data,ls='--',lw=3,label='h5')
        plt.show()

    # test against 
    if not np.allclose(ktot,ktot_ref):
        print "species kinetic energy estimator test failed!"
        import matplotlib.pyplot as plt
        fig,ax = plt.subplots(1,1)
        ax.set_xlabel('simulation time (ha$^{-1}$)',fontsize=16)
        ax.set_ylabel('total kinetic (ha)',fontsize=16)
        ax.plot(sim_time,ktot,ls='--',lw=3,label='')
        ax.plot(sim_time,ktot_ref,lw=1,c='k')
        plt.show()

