#!/usr/bin/env python
import numpy as np
import pandas as pd
import h5py

if __name__ == '__main__':

    natom = 1
    ndim  = 3

    fdat  = 'bcc.s000.scalar.dat'
    fstat = 'bcc.s000.stat.h5'

    # get particle-averaged latdev from scalar.dat
    sdata = np.loadtxt(fdat)
    with open(fdat,'r') as f:
        header = f.readline().strip('#').split()
    df = pd.DataFrame(sdata,columns=header)
    nblock = sdata.shape[0]

    # get particle-resolved latdev from stat.dat
    fp = h5py.File(fstat)
    latdev = fp['latdev/value'].value
    latdir = latdev.reshape(nblock,natom,ndim).mean(axis=1)
    lat_cols = [col for col in df.columns if col.startswith('latdev')]
    slatdir  = df.loc[:,lat_cols].values

    if not np.allclose(latdir,slatdir):
        print "lattice deviation estimator test failed"

        import matplotlib.pyplot as plt
        fig,ax = plt.subplots(1,1)
        ax.plot(latdir,ls='--',lw=2)
        ax.plot(slatdir)
        plt.show()
