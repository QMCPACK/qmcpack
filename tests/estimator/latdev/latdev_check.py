#! /usr/bin/env python3

import numpy as np
import h5py

def print_fail_2d(a1_name, a1, a2_name, a2):
  close = np.isclose(a1, a2)
  print('  Index  %s  %s  Difference'%(a1_name, a2_name))
  for i in range(close.shape[0]):
    for j in range(close.shape[1]):
      if not close[i,j]:
        print(' ',i,j,a1[i,j],a2[i,j],abs(a1[i,j]-a2[i,j]))


if __name__ == '__main__':

    natom = 1
    ndim  = 3

    fdat  = 'bcc.s000.scalar.dat'
    fstat = 'bcc.s000.stat.h5'

    # get particle-averaged latdev from scalar.dat
    sdata = np.loadtxt(fdat)
    with open(fdat,'r') as f:
        file_data = f.readlines()

    # Gather columns into list, removing the '#' from the beginning
    cols = file_data[0].strip('#').split()

    # Initialize empty dictionary
    data = {}

    # Initialize columns
    for col in cols:
        data[col] = []

    # Add data from file to columns
    for row in file_data[1:]:
        row_data = row.split()
        for col, datum in zip(cols, row_data):
            data[col] += [float(datum)]


    nblock = sdata.shape[0]

    # get particle-resolved latdev from stat.dat
    fp = h5py.File(fstat)
    # The trailing [:] converts the Dataset to numpy array
    latdev = fp['latdev/value'][:]
    latdir = latdev.reshape(nblock,natom,ndim).mean(axis=1)
    lat_cols = [col for col in list(data.keys()) if col.startswith('latdev')]

    # Pull relevant data (the columns in lat_cols) out of the dictionary that was read in
    slatdir = []
    for col in lat_cols:
        slatdir += [data[col]]

    # slatdir is originally of shape (3, 100), while latdir is of shape (100, 3)
    # the transpose operation makes slatdir have shape (100, 3)
    slatdir = np.array(slatdir, dtype=np.float64).transpose()

    passed = True
    # check h5 against scalar.dat
    if not np.allclose(latdir,slatdir):
        print("lattice deviation estimator test failed - the values in the HDF file do not match the values in the .scalar.dat file")
        print_fail_2d("h5 latdev",latdir,"scalar.dat latdev",slatdir)
        passed = False

        #import matplotlib.pyplot as plt
        #fig,ax = plt.subplots(1,1)
        #ax.plot(latdir,ls='--',lw=2)
        #ax.plot(slatdir)
        #plt.show()
    # end if

    if passed:
        exit(0)
    else:
        exit(1)
    # end if

# end __main__
