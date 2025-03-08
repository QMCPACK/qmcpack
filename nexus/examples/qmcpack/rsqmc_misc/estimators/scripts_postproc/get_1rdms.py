#!/usr/bin/env python

import numpy as np
import h5py
import sys

class OneReducedDensityMatrix:
    def __init__(self):
        self.hname     = 'none'   # path/to/qmcpack.stat.h5
        self.is_valid  = False
        self.rdm_u     = 0        # up   channel
        self.rdm_u_err = 0
        self.rdm_d     = 0        # down channel
        self.rdm_d_err = 0
    # End function

    # Read a qmcpack.stat.h5 file and read the average 1rdms into numpy arrays
    def read_qmcpack_stath5(self, hname):
        # Open the file, extract the arrays
        self.hname = hname
        f = h5py.File(self.hname, 'r')
        uparr = f['OneBodyDensityMatrices/number_matrix/u/value']
        dwarr = f['OneBodyDensityMatrices/number_matrix/d/value']

        # Report details, set up the numpy arrays
        dims = uparr.shape
        nblk = dims[0] # Number of QMC blocks
        dim1 = dims[1] # Number of basis states
        dim2 = dims[2] # Number of basis states
        self.rdm_u     = np.zeros((dim1, dim2), dtype=float)
        self.rdm_d     = np.zeros((dim1, dim2), dtype=float)
        self.rdm_u_err = np.zeros((dim1, dim2), dtype=float)
        self.rdm_d_err = np.zeros((dim1, dim2), dtype=float)

        # Assign the arrays - average over blocks
        for i in range(dim1):
            for j in range(dim2):
                self.rdm_u[i,j]     = np.mean(uparr[:,i,j])
                self.rdm_d[i,j]     = np.mean(dwarr[:,i,j])
                self.rdm_u_err[i,j] = np.std(uparr[:,i,j], ddof=1)
                self.rdm_d_err[i,j] = np.std(dwarr[:,i,j], ddof=1)

        self.rdm_u_err /= np.sqrt(nblk - 1)
        self.rdm_d_err /= np.sqrt(nblk - 1)

        # Enforce symmetry
        self.rdm_u     = 0.5*(self.rdm_u + self.rdm_u.T)
        self.rdm_d     = 0.5*(self.rdm_d + self.rdm_d.T)
        self.rdm_u_err = 0.5*(self.rdm_u_err + self.rdm_u_err.T)
        self.rdm_d_err = 0.5*(self.rdm_d_err + self.rdm_d_err.T)

        print("Read summary:")
        print("Norbs  = {:4d}".format(dim1))
        print("Nblocks= {:4d}".format(nblk))
        print("Trace of 1rdms:")
        print("   spin 0 = {:12.8f} +/- {:12.8f}".format(np.trace(self.rdm_u), np.trace(self.rdm_u_err)))
        print("   spin 1 = {:12.8f} +/- {:12.8f}".format(np.trace(self.rdm_d), np.trace(self.rdm_d_err)))
        
        self.is_valid = True
    # End function


    # Print the array
    def pprint(self, channel):
        if channel == "up":
            for i in range(self.rdm_u.shape[0]):
                txt = ""
                for j in range(self.rdm_u.shape[1]):
                    txt += "{:12.8e} ".format(self.rdm_u[i,j])
                print(txt)
        elif channel == "down":
            for i in range(self.rdm_d.shape[0]):
                txt = ""
                for j in range(self.rdm_d.shape[1]):
                    txt += "{:12.8e} ".format(self.rdm_d[i,j])
                print(txt)
    # End function


    # Write the array and errors
    def write(self, channel):
        if channel == "up":
            fa  = open("1rdm_"+channel+".dat", mode="w")
            fe  = open("1rdm_"+channel+"_err.dat", mode="w")
            arr = np.copy(self.rdm_u)
            err = np.copy(self.rdm_u_err)
        else:
            fa  = open("1rdm_"+channel+".dat", mode="w")
            fe  = open("1rdm_"+channel+"_err.dat", mode="w")
            arr = np.copy(self.rdm_d)
            err = np.copy(self.rdm_d_err)
            
        print("Writing file: ", fa)
        print("Writing file: ", fe)
        for i in range(self.rdm_u.shape[0]):
            txta = ""
            txte = ""
            for j in range(self.rdm_u.shape[1]):
                txta += "{:12.8e} ".format(arr[i,j])
                txte += "{:12.8e} ".format(err[i,j])
            fa.write(txta + "\n")
            fe.write(txte + "\n")
        fa.close()
        fe.close()
    # End function
# End class
        


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: ./get_1rdms.py path/to/qmcpack.stat.h5")
        sys.exit()
    
    print("Starting job...")
    rdms = OneReducedDensityMatrix()
    rdms.read_qmcpack_stath5(sys.argv[1])
    if rdms.is_valid:
        rdms.write("up")
        rdms.write("down")
    else:
        print("It failed!")
        sys.exit()
    
    print("Done.")
