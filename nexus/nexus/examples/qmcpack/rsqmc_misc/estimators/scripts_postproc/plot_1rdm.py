#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

def plot_matrix(arr, figname="arr"):
    print("Creating figure: ", figname)
    myarr = np.copy(arr)
    for i in range(myarr.shape[0]):
        for j in range(myarr.shape[1]):
            myarr[i,j] = np.log10(abs(arr[i,j]))
    plt.imshow(myarr, cmap=plt.cm.viridis, origin='upper', vmin=-6, vmax=0)
    plt.colorbar(label=r"log$_{10}~|\rho|$")

    ax = plt.gca()

    ax.set_title("One Body Reduced Density Matrix")
    ax.set_xlabel("Orbital Index")
    ax.set_ylabel("Orbital Index")
    
    ax.set_xticks(np.arange(0,arr.shape[0],5))   # major ticks
    ax.set_yticks(np.arange(0,arr.shape[0],5))
    ax.set_xticklabels(np.arange(1,arr.shape[0]+1,5))
    ax.set_yticklabels(np.arange(1,arr.shape[1]+1,5))
    ax.set_xticks(np.arange(-0.5,arr.shape[0],1), minor=True)
    ax.set_yticks(np.arange(-0.5,arr.shape[1],1), minor=True)
    ax.tick_params(axis='both', which="both", length=0)
    plt.grid(False)
    plt.grid(which="minor", color="k", linestyle="-", linewidth=1)
    
    plt.savefig(figname + ".pdf", bbox_inches="tight")
    plt.show()
    plt.close()
    return 0

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: plot_rdm path/to/rdm_spin.dat")
        sys.exit()

    rdm = np.loadtxt(sys.argv[1], unpack=False, dtype=float)
    plot_matrix(rdm, "rdm_grid")
