
# Generate a constructor for a skew symmetric matrix from a list of parameters
# Output is Python code to be used in further scripts

# Autodiff expects matrices to be constructed and not updated.  This makes it
# difficult to use autodiff with programmatically constructed rotation matrices.

import numpy as np


# Rotation indices for a ground state
# This should produce the same list of indices as
# QMCWaveFunctions/RotatedSPOs.cpp::createRotationIndices
def get_rotation_indices(nel, nmo):
    rot_ind = list()
    for i in range(nel):
        for j in range(nel,nmo):
            rot_ind.append( (i,j) )

    return rot_ind

# Rotation indices including one excited state.
# This implementation is specific to the first excited state being
# included in one of the determinants.
# Should produce the same list of indices as
# QMCWaveFunctions/Fermion/MultiDiracDeterminant.cpp::buildOptVariables
# for this problem.
def get_rotation_indices_ex(nel, nmo):
    rot_ind = list()
    for i in range(nel+1):
        for j in range(i+1,nmo):
            rot_ind.append( (i,j) )

    return rot_ind

# Full rotation indices (for global rotation)
# This should produce the same list of indices as
# QMCWaveFunctions/RotatedSPOs.cpp::createRotationIndicesFull
def get_full_rotation_indices(nel, nmo):
    rot_ind = list()
    for i in range(nel):
        for j in range(nel,nmo):
            rot_ind.append( (i,j) )

    for i in range(nel):
        for j in range(i+1,nel):
            rot_ind.append( (i,j) )

    for i in range(nel,nmo):
        for j in range(i+1,nmo):
            rot_ind.append( (i,j) )

    return rot_ind


# Populate a skew symmetric matrix with corresponding indices
def construct_antisym(nmo,rot_ind):
   rot = np.zeros((nmo,nmo),dtype=np.int64)
   for idx in range(len(rot_ind)):
        p,q = rot_ind[idx]
        rot[q,p] = idx+1
        rot[p,q] = -(idx+1)

   return rot

# Convert the output of construct_antisym to a numpy array constructor
def print_anti(rot_mat):
    print("np.array([")
    for i in range(rot_mat.shape[0]):
        print("[ ",end='')
        for j in range(rot_mat.shape[1]):
            idx = np.abs(rot_mat[i,j])
            sign = ""
            if rot_mat[i,j] < 0:
                sign = "-"
            if idx == 0:
                print(" 0",end='')
            else:
                print(" {}p[{}]".format(sign,idx-1),end='')
            if j != rot_mat.shape[1]-1:
                print(",",end='')
        print("]",end='')
        if i != rot_mat.shape[0]-1:
            print(",")
    print("])")


if __name__ == "__main__"
    # For the Be problem
    nel = 2
    nmo = 7

    # For rot_be_sto_wf.py
    #rot_ind = get_rotation_indices(nel,nmo)

    # For rot_multi_be_sto_wf.py
    rot_ind = get_rotation_indices_ex(nel, nmo)

    #rot_ind = get_full_rotation_indices(nel,nmo)
    rot_mat = construct_antisym(nmo, rot_ind)
    print_anti(rot_mat)


