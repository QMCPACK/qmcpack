#!/usr/bin/env python3

#Tool to edit old Multideterminant and PBC HDF5 for LCAO (Can be extended to any HDF5 file)
#Takes as argument the name of the HDF5 file 
#Will transform the group KPTS_0 into Super_twist 
#Will change the Value of SpinUnResticted to SpinRestricted



import h5py
import sys

file_name = sys.argv[1]
f = h5py.File(file_name, 'r+')     # open the file


#Chenge the name of KPTS_0 to Super_Twist
#f.move('KPTS_0','Super_Twist')

# Change the SpinUnrestricted to SpinRestricted tag 
group = f['parameters']
#group.move('SpinUnResticted','SpinRestricted')

# fix the value of SpinRestricted to !SpinUnrestricted

data = group['SpinRestricted']  # load the data
if (data == False):
  new_value = True 
else:
  new_value = False

#data[...] = new_value  # assign new values to data
data[...] = True  # assign new values to data
f.close()  # close the file

