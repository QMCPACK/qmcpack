
# Create HDF file with harmonic potential on a grid

import h5py
import numpy as np

f = h5py.File('sho.h5','w')

start = -5.0
end = 5.0
num = 22
E = 0.3
m = 5.0
l = 1. / np.sqrt(m * E)
k = E / (l**2)
data = np.zeros( (num, num, num) )

delta = (end-start)/(num-1)

for ix in range(num):
  for iy in range(num):
    for iz in range(num):
      x = delta*ix + start
      y = delta*iy + start
      z = delta*iz + start
      r2 = x*x + y*y + z*z
      val = 0.5 * k * r2  
      # print x,y,z,r2,val
      data[ix,iy,iz] = val

f.create_dataset("pot_data", data=data)
