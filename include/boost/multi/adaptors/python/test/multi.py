# /home/correaa/boost-multi/.venv/bin/python /home/correaa/boost-multi/include/boost/multi/adaptors/python/test/multi.py
import cppyy  # automatic C++ bindings provided by cppyy
import numpy as np  # to demonstrate interoperability
import os

# bind Multi to Python
# cppyy.add_include_path('/home/correaa/boost-multi/include')
cppyy.add_include_path('./include')
cppyy.include("boost/multi/array.hpp")
cppyy.include("boost/multi/io.hpp")

# short name
multi = cppyy.gbl.boost.multi

a2d = multi.array['double', 2]();  # create an empty array
a2d.assign([[1.0, 2.0],[3.0, 4.0]])  # assign elements
print(a2d)

# create an array
a2d = multi.array['double', 2]([[1.0, 2.0],[4.0, 5.0]])
print(a2d)

# create a numpy view to the array
npa2d = np.frombuffer(a2d.base(), dtype=np.float64, count=a2d.num_elements()).reshape(a2d.sizes().get[0](), a2d.sizes().get[1]())
print(npa2d)

# create an array_ref to a numpy array
a2d_ref = multi.array_ref['double', 2](multi.extensions_t[2](*npa2d.shape), npa2d)
print(a2d_ref)

# modify element
a2d_ref[1][1] = 999.9

# show changes
print(npa2d)
print(a2d)
