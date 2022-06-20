# ESHDF Python class

ESHDF (Electronic Structure Hierarchical Data Format) is the file
specification used to provide orbital information for QMCPACK. It is based on
the Hierarchical Data Format version 5 (HDF5; 
https://www.hdfgroup.org/solutions/hdf5/).

A Python class for ESHDF is defined in ´eshdf.py´, using ´h5py´ Python
library. The class hierarchy is intended for reading, managing and writing of
ESHDF data.

The classes are inherited as follows: `Base -> Representation -> Application.`
The Base class ´EshdfFileBase´ contains virtual implementation of the core
methods. The Representation classes add methods and variables for managing
specific choices of the representation. The Application class implements
bindings to a specific application.


## Currently implemented classes

Base:
- EshdfFileBase
Representation:
- EshdfFilePw (plane-wave)
Application:
- EshdfFilePwGpaw (´Eshdf_gpaw.py´)

This is work-in-progress, and other features are coming up.


# GPAW4QMCPACK

Converter of electronic orbitals from GPAW (https://wiki.fysik.dtu.dk/gpaw/)
is implemented in ´gpaw4qmcpack.py´, using the ESHDF Python class.

The converter is run on the command line. It requires a GPAW restart file
(.gpw) with orbitals included with ´mode=all´, and produces an ESHDF file (.h5)
that can be used in QMCPACK.

If the MPC correction will be needed, the charge density information can be 
added with `-d`, but it must be computed on the runtime. Therefore, the 
density is not included by default.

This is work-in-progress, and other features and documentation are coming up.
