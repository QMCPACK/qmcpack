#!/bin/sh

# Attempt to automatically download and patch Quantum-Espresso for pw2qmcpack converter 
# Patch developed by William Parker / Argonne National Lab
# This simple script, Paul Kent / Oak Ridge National LaB

codename=qe-6.4.1
untarname=q-e-qe-6.4.1
archivename=${codename}.tar.gz
if [ ! -e ${archivename} ]; then
echo --- Did not find ${archivename} in current directory.
# Full URL of espresso download link obtained from qe-forge on 22 Feb 2018
# Will need to be updated between versions
qeurl=https://github.com/QEF/q-e/archive/${codename}.tar.gz
echo --- Attempting to download. This can take several minutes. 
wget --no-verbose ${qeurl}
else
echo --- Found and using ${archivename} in current directory.
fi

if [ ! -e ${archivename} ]; then
echo --- ERROR: Could not find ${archivename}
echo --- Something went wrong... possibly a bad URL for the file download or you are offline
echo --- Please advise QMCPACK Developers via Google Groups if this problem persists
exit
fi

if [ -e ${codename} ]; then
echo --- ERROR: folder ${codename} already exist! Could not unpack ${archivename}!
exit
fi

echo --- Unpacking
tar xvzf ${archivename}
mv $untarname $codename
if [ ! -e ${codename}/PW/src/Makefile ]; then
echo --- ERROR: Could not find PW/src/Makefile
echo --- Something went wrong... probably a failure to download the full archive.
echo --- Check ${archivename}. Delete if a partial download and retry.
echo --- Also check $qeurl is valid - perhaps the files have moved online.
echo --- Please advise QMCPACK Developers via Google Groups if this problem persists
exit
fi

cd ${codename}
patch -f -p1 -i ../add_pw2qmcpack_to_${codename}.diff
cd ..
if [ -e $codename/PP/src/pw2qmcpack.f90 ]; then
echo --- SUCCESS: ${codename} patched for pw2qmcpack converter
echo There are two ways to build
echo "1) if your system already has HDF5 installed with Fortran, use the --with-hdf5 configuration option."
echo "   Example configure: ./configure --with-hdf5=path_to_hdf5_include_and_lib"
echo "   Currently HDF5 support in QE itself is preliminary. To enable use of pw2qmcpack"
echo "   but use the old non-HDF5 I/O within QE, replace '-D__HDF5' with '-D__HDF5_C' in make.inc."
echo "2) if your system has HDF5 with C only, manually edit make.inc by adding '-D__HDF5_C -DH5_USE_16_API'"
echo "   in 'DFLAGS' and provide include and library path in 'IFLAGS' and 'HDF5_LIB'"
else
echo --- ERROR: Could not find PP/src/pw2qmcpack.f90 after patching
echo --- Probably the patch is missing or the archive has been updated.
echo --- Please advise QMCPACK Developers via Google Groups.
fi    
