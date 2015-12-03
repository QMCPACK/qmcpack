#!/bin/sh

# Attempt to automatically download and patch Quantum-Espresso for pw2qmcpack converter 
# Patch developed by William Parker / Argonne National Lab
# This simple script, Paul Kent / Oak Ridge National LaB

codename=espresso-5.2.1
archivename=${codename}.tar.gz
if [ ! -e ${archivename} ]; then
echo --- Did not find ${archivename} in current directory.
# Full URL of espresso download link obtained from qe-forge on 29 September 2014
# Will need to be updated between versions
qeurl=http://qe-forge.org/gf/download/frsrelease/199/855/${archivename}
echo --- Attempting to download. This can take several minutes. 
wget --no-verbose ${qeurl}
else
echo --- Found and using ${archivename} in current directory.
fi
if [ ! -e ${archivename} ]; then
echo --- ERROR: Could not find ${archivename}
echo --- Something went wrong... possibly a bad URL for the file download or you are offline
echo --- Please advise QMCPACK Developers via Google Groups if this problem persists
else
echo --- Unpacking
tar xvzf ${archivename}
if [ ! -e ${codename}/PW/src/Makefile ]; then
echo --- ERROR: Could not find PW/src/Makefile
echo --- Something went wrong... probably a failure to download the full archive.
echo --- Check ${archivename}. Delete if a partial download and retry.
echo --- Also check $qeurl is valid - perhaps the files have moved online.
echo --- Please advise QMCPACK Developers via Google Groups if this problem persists
else
cd ${codename}
patch -p1 -i ../add_pw2qmcpack_to_${codename}.diff
cd ..
echo --- SUCCESS: ${codename} patched for pw2qmcpack converter
echo "--- Configure using ./configure --with-hdf5 HDF5_DIR=(HDF5 base directory)"
echo --- Add platform specific options as needed
fi

fi

