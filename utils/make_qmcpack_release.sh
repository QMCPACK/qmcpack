#!/bin/sh

#
# Naive script to create a qmcpack release tarball in $HOME via svn export
#
# /tmp used for workspace
#
# Tarballs are labeled by svn revision and YearMonthDay
# svn revision is parsed from svn export log (brittle hack)
#
# Uses gnutar on Mac (from macports etc. OS X 10.9+)
#
# No files are modified with exported revision number. Should update QMCPLUSPLUS_RELEASE, 
# QMCPLUSPLUS_BRANCH, QMCPLUSPLUS_LAST_CHANGED_DATE via CMakeLists.txt at svn detection
#
# P.Kent 2014/07
#

if [ ! -e /tmp/_tmp_qmcpack_release_work_dir ]; then
echo Using /tmp/_tmp_qmcpack_release_work_dir for temporary svn export workspace
mkdir /tmp/_tmp_qmcpack_release_work_dir
cd /tmp/_tmp_qmcpack_release_work_dir
svn export https://subversion.assembla.com/svn/qmcdev/trunk qmcpack >export.log
if [ -e qmcpack/utils/energy.pl ]; then

# Exclude this script from the release tarball
rm -f qmcpack/utils/make_qmcpack_release.sh

therev=`tail -1 export.log|sed -e 's/Exported revision //g' -e 's/\.//g'`
echo Parsed revision $therev
thedatedashed=`date "+%Y-%m-%d"`
thedate=`date "+%Y%m%d"`
echo Date is $thedate

# Fragile editing of version and dates in CMakeLists.txt
# Requires an exact match
cd qmcpack
cp -p -f CMakeLists.txt CMakeLists.txt_orig
oldrevstr='QMCPLUSPLUS_BRANCH "00000"'
newrevstr='QMCPLUSPLUS_BRANCH "'$therev'"'
olddatestr='QMCPLUSPLUS_LAST_CHANGED_DATE "2001-01-01"'
newdatestr='QMCPLUSPLUS_LAST_CHANGED_DATE "'$thedatedashed'"'
sed -e "s/$oldrevstr/$newrevstr/g" -e "s/$olddatestr/$newdatestr/g" CMakeLists.txt_orig >CMakeLists.txt
#diff CMakeLists.txt_orig CMakeLists.txt
rm -f CMakeLists.txt_orig
cd ..

OS=`uname`
if [ "$OS" == "Darwin" ]; then
echo Using Mac gnutar workaround to avoid extended attributes. Install from macports
gnutar cvf qmcpack_rev${therev}_${thedate}.tar qmcpack
else
tar cvf qmcpack_rev${therev}_${thedate}.tar qmcpack
fi
gzip qmcpack_rev${therev}_${thedate}.tar
md5 qmcpack_rev${therev}_${thedate}.tar.gz > qmcpack_rev${therev}_${thedate}.tar.gz.md5
cat qmcpack_rev${therev}_${thedate}.tar.gz.md5
mv qmcpack_rev${therev}_${thedate}.tar.gz.md5 qmcpack_rev${therev}_${thedate}.tar.gz $HOME
echo Release tar.gz and md5 checksum in $HOME
ls -l $HOME/qmcpack_rev${therev}_${thedate}.tar.gz
ls -l $HOME/qmcpack_rev${therev}_${thedate}.tar.gz.md5
# Tidy workspace
rm -r -f /tmp/_tmp_qmcpack_release_work_dir
else
echo No utils/energy.pl file found
echo Sanity check failed. Likely assembla export or authentication problem
echo ABORT
exit 1
fi
else
echo /tmp/_tmp_qmcpack_release_work_dir already exists
echo Delete or rename to proceed
echo ABORT
exit 1
fi

