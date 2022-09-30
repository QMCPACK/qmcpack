module unload cray-libsci
module load cray-hdf5-parallel
module load gcc/9.3.0
module load cmake/3.20.4
module load intel/19.1.2.254

module list >& load_modules.txt

export CC=cc
export CXX=CC
export BOOST_ROOT=/soft/libraries/boost/1.64.0/intel
export CRAYPE_LINK_TYPE=dynamic

#TYPE=RelWithDebInfo
TYPE=Release
Compiler=Intel

if [[ $# -eq 0 ]]; then
  source_folder=`pwd`
elif [[ $# -eq 1 ]]; then
  source_folder=$1
else
  source_folder=$1
  install_folder=$2
fi


CURRENT_FOLDER=`pwd`

for name in real real_MP cplx cplx_MP
do

CMAKE_FLAGS="-D CMAKE_SYSTEM_NAME=CrayLinuxEnvironment -D CMAKE_BUILD_TYPE=$TYPE -D MPIEXEC_EXECUTABLE=/bin/sh -D MPIEXEC_NUMPROC_FLAG=$source_folder/tests/scripts/aprunhelper.sh"

if [[ $name == *"cplx"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D QMC_COMPLEX=1"
fi

if [[ $name == *"_MP"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D QMC_MIXED_PRECISION=1"
fi

folder=build_KNL_${Compiler}_${name}

if [[ -v install_folder ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DCMAKE_INSTALL_PREFIX=$install_folder/$folder"
fi

echo "**********************************"
echo "$folder"
echo "$CMAKE_FLAGS"
echo "**********************************"
mkdir $folder
cd $folder
if [ ! -f CMakeCache.txt ] ; then
  cmake $CMAKE_FLAGS $source_folder
fi
if [[ -v install_folder ]]; then
  make -j16 install && chmod -R -w $install_folder/$folder
else
  make -j16
fi
cd ..

echo
done
