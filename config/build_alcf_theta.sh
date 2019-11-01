module unload cray-libsci
module load cray-hdf5-parallel
module load gcc
module load cmake

export CC=cc
export CXX=CC
export BOOST_ROOT=/soft/libraries/boost/1.64.0/intel
export CRAYPE_LINK_TYPE=dynamic

#TYPE=RelWithDebInfo
TYPE=Release
Compiler=Intel

for name in real_AoS_legacy real_MP_AoS_legacy cplx_AoS_legacy cplx_MP_AoS_legacy \
            real real_MP cplx cplx_MP
do

CMAKE_FLAGS="-D CMAKE_BUILD_TYPE=$TYPE"

if [[ $name == *"cplx"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D QMC_COMPLEX=1"
fi

if [[ $name == *"_AoS_legacy"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D ENABLE_SOA=0"
fi

if [[ $name == *"_MP"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D QMC_MIXED_PRECISION=1"
fi

folder=build_KNL_${Compiler}_${name}
echo "**********************************"
echo "$folder"
echo "$CMAKE_FLAGS"
echo "**********************************"
mkdir $folder
cd $folder
if [ ! -f CMakeCache.txt ] ; then
cmake $CMAKE_FLAGS ..
fi
make -j32
cd ..

echo
done
