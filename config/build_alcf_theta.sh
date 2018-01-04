module unload cray-libsci
module load cray-hdf5-parallel

export CC=cc
export CXX=CC
export BOOST_ROOT=/soft/libraries/boost/1.64.0/intel
export CRAYPE_LINK_TYPE=dynamic

#TYPE=RelWithDebInfo
TYPE=Release
Compiler=Intel

for name in KNL_${Compiler}_real_SoA KNL_${Compiler}_real_SoA_MP \
            KNL_${Compiler}_cplx_SoA KNL_${Compiler}_cplx_SoA_MP \
            KNL_${Compiler}_real KNL_${Compiler}_real_MP \
            KNL_${Compiler}_cplx KNL_${Compiler}_cplx_MP
do

CMAKE_FLAGS="-D CMAKE_BUILD_TYPE=$TYPE"

if [[ $name == *"_MP"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D QMC_MIXED_PRECISION=1"
fi

if [[ $name == *"_cplx"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D QMC_COMPLEX=1"
fi

if [[ $name == *"_SoA"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D ENABLE_SOA=1"
fi

folder=build_${name}
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
