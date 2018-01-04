#!/bin/bash

Compiler=Clang++11

for name in real_SoA real_SoA_MP cplx_SoA cplx_SoA_MP \
            real real_MP cplx cplx_MP
do

CMAKE_FLAGS="-D CMAKE_TOOLCHAIN_FILE=../config/BGQ_${Compiler}_ToolChain.cmake"

if [[ $name == *"cplx"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D QMC_COMPLEX=1"
fi

if [[ $name == *"_SoA"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D ENABLE_SOA=1"
fi

if [[ $name == *"_MP"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D QMC_MIXED_PRECISION=1"
fi

folder=build_${Compiler}_${name}
echo "**********************************"
echo "$folder"
echo "$CMAKE_FLAGS"
echo "**********************************"
mkdir $folder
cd $folder
if [ ! -f CMakeCache.txt ] ; then
cmake $CMAKE_FLAGS ..
cmake $CMAKE_FLAGS ..
fi
make -j24
cd ..

echo
done
