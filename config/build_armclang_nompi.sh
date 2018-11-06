# This configuration was tested with Arm Compiler for HPC 18.4

# ARMPL_DIR is the path to the ARM Performance Libraries
#export ARMPL_DIR=

cmake -DCMAKE_C_COMPILER=armclang -DCMAKE_CXX_COMPILER=armclang++ -DQMC_MPI=0 \
      -DCMAKE_C_COMPILER_ID=Clang -DCMAKE_CXX_COMPILER_ID=Clang -DCMAKE_CXX_COMPILER_VERSION=5.0 -DCMAKE_CXX_STANDARD_COMPUTED_DEFAULT=98 \
      -DCMAKE_C_FLAGS="-march=armv8-a" -DCMAKE_CXX_FLAGS="-march=armv8-a" \
      -DLAPACK_LIBRARIES="-L$ARMPL_DIR/lib -larmpl_mp" \
      -DFFTW_INCLUDE_DIR="$ARMPL_DIR/include" \
      -DFFTW_LIBRARIES="$ARMPL_DIR/lib/libarmpl_mp.a" \
      ..
