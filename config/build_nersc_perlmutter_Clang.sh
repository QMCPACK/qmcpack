#!/bin/bash
# This recipe is intended for NERSC's Perlmutter supercomputer https://docs.nersc.gov/systems/perlmutter
# It builds all the variants of QMCPACK in the current directory
# This file was originally written by Ye Luo
# Last revision: Apr 24th 2026 by Brock Dyer

set -e

usage() {
  cat << EOF >&2
Usage: $0 [-s <source_dir>] [-i <install_dir>] [-e <espresso_dir>] [-p <python_env>] [-t <cpu|gpu>] [-fmrc] [-h]

-s <source_dir>      Location of the QMCPACK source code. (Defaults to working directory)
-i <install_dir>     Directory to install to (affects -DCMAKE_INSTALL_PREFIX).
-e <espresso_dir>    Location of the bin directory of a Quantum ESPRESSO build with pw2qmcpack.x (sets -DQE_BIN).
-p <python_env>      Location of a Python environment with Nexus installed (used to set -DPython3_EXECUTABLE).
-t <cpu|gpu>         Target CPU or GPU variants only, if not specified then build all (gpu sets -DQMC_GPU_ARCHS=sm_80).
-j <jobs>            Number of jobs to build with. Defaults to 16.
-f                   Only build Full Precision variants (disables Mixed Precision variants).
-m                   Only build Mixed Precision variants (disables Full Precision variants).
-r                   Only build Real variants (Disables Complex variants).
-c                   Only build Complex variants (Disables Real variants).
-h                   Display this help and exit.

This script will build all currently available variants of QMCPACK unless otherwise specified.

The current build variants are:

- CPU
  - Real
  - Real w/ Mixed Precision
  - Complex
  - Complex w/ Mixed Precision
- GPU
  - Real
  - Real w/ Mixed Precision
  - Complex
  - Complex w/ Mixed Precision

EOF
  exit 1
}


source_dir=`pwd`
build_cpu_gpu="ALL"
num_jobs=16
full_precision=1
mixed_precision=1
real_variants=1
complex_variants=1
while getopts s:i:e:p:t:j:fmrch o; do
  case $o in
    (s) source_dir="$(realpath $OPTARG)";;
    (i) install_dir="$(realpath $OPTARG)";;
    (e) qe_bin_dir="$(realpath $OPTARG)";;
    (p) py_venv="$(realpath $OPTARG)";;
    (t) build_cpu_gpu="$(echo $OPTARG | tr '[:lower:]' '[:upper:]')";;
    (j) num_jobs=$OPTARG;;
    (f) mixed_precision=0;;  # Full precision only
    (m) full_precision=0;;   # Mixed precision only
    (r) complex_variants=0;; # Real only
    (c) real_variants=0;;    # Complex only
    (h) usage;;
    (*) usage
  esac
done

if (( $real_variants == 0 && $complex_variants == 0 )); then
  echo -e "\nBoth real and complex variants are set to false, will not build!\n"
  sleep 1
  usage
elif (( $mixed_precision == 0 && $full_precision == 0 )); then
  echo -e "\nBoth Mixed and Full Precision variants are set to false, will not build!\n"
  sleep 1
  usage
fi

##########################################
##     Print user defined variables     ##
##########################################
echo -e "=---------------------------------------------------------------------------------=\n"
if [[ -f $source_dir/CMakeLists.txt ]]; then
  echo "Source folder set to '$source_dir'"
else
  echo "Source directory $source_dir doesn't contain CMakeLists.txt."
  echo "Pass QMCPACK source directory with '-s <source_dir>'"
  exit 1
fi

if [[ -v install_dir ]]; then
  echo "Install folder set to '$install_dir'"
else
  echo No install folder provided.
fi

if [[ -v qe_bin_dir ]]; then
  if [[ -e $qe_bin_dir/pw2qmcpack.x ]]; then
    echo "Quantum ESPRESSO bin directory set to '$qe_bin_dir'"
  else
    unset qe_bin_dir
    echo "Quantum ESPRESSO bin directory does not contain 'pw2qmcpack.x'"
    echo "  -> Will not add QE tests."
  fi
else
  echo No Quantum ESPRESSO bin directory provided
  echo "  -> Will not add QE tests."
fi
if [[ -v py_venv ]]; then
  if [[ ! -e $py_venv ]]; then
    echo "Can not find Python environment at path '$py_venv'"
    echo Exiting...
    exit 1
  else
    user_python=1
    echo "Using Python environment at '$py_venv'"
  fi
else
  user_python=0
  py_venv=$source_dir/.qmcvenv
  echo No Python environment provided
  echo "  -> Will check for one at '$py_venv'"
fi
echo
if [[ $build_cpu_gpu == "ALL" ]]; then
  echo "Building CPU and GPU variants"
else
  echo "Building $build_cpu_gpu variants"
fi

if (( $full_precision == 1 )); then
  echo "Building Full Precision variants : True"
else
  echo "Building Full Precision variants : False"
fi
if (( $mixed_precision == 1 )); then
  echo "Building Mixed Precision variants: True"
else
  echo "Building Mixed Precision variants: False"
fi
if (( $real_variants == 1 )); then
  echo "Building Real variants   : True"
else
  echo "Building Real variants   : False"
fi
if (( $complex_variants == 1 )); then
  echo "Building Complex variants: True"
else
  echo "Building Complex variants: False"
fi
echo
echo "To run QMCPACK built by this script,"
echo "you must load the same modules used during the build"
echo
echo "If you will be using the GPU-accelerated builds of QMCPACK,"
echo "you must set the following environment variable:"
echo "(bash)  export MPICH_GPU_SUPPORT_ENABLED=0"
echo "(zsh)   export MPICH_GPU_SUPPORT_ENABLED=0"
echo "(csh)   setenv MPICH_GPU_SUPPORT_ENABLED 0"
echo "(fish)  set -gx MPICH_GPU_SUPPORT_ENABLED 0"
echo
echo -e "\n=---------------------------------------------------------------------------------=\n"


##########################################
## Set up Python environment with Nexus ##
##########################################
echo -e "\n=---------------------------------------------------------------------------------="
echo -e "\nSetting up Python\n"

if [[ ! -e $py_venv ]]; then
  echo "Virtual environment not found at '$py_venv'"
  echo "  -> Creating Python virtual environment"
  module load python/3.13-26.1.0
  echo "$ python -m venv $py_venv"
  python -m venv $py_venv
  echo
fi

echo -e "Activating virtual environment"
echo -e "$ source $py_venv/bin/activate\n"
source $py_venv/bin/activate

python_executable="$(which python)"
if [[ $python_executable != "$py_venv/bin/python" ]]; then
  echo -e "\nPython executable is at $python_executable"
  echo "Expected $py_venv/bin/python"
  echo Virtual environment activation failed, exiting...
  exit
else
  echo "Python executable is at $python_executable"
fi
echo -e "\nChecking for Nexus installation..."
installed_packages="$(pip freeze)"
found_nexus=0
for line in $installed_packages; do
  if [[ $line == "nexus" ]]; then
    found_nexus=1
    echo Nexus installation found, will not re-install
    break
  fi
done

if [[ $found_nexus == 0 ]]; then
  if [[ $user_python == 1 ]]; then
    echo Nexus not found, ensure you have Nexus installed!
    exit 1
  fi
  cd $source_dir/nexus
  echo -e "\nNexus not found, installing Nexus and its dependencies into virtual environment..."
  echo -e "$ pip install .[full] --group dev\n"
  pip install .[full] --group dev
  cd -
fi

echo -e "\nDone setting up Python for build.\n"
echo -e   "=---------------------------------------------------------------------------------=\n"

if [[ $user_python == 0 ]]; then
  echo "To access Nexus analysis tools such as 'qmca', re-activate the virtual environment:"
  echo "  $ source $py_venv/bin/activate"
fi

##########################################
##    Assemble list of Build Targets    ##
##########################################

build_targets=""
if [[ $build_cpu_gpu == "CPU" ]] || [[ $build_cpu_gpu == "ALL" ]]; then
  if (( $real_variants == 1 )); then
    if (( $full_precision == 1 )); then
      build_targets="$build_targets cpu_real"
    fi
    if (( $mixed_precision == 1 )); then
      build_targets="$build_targets cpu_real_MP"
    fi
  fi
  if (( $complex_variants == 1 )); then
    if (( $full_precision == 1 )); then
      build_targets="$build_targets cpu_cplx"
    fi
    if (( $mixed_precision == 1 )); then
      build_targets="$build_targets cpu_cplx_MP"
    fi
  fi
fi
if [[ $build_cpu_gpu == "GPU" ]] || [[ $build_cpu_gpu == "ALL" ]]; then
  if (( $real_variants == 1 )); then
    if (( $full_precision == 1 )); then
      build_targets="$build_targets gpu_real"
    fi
    if (( $mixed_precision == 1 )); then
      build_targets="$build_targets gpu_real_MP"
    fi
  fi
  if (( $complex_variants == 1 )); then
    if (( $full_precision == 1 )); then
      build_targets="$build_targets gpu_cplx"
    fi
    if (( $mixed_precision == 1 )); then
      build_targets="$build_targets gpu_cplx_MP"
    fi
  fi
fi

#########################################
##        Load required modules        ##
#########################################

module load PrgEnv-gnu
module load cray-libsci
CRAY_LIBSCI_LIB=$CRAY_LIBSCI_PREFIX_DIR/lib/libsci_gnu_mp.so
module unload PrgEnv-gnu
module load craype cray-mpich
module load cray-fftw
module load cray-hdf5-parallel
module load cmake
export BOOST_ROOT=/global/common/software/m2113/opt/boost/1_82_0

TYPE=Release
Machine=perlmutter

for name in $build_targets; do

  CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=$TYPE -DBLAS_LIBRARIES=$CRAY_LIBSCI_LIB -DPython3_EXECUTABLE=$python_executable"

  if [[ $name == *"cplx"* ]]; then
    CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_COMPLEX=ON"
  fi

  if [[ $name == *"_MP"* ]]; then
    CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_MIXED_PRECISION=ON"
  fi

  if [[ $name == *"gpu"* ]]; then
    Compiler=Clang21
    CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_GPU_ARCHS=sm_80"
    echo -e "\n=---------------------------------------------------------------------------------=\n"
    reload_needed=1
    gpu_mod_check="$(module list gpu/1.0 2>&1)"
    for line in $gpu_mod_check; do
      if [[ $line == "None" ]]; then
        echo -e "Required Modules already loaded.\n"
        reload_needed=0
      fi
    done

    if [[ $reload_needed == 1 ]]; then
      echo "Loading the following modules:"
      echo "  - gpu/1.0"
      echo "  - PrgEnv-llvm"
      echo "  - llvm/21.1.4"
      echo "  - openmpi/5.0.7"
      echo

      module load gpu/1.0
      module load PrgEnv-llvm
      module load llvm/21.1.4
      module load openmpi/5.0.7
    fi

    C_compiler=`which mpicc`
    CXX_compiler=`which mpicxx`
    export MPICH_CC=$C_compiler
    export MPICH_CXX=$CXX_compiler

    echo -e "\n=---------------------------------------------------------------------------------=\n"
    echo -e "Clang version information:\n"
    echo "$ clang -v"
    clang -v
    echo -e "\n=---------------------------------------------------------------------------------=\n"
  else
    Compiler=GCC
    echo -e "\n=---------------------------------------------------------------------------------=\n"

    reload_needed=1
    cpu_mod_check="$(module list cpu/1.0 2>&1)"
    for line in $cpu_mod_check; do
      if [[ $line == "None" ]]; then
        echo -e "Required Modules already loaded.\n"
        reload_needed=0
      fi
    done

    if [[ $reload_needed == 1 ]]; then
      echo Loading the following modules:
      echo "  - cpu/1.0"
      echo "  - PrgEnv-gnu"
      echo "  - craype"
      echo "  - cray-mpich"
      echo

      module load cpu/1.0
      module load PrgEnv-gnu
      module load craype cray-mpich
    fi
    C_compiler=`which cc`
    CXX_compiler=`which CC`
  fi

  if [[ -v qe_bin_dir ]]; then
    CMAKE_FLAGS="$CMAKE_FLAGS -DQE_BIN=$qe_bin_dir"
  fi

  folder=build_${Machine}_${Compiler}_${name}
  if [[ -v install_dir ]]; then
    build_dir="$install_dir/$folder"
    CMAKE_FLAGS="$CMAKE_FLAGS -DCMAKE_INSTALL_PREFIX=$build_dir"
  else
    build_dir=$folder
  fi

  echo -e "=---------------------------------------------------------------------------------="
  echo -e "\nUsing flags:"
  for flag in $CMAKE_FLAGS; do
    echo "  $flag"
  done
  echo "  -DCMAKE_C_COMPILER=$C_compiler"
  echo "  -DCMAKE_CXX_COMPILER=$CXX_compiler"
  echo "  -S $source_dir"

  echo -e   "\nInstalling to $build_dir\n"
  echo -e   "=---------------------------------------------------------------------------------=\n"

  if [ ! -d $build_dir ]; then
    mkdir $build_dir
  fi
  cd $build_dir

  if [ ! -f CMakeCache.txt ] ; then
    echo -e "Generating project build system...\n"
    cmake $CMAKE_FLAGS -DCMAKE_C_COMPILER=$C_compiler -DCMAKE_CXX_COMPILER=$CXX_compiler -S $source_dir
  fi

  echo "Building with the following command:"
  echo -e "  $ cmake --build . --parallel $num_jobs\n"
  cmake --build . --parallel $num_jobs

  cd $source_dir

  echo -e "Done building configuration $folder\n"
  echo -e   "=---------------------------------------------------------------------------------=\n"

done

echo -e "Finished building all variants!"
