#!/bin/bash

set -x
HOST_NAME=$(hostname -s)

case "$1" in 

  # Configure qmcpack using cmake out-of-source builds 
  configure)
    
    if [[ "$HOST_NAME" =~ (sulfur) || "$HOST_NAME" =~ (nitrogen) ]]
    then
      echo "Use recent cmake v3.21.3"
      export PATH=/opt/cmake-3.21.3-linux-x86_64/bin:$PATH
      # Make current environment variables available to subsequent steps, 
      # e.g. ctest
      echo "PATH=$PATH" >> $GITHUB_ENV
    fi
    
    if [ -d ${GITHUB_WORKSPACE}/../qmcpack-build ]
    then
      echo "Found existing out-of-source build directory ${GITHUB_WORKSPACE}/../qmcpack-build, removing"
      rm -fr ${GITHUB_WORKSPACE}/../qmcpack-build
    fi
    
    echo "Creating new out-of-source build directory ${GITHUB_WORKSPACE}/../qmcpack-build"
    cd ${GITHUB_WORKSPACE}/..
    mkdir qmcpack-build
    cd qmcpack-build
    
    # MPI or not configuration
    if [[ "${GH_JOBNAME}" =~ (-NoMPI) ]] ; then
      echo 'Configure for build without MPI'
      CMAKE_OPTIONS="-DQMC_MPI=OFF"
    else
      echo 'Configure for MPI build'
      CMAKE_OPTIONS="-DQMC_MPI=ON"
    fi

    # OpenMP or not configuration
    if [[ "${GH_JOBNAME}" =~ (-NoOMP) ]] ; then
      echo 'Configure for build without OpenMP'
      CMAKE_OPTIONS="$CMAKE_OPTIONS -DQMC_OMP=OFF"
    else
      echo 'Configure for OpenMP build'
      CMAKE_OPTIONS="$CMAKE_OPTIONS -DQMC_OMP=ON"
    fi

    # build type configuration
    if [[ "${GH_JOBNAME}" =~ (-Debug) ]] ; then
      echo 'Configure for Debug build'
      CMAKE_OPTIONS="$CMAKE_OPTIONS -DCMAKE_BUILD_TYPE=Debug"
    else
      echo 'Configure for RelWithDebInfo build'
      CMAKE_OPTIONS="$CMAKE_OPTIONS -DCMAKE_BUILD_TYPE=RelWithDebInfo"
    fi

    # Real or Complex configuration
    case "${GH_JOBNAME}" in
      *"Real"*)
        echo 'Configure for real build -DQMC_COMPLEX=OFF'
        CMAKE_OPTIONS="$CMAKE_OPTIONS -DQMC_COMPLEX=OFF"
      ;;
      *"Complex"*)
        echo 'Configure for complex build -DQMC_COMPLEX=ON'
        CMAKE_OPTIONS="$CMAKE_OPTIONS -DQMC_COMPLEX=ON"
      ;; 
    esac

    if [[ "${GH_JOBNAME}" =~ (-CUDA) ]]
    then
      if [[ "${GH_JOBNAME}" =~ (-Offload) ]]
      then
        echo "Set PATH to cuda-11.2 to be associated with the C and C++ compilers"
        export PATH=/usr/local/cuda-11.2/bin:$PATH
        echo "Set CUDACXX CMake environment variable to nvcc cuda 11.2 location due to a regression bug in 11.6"
        export CUDACXX=/usr/local/cuda-11.2/bin/nvcc

        # Make current environment variables available to subsequent steps
        echo "PATH=$PATH" >> $GITHUB_ENV
        echo "CUDACXX=/usr/local/cuda-11.2/bin/nvcc" >> $GITHUB_ENV

      else
        echo "Set CUDACXX CMake environment variable to nvcc 11.8"
        export CUDACXX=/usr/local/cuda-11.8/bin/nvcc

        # Make current environment variables available to subsequent steps
        echo "CUDACXX=/usr/local/cuda-11.8/bin/nvcc" >> $GITHUB_ENV
      fi
    fi 

    # Sanitizer
    case "${GH_JOBNAME}" in
      *"ASan"*)
        echo 'Configure for address sanitizer including leak sanitizer (lsan) -DENABLE_SANITIZER=asan'
        IS_SANITIZER=asan
      ;;
      *"UBSan"*)
        echo 'Configure for undefined behavior sanitizer -DENABLE_SANITIZER=ubsan'
        IS_SANITIZER=ubsan
      ;; 
      *"TSan"*)
        echo 'Configure for thread sanitizer -DENABLE_SANITIZER=tsan'
        IS_SANITIZER=tsan
      ;;
    esac

    # Mixed or Non-Mixed (default, full) precision
    if [[ "${GH_JOBNAME}" =~ (-Mixed) ]] ; then
      echo 'Configure for mixed precision build -DQMC_MIXED_PRECISION=ON'
      CMAKE_OPTIONS="$CMAKE_OPTIONS -DQMC_MIXED_PRECISION=ON"
    else
      CMAKE_OPTIONS="$CMAKE_OPTIONS -DQMC_MIXED_PRECISION=OFF"
    fi

    if [[ "$CONTAINER_OS" =~ (centos) ]]
    then
       module avail
       module load mpi/openmpi-x86_64
       module list
    fi
    
    case "${GH_JOBNAME}" in
      *"macOS-GCC14"*"-Real"*)
        echo 'Configure for building on macOS using gcc14'
        cmake -GNinja $CMAKE_OPTIONS \
              -DCMAKE_C_COMPILER=gcc-14 \
              -DCMAKE_CXX_COMPILER=g++-14 \
              -DCMAKE_EXE_LINKER_FLAGS="-Wl,-ld_classic" \
              ${GITHUB_WORKSPACE}
      ;;
      *"GCC9"*"-CUDA-AFQMC"*)
        echo 'Configure for building with CUDA and AFQMC, need built-from-source OpenBLAS due to bug in rpm'
        cmake -GNinja $CMAKE_OPTIONS \
              -DCMAKE_C_COMPILER=/usr/lib64/openmpi/bin/mpicc \
              -DCMAKE_CXX_COMPILER=/usr/lib64/openmpi/bin/mpicxx \
              -DMPIEXEC_EXECUTABLE=/usr/lib64/openmpi/bin/mpirun \
              -DBUILD_AFQMC=ON \
              -DQMC_GPU=cuda \
              -DQMC_GPU_ARCHS=sm_70 \
              -DCMAKE_PREFIX_PATH="/opt/OpenBLAS/0.3.18" \
              -DQMC_DATA=$QMC_DATA_DIR \
              ${GITHUB_WORKSPACE}
      ;;
      *"GCC9"*"-MKL-"*)
        echo 'Configure for building with GCC and Intel MKL'

        source /opt/intel2020/mkl/bin/mklvars.sh intel64

        cmake -GNinja $CMAKE_OPTIONS \
              -DBLA_VENDOR=Intel10_64lp \
              -DQMC_DATA=$QMC_DATA_DIR \
              ${GITHUB_WORKSPACE}
      ;;
      *"GCC"*"-Sandbox"*)
        echo 'Configure for enabling sandbox (minimal) only option with gcc'
        cmake -GNinja $CMAKE_OPTIONS \
              -DCMAKE_C_COMPILER=gcc \
              -DCMAKE_CXX_COMPILER=g++ \
              -DQMC_BUILD_SANDBOX_ONLY=ON \
              ${GITHUB_WORKSPACE}
      ;;
      *"GCC"*"-Gcov"*)
        echo 'Configure for code coverage with gcc and gcovr -DENABLE_GCOV=TRUE and upload reports to Codecov'
        cmake -GNinja $CMAKE_OPTIONS \
              -DMPI_C_COMPILER=mpicc \
              -DMPI_CXX_COMPILER=mpicxx \
              -DENABLE_GCOV=TRUE \
              -DENABLE_PYCOV=TRUE \
              ${GITHUB_WORKSPACE}
      ;;
      *"GCC"*"-Werror"*)
        echo 'Configure for building with gcc -Werror flag enabled'
        cmake -GNinja $CMAKE_OPTIONS \
              -DCMAKE_C_COMPILER=gcc \
              -DCMAKE_CXX_COMPILER=g++ \
              -DCMAKE_CXX_FLAGS=-Werror \
              ${GITHUB_WORKSPACE}
      ;;
      *"GCC"*)
        echo 'Configure for disabling OpenMP with QMC_OMP=0'
        cmake -GNinja $CMAKE_OPTIONS \
              -DCMAKE_C_COMPILER=gcc \
              -DCMAKE_CXX_COMPILER=g++ \
              ${GITHUB_WORKSPACE}
      ;;
      *"Clang"*"San"*) # Sanitize with clang compilers
        cmake -GNinja $CMAKE_OPTIONS \
              -DCMAKE_C_COMPILER=clang \
              -DCMAKE_CXX_COMPILER=clang++ \
              -DENABLE_SANITIZER=$IS_SANITIZER \
              ${GITHUB_WORKSPACE}
      ;;
      *"Clang16"*"-Offload"*)
        echo 'Configure for building OpenMP offload with clang16 on x86_64 target'
        cmake -GNinja $CMAKE_OPTIONS \
              -DCMAKE_C_COMPILER=clang-16 \
              -DCMAKE_CXX_COMPILER=clang++-16 \
              -DQMC_GPU=openmp \
              -DOFFLOAD_TARGET=x86_64-pc-linux-gnu \
              -DUSE_OBJECT_TARGET=ON \
              ${GITHUB_WORKSPACE}
      ;;
      *"Clang15"*"-CUDA-AFQMC-Offload"*)
        echo "Configure for building with CUDA and AFQMC using OpenMP offload on x86_64 " \
              "with latest llvm, need built-from-source OpenBLAS due to bug in rpm"

        # todo: update to llvm 15 release, currently using release candidate
        export OMPI_CC=/opt/llvm/15.0.0/bin/clang
        export OMPI_CXX=/opt/llvm/15.0.0/bin/clang++
        
        # Make current environment variables available to subsequent steps
        echo "OMPI_CC=/opt/llvm/15.0.0/bin/clang" >> $GITHUB_ENV
        echo "OMPI_CXX=/opt/llvm/15.0.0/bin/clang++" >> $GITHUB_ENV

        # Confirm that cuda 11.2 gets picked up by the compiler
        /opt/llvm/15.0.0/bin/clang++ -v

        cmake -GNinja $CMAKE_OPTIONS \
              -DCMAKE_C_COMPILER=/usr/lib64/openmpi/bin/mpicc \
              -DCMAKE_CXX_COMPILER=/usr/lib64/openmpi/bin/mpicxx \
              -DMPIEXEC_EXECUTABLE=/usr/lib64/openmpi/bin/mpirun \
              -DBUILD_AFQMC=ON \
              -DQMC_GPU="cuda;openmp" \
              -DQMC_GPU_ARCHS=sm_70 \
              -DCMAKE_PREFIX_PATH="/opt/OpenBLAS/0.3.18" \
              -DQMC_DATA=$QMC_DATA_DIR \
              ${GITHUB_WORKSPACE}
      ;;
      *"Intel21"*"-CUDA-AFQMC"*)
        echo "Configure for building with CUDA and AFQMC  " \
             "with Intel classic compiler in OneAPI 2021 (to be deprecated in 2023), " \
             "need built-from-source OpenBLAS due to bug in rpm"
        
        source /opt/intel/oneapi/setvars.sh

        export OMPI_CC=/opt/intel/oneapi/compiler/2023.0.0/linux/bin/intel64/icc
        export OMPI_CXX=/opt/intel/oneapi/compiler/2023.0.0/linux/bin/intel64/icpc
        
        # Make current environment variables available to subsequent steps
        echo "OMPI_CC=/opt/intel/oneapi/compiler/2023.0.0/linux/bin/intel64/icc" >> $GITHUB_ENV
        echo "OMPI_CXX=/opt/intel/oneapi/compiler/2023.0.0/linux/bin/intel64/icpc" >> $GITHUB_ENV

        cmake -GNinja $CMAKE_OPTIONS \
              -DCMAKE_C_COMPILER=/usr/lib64/openmpi/bin/mpicc \
              -DCMAKE_CXX_COMPILER=/usr/lib64/openmpi/bin/mpicxx \
              -DMPIEXEC_EXECUTABLE=/usr/lib64/openmpi/bin/mpirun \
              -DCMAKE_C_FLAGS="-diag-disable=10441" \
              -DCMAKE_CXX_FLAGS="-diag-disable=10441" \
              -DCMAKE_CUDA_FLAGS="-diag-disable=10441" \
              -DBUILD_AFQMC=ON \
              -DQMC_GPU=cuda \
              -DQMC_GPU_ARCHS=sm_70 \
              -DCMAKE_PREFIX_PATH="/opt/OpenBLAS/0.3.18" \
              -DQMC_DATA=$QMC_DATA_DIR \
              ${GITHUB_WORKSPACE}
      ;;
    esac
    ;;

  # Build using ninja (~ 25 minutes on GitHub-hosted runner)
  build)
    # CUDA toolchain can be used implicitly by the compiler. Double check the location.
    if [[ "${GH_JOBNAME}" =~ (CUDA) ]]
    then
      which nvcc
    fi

    cd ${GITHUB_WORKSPACE}/../qmcpack-build
    ninja
    ;;

  # Run deterministic tests
  test)
    
    # Run only deterministic tests (reasonable for CI) by default
    case "${GH_JOBNAME}" in
      *"macOS-GCC14"*"-Real"*)
        TEST_LABEL="-L deterministic -E deterministic-unit_test_estimators"
        # estimator test bus error on mac only
      ;;
      *)  
        TEST_LABEL="-L deterministic"
      ;;  
    esac  

    cd ${GITHUB_WORKSPACE}/../qmcpack-build
    
    # Enable oversubscription in OpenMPI
    if [[ "${GH_JOBNAME}" =~ (-MPI-) ]]
    then
      echo "Enabling OpenMPI oversubscription"
      export OMPI_MCA_rmaps_base_oversubscribe=1
      export OMPI_MCA_hwloc_base_binding_policy=none
      
      if [[ "$HOST_NAME" =~ (sulfur) || "$HOST_NAME" =~ (nitrogen) ]]
      then
        echo "Set the management layer to ucx"
        export OMPI_MCA_pml=ucx
      fi
    fi 
    
    if [[ "${GH_JOBNAME}" =~ (Clang16-NoMPI-Offload) ]]
    then
       export KMP_TEAMS_THREAD_LIMIT=1
       # Run only unit tests (reasonable for CI)
       TEST_LABEL="-L unit"
    fi

    if [[ "${GH_JOBNAME}" =~ (CUDA) ]]
    then
      if [[ "${GH_JOBNAME}" =~ (-Offload) ]]
      then
        export LD_LIBRARY_PATH=/usr/local/cuda-11.2/lib64:${LD_LIBRARY_PATH}
      else
        export LD_LIBRARY_PATH=/usr/local/cuda-11.8/lib64:${LD_LIBRARY_PATH}
      fi
    fi

    if [[ "${GH_JOBNAME}" =~ (AFQMC) ]]
    then
       # Avoid polluting the stderr output with libfabric error message
       export OMPI_MCA_btl=self
    fi
    
    if [[ "${GH_JOBNAME}" =~ (Offload) ]]
    then
       # Clang helper threads used by target nowait is very broken. Disable this feature
       export LIBOMP_USE_HIDDEN_HELPER_TASK=0
    fi

    if [[ "${GH_JOBNAME}" =~ (AFQMC-Offload) ]]
    then
       export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:${LD_LIBRARY_PATH}
    fi

    if [[ "${GH_JOBNAME}" =~ (Intel21) ]]
    then
       source /opt/intel/oneapi/setvars.sh
    fi

    if [[ "${GH_JOBNAME}" =~ (MKL) ]]
    then 
       source /opt/intel2020/mkl/bin/mklvars.sh intel64
    fi

    # Add ctest concurrent parallel jobs 
    # 4 for Linux and macOS GitHub Actions free runners
    ctest --output-on-failure $TEST_LABEL -j 4
    ;;
  
  # Generate coverage reports
  coverage)
    cd ${GITHUB_WORKSPACE}/../qmcpack-build
    # filter unreachable branches with gcovr
    # see https://gcovr.com/en/stable/faq.html#why-does-c-code-have-so-many-uncovered-branches
    gcovr --exclude-unreachable-branches --exclude-throw-branches --root=${GITHUB_WORKSPACE}/.. --xml-pretty -o coverage.xml
    du -hs coverage.xml
    #cat coverage.xml
    python3-coverage combine nexus/nexus/tests/.coverage*
    du -hs .coverage
    python3-coverage report
    python3-coverage xml -o python_coverage.xml
    du -hs python_coverage.xml
    # Debug only: overwrite gcov file with python coverage to test codecov.io
    #mv python_coverage.xml coverage.xml
    ;;
  
  # Install the library (not triggered at the moment)
  install)
    cd ${GITHUB_WORKSPACE}/../qmcpack-build
    ninja install
    ;;

  rebase)
    source external_codes/github_actions/auto-rebase.sh
    ;;
  
  pull-rebase)
    source external_codes/github_actions/trigger-rebase-on-push.sh
    ;;

  *)
    echo " Invalid step" "$1"
    exit -1
    ;;
esac
