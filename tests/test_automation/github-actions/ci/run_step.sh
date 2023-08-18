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
    
    # Real or Complex configuration
    case "${GH_JOBNAME}" in
      *"Real"*)
        echo 'Configure for real build -DQMC_COMPLEX=0'
        IS_COMPLEX=0
      ;;
      *"Complex"*)
        echo 'Configure for complex build -DQMC_COMPLEX=1'
        IS_COMPLEX=1
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

    # Mixed or Non-Mixed (default, full) precision, used with GPU code
    case "${GH_JOBNAME}" in
      *"Mixed"*)
        echo 'Configure for mixed precision build -DQMC_MIXED_PRECISION=1'
        IS_MIXED_PRECISION=1
      ;; 
      *)
        IS_MIXED_PRECISION=0
      ;;
    esac

    if [[ "$CONTAINER_OS" =~ (centos) ]]
    then
      # use spack
      export PATH=/opt/rh/gcc-toolset-11/root/bin/:/opt/view:/opt/view/bin:$PATH
      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`which gcc|sed 's/bin\/gcc/lib64/g'`
      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/view/lib
      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/view/include
      export FFTW_HOME=/opt/view
      export LibXml2_ROOT=/opt/view
      export HDF5_ROOT=/opt/view
      export BOOST_ROOT=/opt/view


      # Make current environment variables available to subsequent steps
      echo "PATH=/opt/rh/gcc-toolset-11/root/bin/:/opt/view:/opt/view/bin:$PATH" >> $GITHUB_ENV
      echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`which gcc|sed 's/bin\/gcc/lib64/g'`" >> $GITHUB_ENV
      echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/view/lib" >> $GITHUB_ENV
      echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/view/include" >> $GITHUB_ENV
      echo "FFTW_HOME=/opt/view" >> $GITHUB_ENV
      echo "LibXml2_ROOT=/opt/view" >> $GITHUB_ENV
      echo "HDF5_ROOT=/opt/view" >> $GITHUB_ENV
      echo "BOOST_ROOT=/opt/view" >> $GITHUB_ENV
    fi
    
    case "${GH_JOBNAME}" in
      *"GCC9-NoMPI-Debug-"*|*"GCC11-NoMPI-Debug-"*)
        echo 'Configure for debug mode to capture asserts with gcc'
        cmake -GNinja \
              -DCMAKE_C_COMPILER=gcc \
              -DCMAKE_CXX_COMPILER=g++ \
              -DQMC_MPI=0 \
              -DCMAKE_BUILD_TYPE=Debug \
              ${GITHUB_WORKSPACE}
      ;;
      *"GCC9-NoMPI-NoOMP-"*|*"GCC11-NoMPI-NoOMP-"*)
        echo 'Configure for disabling OpenMP with QMC_OMP=0'
        cmake -GNinja \
              -DCMAKE_C_COMPILER=gcc \
              -DCMAKE_CXX_COMPILER=g++ \
              -DQMC_MPI=0 \
              -DQMC_OMP=0 \
              -DQMC_COMPLEX=$IS_COMPLEX \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
              ${GITHUB_WORKSPACE}
      ;;
      *"GCC9-NoMPI-Sandbox-"*|*"GCC11-NoMPI-Sandbox-"*)
        echo 'Configure for enabling sandbox (minimal) only option with gcc'
        cmake -GNinja \
              -DCMAKE_C_COMPILER=gcc \
              -DCMAKE_CXX_COMPILER=g++ \
              -DQMC_MPI=0 \
              -DQMC_BUILD_SANDBOX_ONLY=ON \
              -DQMC_COMPLEX=$IS_COMPLEX \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
              ${GITHUB_WORKSPACE}
      ;;
      *"GCC9-MPI-Gcov-"*)
        echo 'Configure for code coverage with gcc and gcovr -DENABLE_GCOV=TRUE and upload reports to Codecov'
        cmake -GNinja \
              -DMPI_C_COMPILER=mpicc \
              -DMPI_CXX_COMPILER=mpicxx \
              -DENABLE_GCOV=TRUE \
              -DQMC_COMPLEX=$IS_COMPLEX \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
              ${GITHUB_WORKSPACE}
      ;;
      *"GCC11-NoMPI-Werror-"*)
        echo 'Configure for building with gcc -Werror flag enabled'
        cmake -GNinja \
              -DCMAKE_C_COMPILER=gcc \
              -DCMAKE_CXX_COMPILER=g++ \
              -DQMC_MPI=0 \
              -DCMAKE_CXX_FLAGS=-Werror \
              -DQMC_COMPLEX=$IS_COMPLEX \
              -DQMC_MIXED_PRECISION=$IS_MIXED_PRECISION \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
              ${GITHUB_WORKSPACE}
      ;;
      *"San-"*) # Sanitize with clang compilers
        cmake -GNinja \
              -DCMAKE_C_COMPILER=clang \
              -DCMAKE_CXX_COMPILER=clang++ \
              -DQMC_MPI=0 \
              -DENABLE_SANITIZER=$IS_SANITIZER \
              -DQMC_COMPLEX=$IS_COMPLEX \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
              ${GITHUB_WORKSPACE}
      ;;
      *"Clang16-NoMPI-Offload-Real"*)
        echo 'Configure for building OpenMP offload with clang16 on x86_64 target'
        cmake -GNinja \
              -DCMAKE_C_COMPILER=clang-16 \
              -DCMAKE_CXX_COMPILER=clang++-16 \
              -DQMC_MPI=0 \
              -DENABLE_OFFLOAD=ON \
              -DOFFLOAD_TARGET=x86_64-pc-linux-gnu \
              -DUSE_OBJECT_TARGET=ON \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
              ${GITHUB_WORKSPACE}
      ;;
      *"Clang15-MPI-CUDA-AFQMC-Offload"*)
        echo "Configure for building with ENABLE_CUDA and AFQMC using OpenMP offload on x86_64 " \
              "with latest llvm, need built-from-source OpenBLAS due to bug in rpm"

        # todo: update to llvm 15 release, currently using release candidate
        export OMPI_CC=/opt/llvm/15.0.0/bin/clang
        export OMPI_CXX=/opt/llvm/15.0.0/bin/clang++
        
        # Make current environment variables available to subsequent steps
        echo "OMPI_CC=/opt/llvm/15.0.0/bin/clang" >> $GITHUB_ENV
        echo "OMPI_CXX=/opt/llvm/15.0.0/bin/clang++" >> $GITHUB_ENV

        # Confirm that cuda 11.2 gets picked up by the compiler
        /opt/llvm/15.0.0/bin/clang++ -v

        cmake -GNinja \
              -DCMAKE_C_COMPILER=/usr/lib64/openmpi/bin/mpicc \
              -DCMAKE_CXX_COMPILER=/usr/lib64/openmpi/bin/mpicxx \
              -DMPIEXEC_EXECUTABLE=/usr/lib64/openmpi/bin/mpirun \
              -DBUILD_AFQMC=ON \
              -DENABLE_CUDA=ON \
              -DQMC_GPU_ARCHS=sm_70 \
              -DENABLE_OFFLOAD=ON \
              -DCMAKE_PREFIX_PATH="/opt/OpenBLAS/0.3.18" \
              -DQMC_COMPLEX=$IS_COMPLEX \
              -DQMC_MIXED_PRECISION=$IS_MIXED_PRECISION \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
              -DQMC_DATA=$QMC_DATA_DIR \
              ${GITHUB_WORKSPACE}
      ;;
      *"Intel21-MPI-CUDA-AFQMC"*)
        echo "Configure for building with ENABLE_CUDA and AFQMC  " \
             "with Intel classic compiler in OneAPI 2021 (to be deprecated in 2023), " \
             "need built-from-source OpenBLAS due to bug in rpm"
        
        source /opt/intel/oneapi/setvars.sh

        export OMPI_CC=/opt/intel/oneapi/compiler/2023.0.0/linux/bin/intel64/icc
        export OMPI_CXX=/opt/intel/oneapi/compiler/2023.0.0/linux/bin/intel64/icpc
        
        # Make current environment variables available to subsequent steps
        echo "OMPI_CC=/opt/intel/oneapi/compiler/2023.0.0/linux/bin/intel64/icc" >> $GITHUB_ENV
        echo "OMPI_CXX=/opt/intel/oneapi/compiler/2023.0.0/linux/bin/intel64/icpc" >> $GITHUB_ENV

        cmake -GNinja \
              -DCMAKE_C_COMPILER=/usr/lib64/openmpi/bin/mpicc \
              -DCMAKE_CXX_COMPILER=/usr/lib64/openmpi/bin/mpicxx \
              -DMPIEXEC_EXECUTABLE=/usr/lib64/openmpi/bin/mpirun \
              -DCMAKE_C_FLAGS="-diag-disable=10441" \
              -DCMAKE_CXX_FLAGS="-diag-disable=10441" \
              -DCMAKE_CUDA_FLAGS="-diag-disable=10441" \
              -DBUILD_AFQMC=ON \
              -DENABLE_CUDA=ON \
              -DQMC_GPU_ARCHS=sm_70 \
              -DCMAKE_PREFIX_PATH="/opt/OpenBLAS/0.3.18" \
              -DQMC_COMPLEX=$IS_COMPLEX \
              -DQMC_MIXED_PRECISION=$IS_MIXED_PRECISION \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
              -DQMC_DATA=$QMC_DATA_DIR \
              ${GITHUB_WORKSPACE}
      ;;
      *"GCC9-MPI-CUDA-AFQMC"*)
        echo 'Configure for building with ENABLE_CUDA and AFQMC, need built-from-source OpenBLAS due to bug in rpm'
        cmake -GNinja \
              -DCMAKE_C_COMPILER=/usr/lib64/openmpi/bin/mpicc \
              -DCMAKE_CXX_COMPILER=/usr/lib64/openmpi/bin/mpicxx \
              -DMPIEXEC_EXECUTABLE=/usr/lib64/openmpi/bin/mpirun \
              -DBUILD_AFQMC=ON \
              -DENABLE_CUDA=ON \
              -DQMC_GPU_ARCHS=sm_70 \
              -DCMAKE_PREFIX_PATH="/opt/OpenBLAS/0.3.18" \
              -DQMC_COMPLEX=$IS_COMPLEX \
              -DQMC_MIXED_PRECISION=$IS_MIXED_PRECISION \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
              -DQMC_DATA=$QMC_DATA_DIR \
              ${GITHUB_WORKSPACE}
      ;;
      *"GCC9-NoMPI-MKL-"*)
        echo 'Configure for building with GCC and Intel MKL'

        source /opt/intel2020/mkl/bin/mklvars.sh intel64

        cmake -GNinja \
              -DBLA_VENDOR=Intel10_64lp \
              -DQMC_MPI=0 \
              -DQMC_COMPLEX=$IS_COMPLEX \
              -DQMC_MIXED_PRECISION=$IS_MIXED_PRECISION \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
              -DQMC_DATA=$QMC_DATA_DIR \
              ${GITHUB_WORKSPACE}
      ;;
      *"macOS-GCC11-NoMPI-Real"*)
        echo 'Configure for building on macOS using gcc11'
        cmake -GNinja \
              -DCMAKE_C_COMPILER=gcc-11 \
              -DCMAKE_CXX_COMPILER=g++-11 \
              -DQMC_MPI=0 \
              -DQMC_COMPLEX=$IS_COMPLEX \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
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
    TEST_LABEL="-L deterministic"
    
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
       echo "Adding /usr/lib/llvm-12/lib/ to LD_LIBRARY_PATH to enable libomptarget.so"
       export LD_LIBRARY_PATH=/usr/lib/llvm-12/lib/:${LD_LIBRARY_PATH}
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
       export LD_LIBRARY_PATH=/opt/llvm/15.0.0/lib:/usr/lib64/openmpi/lib:${LD_LIBRARY_PATH}
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
    # Default for Linux GitHub Action runners
    CTEST_JOBS="2"
    # Default for macOS GitHub Action runners
    if [[ "${GH_OS}" =~ (macOS) ]]
    then
      CTEST_JOBS="3"
    fi

    if [[ "$HOST_NAME" =~ (sulfur) || "$HOST_NAME" =~ (nitrogen) ]]
    then
      CTEST_JOBS="32"
    fi
    
    ctest --output-on-failure $TEST_LABEL -j $CTEST_JOBS
    ;;
  
  # Generate coverage reports
  coverage)
    cd ${GITHUB_WORKSPACE}/../qmcpack-build
    # filter unreachable branches with gcovr
    # see https://gcovr.com/en/stable/faq.html#why-does-c-code-have-so-many-uncovered-branches
    gcovr --exclude-unreachable-branches --exclude-throw-branches --root=${GITHUB_WORKSPACE}/.. --xml-pretty -o coverage.xml
    du -hs coverage.xml
    cat coverage.xml
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
