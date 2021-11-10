#!/bin/bash

set -x

case "$1" in 

  # Configure qmcpack using cmake out-of-source builds 
  configure)
    
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
    
    case "${GH_JOBNAME}" in
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
      *"Clang12-NoMPI-Offload-Real"*)
        echo 'Configure for building OpenMP offload with clang12 on x86_64 target'
        cmake -GNinja \
              -DCMAKE_C_COMPILER=clang-12 \
              -DCMAKE_CXX_COMPILER=clang++-12 \
              -DQMC_MPI=0 \
              -DENABLE_OFFLOAD=ON \
              -DOFFLOAD_TARGET=x86_64-pc-linux-gnu \
              -DUSE_OBJECT_TARGET=ON \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
              ${GITHUB_WORKSPACE}
      ;;
      *"Clang14Dev-MPI-CUDA-AFQMC-Offload"*)
        echo "Configure for building with ENABLE_CUDA and AFQMC using OpenMP offload on x86_64 " \
              "with llvm development commit 01d59c0de822, need built-from-source OpenBLAS due to bug in rpm"
              # TODO: upgrade to llvm14 clang14 when available
        export OMPI_CC=/opt/llvm/01d59c0de822/bin/clang
        export OMPI_CXX=/opt/llvm/01d59c0de822/bin/clang++
        
        # Make current environment variables available to subsequent steps
        echo "OMPI_CC=/opt/llvm/01d59c0de822/bin/clang" >> $GITHUB_ENV
        echo "OMPI_CXX=/opt/llvm/01d59c0de822/bin/clang++" >> $GITHUB_ENV

        cmake -GNinja \
              -DCMAKE_C_COMPILER=/usr/lib64/openmpi/bin/mpicc \
              -DCMAKE_CXX_COMPILER=/usr/lib64/openmpi/bin/mpicxx \
              -DMPIEXEC_EXECUTABLE=/usr/lib64/openmpi/bin/mpirun \
              -DBUILD_AFQMC=ON \
              -DENABLE_CUDA=ON \
              -DENABLE_OFFLOAD=ON \
              -DUSE_OBJECT_TARGET=ON \
              -DCMAKE_PREFIX_PATH="/opt/OpenBLAS/0.3.18" \
              -DQMC_COMPLEX=$IS_COMPLEX \
              -DQMC_MIXED_PRECISION=$IS_MIXED_PRECISION \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
              ${GITHUB_WORKSPACE}
      ;;
      *"ROCm-Clang13-NoMPI-CUDA2HIP"*)
        echo 'Configure for building CUDA2HIP with clang compilers shipped with ROCM on AMD hardware'
        cmake -GNinja \
              -DCMAKE_C_COMPILER=/opt/rocm/llvm/bin/clang \
              -DCMAKE_CXX_COMPILER=/opt/rocm/llvm/bin/clang++ \
              -DQMC_MPI=0 \
              -DENABLE_CUDA=ON \
              -DQMC_CUDA2HIP=ON \
              -DCMAKE_PREFIX_PATH="/opt/OpenBLAS/0.3.18" \
              -DQMC_COMPLEX=$IS_COMPLEX \
              -DQMC_MIXED_PRECISION=$IS_MIXED_PRECISION \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
              ${GITHUB_WORKSPACE}
      ;;
      *"GCC8-MPI-CUDA-AFQMC"*)
        echo 'Configure for building with ENABLE_CUDA and AFQMC, need built-from-source OpenBLAS due to bug in rpm'
        cmake -GNinja \
              -DCMAKE_C_COMPILER=/usr/lib64/openmpi/bin/mpicc \
              -DCMAKE_CXX_COMPILER=/usr/lib64/openmpi/bin/mpicxx \
              -DMPIEXEC_EXECUTABLE=/usr/lib64/openmpi/bin/mpirun \
              -DBUILD_AFQMC=ON \
              -DENABLE_CUDA=ON \
              -DCMAKE_PREFIX_PATH="/opt/OpenBLAS/0.3.18" \
              -DQMC_COMPLEX=$IS_COMPLEX \
              -DQMC_MIXED_PRECISION=$IS_MIXED_PRECISION \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
              ${GITHUB_WORKSPACE}
      ;;
      *"GCC8-NoMPI-Legacy-CUDA"*)
        echo 'Configure for building with Legacy CUDA'
        cmake -GNinja \
              -DQMC_CUDA=1 \
              -DQMC_MPI=0 \
              -DQMC_COMPLEX=$IS_COMPLEX \
              -DQMC_MIXED_PRECISION=$IS_MIXED_PRECISION \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
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
    fi 
    
    if [[ "${GH_JOBNAME}" =~ (Clang12-NoMPI-Offload) ]]
    then
       echo "Adding /usr/lib/llvm-12/lib/ to LD_LIBRARY_PATH to enable libomptarget.so"
       export LD_LIBRARY_PATH=/usr/lib/llvm-12/lib/:${LD_LIBRARY_PATH}
       # Run only unit tests (reasonable for CI)
       TEST_LABEL="-L unit"
    fi

    if [[ "${GH_JOBNAME}" =~ (CUDA) ]]
    then
       export LD_LIBRARY_PATH=/usr/local/cuda/lib/:/usr/local/cuda/lib64/:${LD_LIBRARY_PATH}
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
       export LD_LIBRARY_PATH=/opt/llvm/01d59c0de822/lib:/usr/lib64/openmpi/lib/:${LD_LIBRARY_PATH}
    fi
    
    ctest --output-on-failure $TEST_LABEL
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

  *)
    echo " Invalid step" "$1"
    exit -1
    ;;
esac
