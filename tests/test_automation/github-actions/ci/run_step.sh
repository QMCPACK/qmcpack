#!/bin/bash

case "$1" in 

  # Configure qmcpack using cmake out-of-source builds 
  configure)
    
    cd ${GITHUB_WORKSPACE}/..
    mkdir qmcpack-build
    cd qmcpack-build
    
    case "${GH_JOBNAME}" in
      *"real"*)
        echo 'Configure for real build -DQMC_COMPLEX=0'
        IS_COMPLEX=0
      ;;
      *"complex"*)
        echo 'Configure for complex build -DQMC_COMPLEX=1'
        IS_COMPLEX=1
      ;; 
    esac
    
    case "${GH_JOBNAME}" in
      # Sanitize with clang compilers
      *"asan"*)
        echo 'Configure for address sanitizer asan including lsan (leaks)'
        CC=clang CXX=clang++ \
        cmake -GNinja -DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_SANITIZER=asan \
                      -DQMC_MPI=0 \
                      -DQMC_COMPLEX=$IS_COMPLEX \
                      ${GITHUB_WORKSPACE}
      ;;
      *"ubsan"*)
        echo 'Configure for undefined behavior sanitizer ubsan'
        CC=clang CXX=clang++ \
        cmake -GNinja -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx \
                      -DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_SANITIZER=ubsan \
                      -DQMC_COMPLEX=$IS_COMPLEX \
                      ${GITHUB_WORKSPACE}
      ;;
      *"tsan"*)
        echo 'Configure for thread sanitizer tsan'
        CC=clang CXX=clang++ \
        cmake -GNinja -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx \
                      -DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_SANITIZER=tsan \
                      -DQMC_COMPLEX=$IS_COMPLEX \
                      ${GITHUB_WORKSPACE}
      ;;
      *"msan"*)
        echo 'Configure for (uninitialized) memory sanitizer msan'
        CC=clang CXX=clang++ \
        cmake -GNinja -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx \
                      -DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_SANITIZER=msan \
                      -DQMC_COMPLEX=$IS_COMPLEX \
                      ${GITHUB_WORKSPACE}
      ;;
      *"coverage"*)
        echo 'Configure for code coverage with gcc and gcovr'
        cmake -GNinja -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx \
                      -DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_GCOV=TRUE \
                      -DQMC_COMPLEX=$IS_COMPLEX \
                      ${GITHUB_WORKSPACE}
      ;;
      *"clang-latest-openmp-offload"*)
        echo 'Configure for building OpenMP offload with clang-12 on x86_64'
        cmake -GNinja -DCMAKE_C_COMPILER=clang-12 -DCMAKE_CXX_COMPILER=clang++-12 \
                      -DENABLE_OFFLOAD=ON -DOFFLOAD_TARGET=x86_64-pc-linux-gnu \
                      -DUSE_OBJECT_TARGET=ON -DQMC_MPI=0 \
                      ${GITHUB_WORKSPACE}
      ;;
      # Configure with default compilers
      *)
        echo 'Configure for default system compilers and options'
        cmake -GNinja -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx \
                      -DQMC_COMPLEX=$IS_COMPLEX \
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
    cd ${GITHUB_WORKSPACE}/../qmcpack-build
    
    # Enable oversubscription in OpenMPI
    if [[ "${GH_JOBNAME}" =~ (openmpi) ]]
    then
      echo "Enabling OpenMPI oversubscription"
      export OMPI_MCA_rmaps_base_oversubscribe=1
      export OMPI_MCA_hwloc_base_binding_policy=none
    fi 
    
    # Run only deterministic tests (reasonable for CI) by default
    TEST_LABEL="-L deterministic"
    
    # Enable ASAN_OPTION=suppression=suppresion_file
    if [[ "${GH_JOBNAME}" =~ (asan) ]]
    then
      TEST_LABEL="-L unit -LE noasan"
    fi
    
    if [[ "${GH_JOBNAME}" =~ (clang-latest-openmp-offload) ]]
    then
       echo "Adding /usr/lib/llvm-12/lib/ to LD_LIBRARY_PATH to enable libomptarget.so"
       export LD_LIBRARY_PATH=/usr/lib/llvm-12/lib/:${LD_LIBRARY_PATH}
       # Run only unit tests (reasonable for CI using openmp-offload)
       TEST_LABEL="-L unit"
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
