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
        cmake -GNinja -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx \
                      -DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_SANITIZER=asan \
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
    
    # Enable ASAN_OPTION=suppression=suppresion_file
    if [[ "${GH_JOBNAME}" =~ (asan) ]]
    then
      echo "Enabling LSAN suppression file config/sanitizers/lsan.supp"
      export LSAN_OPTIONS=suppression=${GITHUB_WORKSPACE}/config/sanitizers/lsan.supp	
    fi
    
    # Run only deterministic tests (reasonable for CI)
    ctest -L deterministic
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
