# CMake generated Testfile for 
# Source directory: /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/Containers/Pools/tests
# Build directory: /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/Containers/Pools/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(deterministic-unit_test_containers_pools "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpiexec" "-n" "1" "--bind-to" "none" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/Containers/Pools/tests/test_containers_pools")
set_tests_properties(deterministic-unit_test_containers_pools PROPERTIES  ENVIRONMENT "OMP_NUM_THREADS=1" LABELS "quality_unknown;deterministic;unit" PROCESSORS "1" PROCESSOR_AFFINITY "TRUE" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/unit_test.cmake;8;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/Containers/Pools/tests/CMakeLists.txt;18;add_unit_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/Containers/Pools/tests/CMakeLists.txt;0;")
