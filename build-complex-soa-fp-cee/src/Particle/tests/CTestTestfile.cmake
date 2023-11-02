# CMake generated Testfile for 
# Source directory: /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/Particle/tests
# Build directory: /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/Particle/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(deterministic-unit_test_particle "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpiexec" "-n" "1" "--bind-to" "none" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/Particle/tests/test_particle")
set_tests_properties(deterministic-unit_test_particle PROPERTIES  ENVIRONMENT "OMP_NUM_THREADS=1" LABELS "quality_unknown;deterministic;unit" PROCESSORS "1" PROCESSOR_AFFINITY "TRUE" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/unit_test.cmake;8;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/Particle/tests/CMakeLists.txt;29;add_unit_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/Particle/tests/CMakeLists.txt;0;")
add_test(deterministic-unit_test_distance_table "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpiexec" "-n" "1" "--bind-to" "none" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/Particle/tests/test_distance_table")
set_tests_properties(deterministic-unit_test_distance_table PROPERTIES  ENVIRONMENT "OMP_NUM_THREADS=1" LABELS "quality_unknown;deterministic;unit" PROCESSORS "1" PROCESSOR_AFFINITY "TRUE" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/unit_test.cmake;8;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/Particle/tests/CMakeLists.txt;43;add_unit_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/Particle/tests/CMakeLists.txt;0;")
