# CMake generated Testfile for 
# Source directory: /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/solids/InSe_a4
# Build directory: /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/tests/solids/InSe_a4
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(deterministic-InSe_a4_slab "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpiexec" "-n" "1" "--bind-to" "none" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/bin/qmcpack_complex" "InSe-S1.xml")
set_tests_properties(deterministic-InSe_a4_slab PROPERTIES  ENVIRONMENT "OMP_NUM_THREADS=1" FAIL_REGULAR_EXPRESSION "ERROR" LABELS "quality_unknown;deterministic;QMCPACK" PASS_REGULAR_EXPRESSION "QMCPACK execution completed successfully" PROCESSORS "1" PROCESSOR_AFFINITY "TRUE" WORKING_DIRECTORY "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/tests/solids/InSe_a4/deterministic-InSe_a4_slab" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;145;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;234;run_qmc_app_no_copy;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/solids/InSe_a4/CMakeLists.txt;4;run_qmc_app;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/solids/InSe_a4/CMakeLists.txt;0;")
