# CMake generated Testfile for 
# Source directory: /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCApp
# Build directory: /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCApp
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(build_output_qmcpack_exists "ls" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/bin/qmcpack_complex")
set_tests_properties(build_output_qmcpack_exists PROPERTIES  LABELS "unit;deterministic" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/unit_test.cmake;54;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCApp/CMakeLists.txt;61;add_test_target_in_output_location;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCApp/CMakeLists.txt;0;")
