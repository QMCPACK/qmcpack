# CMake generated Testfile for 
# Source directory: /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/ppconvert/src
# Build directory: /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools/ppconvert/src
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(build_output_ppconvert_exists "ls" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/bin/ppconvert")
set_tests_properties(build_output_ppconvert_exists PROPERTIES  LABELS "unit;deterministic" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/unit_test.cmake;54;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/ppconvert/src/CMakeLists.txt;13;add_test_target_in_output_location;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/ppconvert/src/CMakeLists.txt;0;")
subdirs("common")
