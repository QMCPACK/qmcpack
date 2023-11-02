# CMake generated Testfile for 
# Source directory: /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/ppconvert/test/Oxygen
# Build directory: /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools/ppconvert/test/Oxygen
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(ppconvert_runs "bash" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools/ppconvert/test/Oxygen/O_test.sh" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/bin/ppconvert")
set_tests_properties(ppconvert_runs PROPERTIES  LABELS "deterministic" WORKING_DIRECTORY "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools/ppconvert/test/Oxygen" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/ppconvert/test/Oxygen/CMakeLists.txt;6;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/ppconvert/test/Oxygen/CMakeLists.txt;0;")
