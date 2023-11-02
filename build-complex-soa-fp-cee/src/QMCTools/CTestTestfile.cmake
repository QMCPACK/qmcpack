# CMake generated Testfile for 
# Source directory: /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools
# Build directory: /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(build_output_convert4qmc_exists "ls" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/bin/convert4qmc")
set_tests_properties(build_output_convert4qmc_exists PROPERTIES  LABELS "unit;deterministic" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/unit_test.cmake;54;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/CMakeLists.txt;66;add_test_target_in_output_location;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/CMakeLists.txt;0;")
add_test(build_output_qmc-extract-eshdf-kvectors_exists "ls" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/bin/qmc-extract-eshdf-kvectors")
set_tests_properties(build_output_qmc-extract-eshdf-kvectors_exists PROPERTIES  LABELS "unit;deterministic" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/unit_test.cmake;54;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/CMakeLists.txt;66;add_test_target_in_output_location;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/CMakeLists.txt;0;")
add_test(build_output_qmc-get-supercell_exists "ls" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/bin/qmc-get-supercell")
set_tests_properties(build_output_qmc-get-supercell_exists PROPERTIES  LABELS "unit;deterministic" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/unit_test.cmake;54;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/CMakeLists.txt;66;add_test_target_in_output_location;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/CMakeLists.txt;0;")
add_test(build_output_qmc-check-affinity_exists "ls" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/bin/qmc-check-affinity")
set_tests_properties(build_output_qmc-check-affinity_exists PROPERTIES  LABELS "unit;deterministic" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/unit_test.cmake;54;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/CMakeLists.txt;66;add_test_target_in_output_location;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/CMakeLists.txt;0;")
add_test(build_output_convertpw4qmc_exists "ls" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/bin/convertpw4qmc")
set_tests_properties(build_output_convertpw4qmc_exists PROPERTIES  LABELS "unit;deterministic" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/unit_test.cmake;54;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/CMakeLists.txt;66;add_test_target_in_output_location;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/CMakeLists.txt;0;")
add_test(build_output_qmcfinitesize_exists "ls" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/bin/qmcfinitesize")
set_tests_properties(build_output_qmcfinitesize_exists PROPERTIES  LABELS "unit;deterministic" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/unit_test.cmake;54;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/CMakeLists.txt;66;add_test_target_in_output_location;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/CMakeLists.txt;0;")
subdirs("ppconvert")
subdirs("tests")
