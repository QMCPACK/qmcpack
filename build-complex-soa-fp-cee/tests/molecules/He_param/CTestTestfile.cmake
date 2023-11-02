# CMake generated Testfile for 
# Source directory: /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param
# Build directory: /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/tests/molecules/He_param
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(short-He_param_grad-1-16 "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpiexec" "-n" "1" "--bind-to" "none" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/bin/qmcpack_complex" "He_param_grad.xml")
set_tests_properties(short-He_param_grad-1-16 PROPERTIES  ENVIRONMENT "OMP_NUM_THREADS=16" FAIL_REGULAR_EXPRESSION "ERROR" LABELS "quality_unknown;QMCPACK;QMCPACK" PASS_REGULAR_EXPRESSION "QMCPACK execution completed successfully" PROCESSORS "16" PROCESSOR_AFFINITY "TRUE" WORKING_DIRECTORY "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/tests/molecules/He_param/short-He_param_grad-1-16" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;145;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;234;run_qmc_app_no_copy;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;497;run_qmc_app;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;8;qmc_run_and_check_custom_scalar;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;0;")
add_test(short-He_param_grad-1-16-jud_0 "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/anaconda3-2021.05-hzr5w5v54sh3moi43gtl2fzdffebm2tf/bin/python3.8" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "He_param_grad.param" "-e" "2" "--name" "jud_0" "--ref-value" "-0.121914" "--ref-error" "0.002")
set_tests_properties(short-He_param_grad-1-16-jud_0 PROPERTIES  DEPENDS "short-He_param_grad-1-16" LABELS "QMCPACK-checking-results;quality_unknown;unstable" WORKING_DIRECTORY "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/tests/molecules/He_param/short-He_param_grad-1-16" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;551;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;8;qmc_run_and_check_custom_scalar;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;0;")
add_test(short-He_param_grad-1-16-jud_1 "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/anaconda3-2021.05-hzr5w5v54sh3moi43gtl2fzdffebm2tf/bin/python3.8" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "He_param_grad.param" "-e" "2" "--name" "jud_1" "--ref-value" "0.069689" "--ref-error" "0.002")
set_tests_properties(short-He_param_grad-1-16-jud_1 PROPERTIES  DEPENDS "short-He_param_grad-1-16" LABELS "QMCPACK-checking-results;quality_unknown;unstable" WORKING_DIRECTORY "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/tests/molecules/He_param/short-He_param_grad-1-16" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;551;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;8;qmc_run_and_check_custom_scalar;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;0;")
add_test(short-He_param_grad-1-16-jud_2 "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/anaconda3-2021.05-hzr5w5v54sh3moi43gtl2fzdffebm2tf/bin/python3.8" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "He_param_grad.param" "-e" "2" "--name" "jud_2" "--ref-value" "0.051412" "--ref-error" "0.0005")
set_tests_properties(short-He_param_grad-1-16-jud_2 PROPERTIES  DEPENDS "short-He_param_grad-1-16" LABELS "QMCPACK-checking-results;quality_unknown;unstable" WORKING_DIRECTORY "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/tests/molecules/He_param/short-He_param_grad-1-16" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;551;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;8;qmc_run_and_check_custom_scalar;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;0;")
add_test(short-He_param_grad-1-16-jud_3 "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/anaconda3-2021.05-hzr5w5v54sh3moi43gtl2fzdffebm2tf/bin/python3.8" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "He_param_grad.param" "-e" "2" "--name" "jud_3" "--ref-value" "0.000812" "--ref-error" "0.000038")
set_tests_properties(short-He_param_grad-1-16-jud_3 PROPERTIES  DEPENDS "short-He_param_grad-1-16" LABELS "QMCPACK-checking-results;quality_unknown;unstable" WORKING_DIRECTORY "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/tests/molecules/He_param/short-He_param_grad-1-16" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;551;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;8;qmc_run_and_check_custom_scalar;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;0;")
add_test(short-He_param_grad_legacy_driver-1-16 "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpiexec" "-n" "1" "--bind-to" "none" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/bin/qmcpack_complex" "He_param_grad_legacy_driver.xml")
set_tests_properties(short-He_param_grad_legacy_driver-1-16 PROPERTIES  ENVIRONMENT "OMP_NUM_THREADS=16" FAIL_REGULAR_EXPRESSION "ERROR" LABELS "quality_unknown;QMCPACK;QMCPACK" PASS_REGULAR_EXPRESSION "QMCPACK execution completed successfully" PROCESSORS "16" PROCESSOR_AFFINITY "TRUE" WORKING_DIRECTORY "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/tests/molecules/He_param/short-He_param_grad_legacy_driver-1-16" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;145;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;234;run_qmc_app_no_copy;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;497;run_qmc_app;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;27;qmc_run_and_check_custom_scalar;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;0;")
add_test(short-He_param_grad_legacy_driver-1-16-jud_0 "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/anaconda3-2021.05-hzr5w5v54sh3moi43gtl2fzdffebm2tf/bin/python3.8" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "He_param_grad_legacy_driver.param" "-e" "2" "--name" "jud_0" "--ref-value" "-0.121914" "--ref-error" "0.002")
set_tests_properties(short-He_param_grad_legacy_driver-1-16-jud_0 PROPERTIES  DEPENDS "short-He_param_grad_legacy_driver-1-16" LABELS "QMCPACK-checking-results;quality_unknown;unstable" WORKING_DIRECTORY "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/tests/molecules/He_param/short-He_param_grad_legacy_driver-1-16" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;551;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;27;qmc_run_and_check_custom_scalar;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;0;")
add_test(short-He_param_grad_legacy_driver-1-16-jud_1 "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/anaconda3-2021.05-hzr5w5v54sh3moi43gtl2fzdffebm2tf/bin/python3.8" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "He_param_grad_legacy_driver.param" "-e" "2" "--name" "jud_1" "--ref-value" "0.069689" "--ref-error" "0.002")
set_tests_properties(short-He_param_grad_legacy_driver-1-16-jud_1 PROPERTIES  DEPENDS "short-He_param_grad_legacy_driver-1-16" LABELS "QMCPACK-checking-results;quality_unknown;unstable" WORKING_DIRECTORY "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/tests/molecules/He_param/short-He_param_grad_legacy_driver-1-16" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;551;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;27;qmc_run_and_check_custom_scalar;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;0;")
add_test(short-He_param_grad_legacy_driver-1-16-jud_2 "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/anaconda3-2021.05-hzr5w5v54sh3moi43gtl2fzdffebm2tf/bin/python3.8" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "He_param_grad_legacy_driver.param" "-e" "2" "--name" "jud_2" "--ref-value" "0.051412" "--ref-error" "0.0005")
set_tests_properties(short-He_param_grad_legacy_driver-1-16-jud_2 PROPERTIES  DEPENDS "short-He_param_grad_legacy_driver-1-16" LABELS "QMCPACK-checking-results;quality_unknown;unstable" WORKING_DIRECTORY "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/tests/molecules/He_param/short-He_param_grad_legacy_driver-1-16" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;551;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;27;qmc_run_and_check_custom_scalar;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;0;")
add_test(short-He_param_grad_legacy_driver-1-16-jud_3 "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/anaconda3-2021.05-hzr5w5v54sh3moi43gtl2fzdffebm2tf/bin/python3.8" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "He_param_grad_legacy_driver.param" "-e" "2" "--name" "jud_3" "--ref-value" "0.000812" "--ref-error" "0.000038")
set_tests_properties(short-He_param_grad_legacy_driver-1-16-jud_3 PROPERTIES  DEPENDS "short-He_param_grad_legacy_driver-1-16" LABELS "QMCPACK-checking-results;quality_unknown;unstable" WORKING_DIRECTORY "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/tests/molecules/He_param/short-He_param_grad_legacy_driver-1-16" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;551;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;27;qmc_run_and_check_custom_scalar;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;0;")
add_test(He_param_grad_load-1-16 "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpiexec" "-n" "1" "--bind-to" "none" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/bin/qmcpack_complex" "He_param_grad_load.xml")
set_tests_properties(He_param_grad_load-1-16 PROPERTIES  ENVIRONMENT "OMP_NUM_THREADS=16" FAIL_REGULAR_EXPRESSION "ERROR" LABELS "quality_unknown;QMCPACK;QMCPACK" PASS_REGULAR_EXPRESSION "QMCPACK execution completed successfully" PROCESSORS "16" PROCESSOR_AFFINITY "TRUE" WORKING_DIRECTORY "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/tests/molecules/He_param/He_param_grad_load-1-16" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;145;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;234;run_qmc_app_no_copy;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;497;run_qmc_app;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;66;qmc_run_and_check_custom_scalar;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;0;")
add_test(He_param_grad_load-1-16-jud_0 "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/anaconda3-2021.05-hzr5w5v54sh3moi43gtl2fzdffebm2tf/bin/python3.8" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "He_param_grad_load.param" "-e" "2" "--name" "jud_0" "--ref-value" "0.00000124" "--ref-error" "0.002")
set_tests_properties(He_param_grad_load-1-16-jud_0 PROPERTIES  DEPENDS "He_param_grad_load-1-16" LABELS "QMCPACK-checking-results;quality_unknown" WORKING_DIRECTORY "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/tests/molecules/He_param/He_param_grad_load-1-16" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;551;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;66;qmc_run_and_check_custom_scalar;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;0;")
add_test(He_param_grad_load-1-16-jud_1 "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/anaconda3-2021.05-hzr5w5v54sh3moi43gtl2fzdffebm2tf/bin/python3.8" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "He_param_grad_load.param" "-e" "2" "--name" "jud_1" "--ref-value" "-0.000273" "--ref-error" "0.0012")
set_tests_properties(He_param_grad_load-1-16-jud_1 PROPERTIES  DEPENDS "He_param_grad_load-1-16" LABELS "QMCPACK-checking-results;quality_unknown" WORKING_DIRECTORY "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/tests/molecules/He_param/He_param_grad_load-1-16" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;551;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;66;qmc_run_and_check_custom_scalar;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;0;")
add_test(He_param_grad_load-1-16-jud_2 "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/anaconda3-2021.05-hzr5w5v54sh3moi43gtl2fzdffebm2tf/bin/python3.8" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "He_param_grad_load.param" "-e" "2" "--name" "jud_2" "--ref-value" "-0.000181" "--ref-error" "0.00082")
set_tests_properties(He_param_grad_load-1-16-jud_2 PROPERTIES  DEPENDS "He_param_grad_load-1-16" LABELS "QMCPACK-checking-results;quality_unknown" WORKING_DIRECTORY "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/tests/molecules/He_param/He_param_grad_load-1-16" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;551;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;66;qmc_run_and_check_custom_scalar;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;0;")
add_test(He_param_grad_load-1-16-jud_3 "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/anaconda3-2021.05-hzr5w5v54sh3moi43gtl2fzdffebm2tf/bin/python3.8" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "He_param_grad_load.param" "-e" "2" "--name" "jud_3" "--ref-value" "0.0004463" "--ref-error" "0.00008")
set_tests_properties(He_param_grad_load-1-16-jud_3 PROPERTIES  DEPENDS "He_param_grad_load-1-16" LABELS "QMCPACK-checking-results;quality_unknown" WORKING_DIRECTORY "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/tests/molecules/He_param/He_param_grad_load-1-16" _BACKTRACE_TRIPLES "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/CMake/macros.cmake;551;add_test;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;66;qmc_run_and_check_custom_scalar;/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/tests/molecules/He_param/CMakeLists.txt;0;")
