# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/cmake-3.23.1-v3pkd6okl3lzh723twjgmlfqs54zqrdv/bin/cmake

# The command to remove a file.
RM = /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/cmake-3.23.1-v3pkd6okl3lzh723twjgmlfqs54zqrdv/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee

# Include any dependencies generated for this target.
include src/QMCTools/tests/CMakeFiles/test_tools.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/QMCTools/tests/CMakeFiles/test_tools.dir/compiler_depend.make

# Include the progress variables for this target.
include src/QMCTools/tests/CMakeFiles/test_tools.dir/progress.make

# Include the compile flags for this target's objects.
include src/QMCTools/tests/CMakeFiles/test_tools.dir/flags.make

src/QMCTools/tests/CMakeFiles/test_tools.dir/test_qmcfstool.cpp.o: src/QMCTools/tests/CMakeFiles/test_tools.dir/flags.make
src/QMCTools/tests/CMakeFiles/test_tools.dir/test_qmcfstool.cpp.o: ../src/QMCTools/tests/test_qmcfstool.cpp
src/QMCTools/tests/CMakeFiles/test_tools.dir/test_qmcfstool.cpp.o: src/QMCTools/tests/CMakeFiles/test_tools.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/QMCTools/tests/CMakeFiles/test_tools.dir/test_qmcfstool.cpp.o"
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools/tests && /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/QMCTools/tests/CMakeFiles/test_tools.dir/test_qmcfstool.cpp.o -MF CMakeFiles/test_tools.dir/test_qmcfstool.cpp.o.d -o CMakeFiles/test_tools.dir/test_qmcfstool.cpp.o -c /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/tests/test_qmcfstool.cpp

src/QMCTools/tests/CMakeFiles/test_tools.dir/test_qmcfstool.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_tools.dir/test_qmcfstool.cpp.i"
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools/tests && /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/tests/test_qmcfstool.cpp > CMakeFiles/test_tools.dir/test_qmcfstool.cpp.i

src/QMCTools/tests/CMakeFiles/test_tools.dir/test_qmcfstool.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_tools.dir/test_qmcfstool.cpp.s"
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools/tests && /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/tests/test_qmcfstool.cpp -o CMakeFiles/test_tools.dir/test_qmcfstool.cpp.s

# Object files for target test_tools
test_tools_OBJECTS = \
"CMakeFiles/test_tools.dir/test_qmcfstool.cpp.o"

# External object files for target test_tools
test_tools_EXTERNAL_OBJECTS =

src/QMCTools/tests/test_tools: src/QMCTools/tests/CMakeFiles/test_tools.dir/test_qmcfstool.cpp.o
src/QMCTools/tests/test_tools: src/QMCTools/tests/CMakeFiles/test_tools.dir/build.make
src/QMCTools/tests/test_tools: src/Message/libcatch_main.a
src/QMCTools/tests/test_tools: src/QMCTools/libfstool.a
src/QMCTools/tests/test_tools: src/QMCApp/libqmc.a
src/QMCTools/tests/test_tools: src/QMCDrivers/libqmcdriver.a
src/QMCTools/tests/test_tools: src/Estimators/libqmcestimators.a
src/QMCTools/tests/test_tools: src/QMCHamiltonians/libqmcham.a
src/QMCTools/tests/test_tools: src/formic/utils/libformic_utils.a
src/QMCTools/tests/test_tools: src/QMCWaveFunctions/libqmcwfs.a
src/QMCTools/tests/test_tools: src/einspline/libeinspline.a
src/QMCTools/tests/test_tools: src/Particle/libqmcparticle.a
src/QMCTools/tests/test_tools: src/Numerics/libqmcnumerics.a
src/QMCTools/tests/test_tools: src/Utilities/libqmcutil.a
src/QMCTools/tests/test_tools: src/io/hdf/libqmcio_hdf.a
src/QMCTools/tests/test_tools: src/Message/libmessage.a
src/QMCTools/tests/test_tools: /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/hdf5-1.10.6-6wt7jwobpoa3iwvoyabbtoilhpsguglf/lib/libhdf5.so.103.2.0
src/QMCTools/tests/test_tools: src/io/OhmmsData/libqmcio_xml.a
src/QMCTools/tests/test_tools: src/Containers/libcontainers.a
src/QMCTools/tests/test_tools: src/Utilities/libcxx_helpers.a
src/QMCTools/tests/test_tools: src/Utilities/libqmcrng.a
src/QMCTools/tests/test_tools: /usr/lib64/libxml2.so
src/QMCTools/tests/test_tools: src/Platforms/libplatform_runtime.a
src/QMCTools/tests/test_tools: src/Platforms/CPU/libplatform_cpu_LA.a
src/QMCTools/tests/test_tools: src/Platforms/CPU/libplatform_cpu_runtime.a
src/QMCTools/tests/test_tools: /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/intel-oneapi-mkl-2021.4.0-uwfuzjcbz7p2ipcvv6wcrr3o5adh3ql2/mkl/2021.4.0/lib/intel64/libmkl_intel_lp64.so
src/QMCTools/tests/test_tools: /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/intel-oneapi-mkl-2021.4.0-uwfuzjcbz7p2ipcvv6wcrr3o5adh3ql2/mkl/2021.4.0/lib/intel64/libmkl_intel_thread.so
src/QMCTools/tests/test_tools: /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/intel-oneapi-mkl-2021.4.0-uwfuzjcbz7p2ipcvv6wcrr3o5adh3ql2/mkl/2021.4.0/lib/intel64/libmkl_core.so
src/QMCTools/tests/test_tools: /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/intel-oneapi-compilers-2021.1.2-lufoj3442adjwmyj2djuozq5aec3ofn2/compiler/2021.1.2/linux/compiler/lib/intel64_lin/libiomp5.so
src/QMCTools/tests/test_tools: src/Platforms/OMPTarget/libplatform_omptarget_LA.a
src/QMCTools/tests/test_tools: src/Platforms/OMPTarget/libplatform_omptarget_runtime.a
src/QMCTools/tests/test_tools: src/Platforms/Host/libplatform_host_runtime.a
src/QMCTools/tests/test_tools: src/QMCTools/tests/CMakeFiles/test_tools.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test_tools"
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_tools.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/QMCTools/tests/CMakeFiles/test_tools.dir/build: src/QMCTools/tests/test_tools
.PHONY : src/QMCTools/tests/CMakeFiles/test_tools.dir/build

src/QMCTools/tests/CMakeFiles/test_tools.dir/clean:
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools/tests && $(CMAKE_COMMAND) -P CMakeFiles/test_tools.dir/cmake_clean.cmake
.PHONY : src/QMCTools/tests/CMakeFiles/test_tools.dir/clean

src/QMCTools/tests/CMakeFiles/test_tools.dir/depend:
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/tests /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools/tests /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools/tests/CMakeFiles/test_tools.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/QMCTools/tests/CMakeFiles/test_tools.dir/depend

