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
include src/Utilities/CMakeFiles/cxx_helpers.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/Utilities/CMakeFiles/cxx_helpers.dir/compiler_depend.make

# Include the progress variables for this target.
include src/Utilities/CMakeFiles/cxx_helpers.dir/progress.make

# Include the compile flags for this target's objects.
include src/Utilities/CMakeFiles/cxx_helpers.dir/flags.make

src/Utilities/CMakeFiles/cxx_helpers.dir/ModernStringUtils.cpp.o: src/Utilities/CMakeFiles/cxx_helpers.dir/flags.make
src/Utilities/CMakeFiles/cxx_helpers.dir/ModernStringUtils.cpp.o: ../src/Utilities/ModernStringUtils.cpp
src/Utilities/CMakeFiles/cxx_helpers.dir/ModernStringUtils.cpp.o: src/Utilities/CMakeFiles/cxx_helpers.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/Utilities/CMakeFiles/cxx_helpers.dir/ModernStringUtils.cpp.o"
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/Utilities && /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/Utilities/CMakeFiles/cxx_helpers.dir/ModernStringUtils.cpp.o -MF CMakeFiles/cxx_helpers.dir/ModernStringUtils.cpp.o.d -o CMakeFiles/cxx_helpers.dir/ModernStringUtils.cpp.o -c /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/Utilities/ModernStringUtils.cpp

src/Utilities/CMakeFiles/cxx_helpers.dir/ModernStringUtils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cxx_helpers.dir/ModernStringUtils.cpp.i"
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/Utilities && /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/Utilities/ModernStringUtils.cpp > CMakeFiles/cxx_helpers.dir/ModernStringUtils.cpp.i

src/Utilities/CMakeFiles/cxx_helpers.dir/ModernStringUtils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cxx_helpers.dir/ModernStringUtils.cpp.s"
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/Utilities && /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/Utilities/ModernStringUtils.cpp -o CMakeFiles/cxx_helpers.dir/ModernStringUtils.cpp.s

# Object files for target cxx_helpers
cxx_helpers_OBJECTS = \
"CMakeFiles/cxx_helpers.dir/ModernStringUtils.cpp.o"

# External object files for target cxx_helpers
cxx_helpers_EXTERNAL_OBJECTS =

src/Utilities/libcxx_helpers.a: src/Utilities/CMakeFiles/cxx_helpers.dir/ModernStringUtils.cpp.o
src/Utilities/libcxx_helpers.a: src/Utilities/CMakeFiles/cxx_helpers.dir/build.make
src/Utilities/libcxx_helpers.a: src/Utilities/CMakeFiles/cxx_helpers.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libcxx_helpers.a"
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/Utilities && $(CMAKE_COMMAND) -P CMakeFiles/cxx_helpers.dir/cmake_clean_target.cmake
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/Utilities && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cxx_helpers.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/Utilities/CMakeFiles/cxx_helpers.dir/build: src/Utilities/libcxx_helpers.a
.PHONY : src/Utilities/CMakeFiles/cxx_helpers.dir/build

src/Utilities/CMakeFiles/cxx_helpers.dir/clean:
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/Utilities && $(CMAKE_COMMAND) -P CMakeFiles/cxx_helpers.dir/cmake_clean.cmake
.PHONY : src/Utilities/CMakeFiles/cxx_helpers.dir/clean

src/Utilities/CMakeFiles/cxx_helpers.dir/depend:
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/Utilities /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/Utilities /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/Utilities/CMakeFiles/cxx_helpers.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/Utilities/CMakeFiles/cxx_helpers.dir/depend

