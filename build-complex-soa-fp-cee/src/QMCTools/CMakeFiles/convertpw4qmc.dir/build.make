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
include src/QMCTools/CMakeFiles/convertpw4qmc.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/QMCTools/CMakeFiles/convertpw4qmc.dir/compiler_depend.make

# Include the progress variables for this target.
include src/QMCTools/CMakeFiles/convertpw4qmc.dir/progress.make

# Include the compile flags for this target's objects.
include src/QMCTools/CMakeFiles/convertpw4qmc.dir/flags.make

src/QMCTools/CMakeFiles/convertpw4qmc.dir/convertpw4qmc.cpp.o: src/QMCTools/CMakeFiles/convertpw4qmc.dir/flags.make
src/QMCTools/CMakeFiles/convertpw4qmc.dir/convertpw4qmc.cpp.o: ../src/QMCTools/convertpw4qmc.cpp
src/QMCTools/CMakeFiles/convertpw4qmc.dir/convertpw4qmc.cpp.o: src/QMCTools/CMakeFiles/convertpw4qmc.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/QMCTools/CMakeFiles/convertpw4qmc.dir/convertpw4qmc.cpp.o"
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools && /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/QMCTools/CMakeFiles/convertpw4qmc.dir/convertpw4qmc.cpp.o -MF CMakeFiles/convertpw4qmc.dir/convertpw4qmc.cpp.o.d -o CMakeFiles/convertpw4qmc.dir/convertpw4qmc.cpp.o -c /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/convertpw4qmc.cpp

src/QMCTools/CMakeFiles/convertpw4qmc.dir/convertpw4qmc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/convertpw4qmc.dir/convertpw4qmc.cpp.i"
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools && /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/convertpw4qmc.cpp > CMakeFiles/convertpw4qmc.dir/convertpw4qmc.cpp.i

src/QMCTools/CMakeFiles/convertpw4qmc.dir/convertpw4qmc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/convertpw4qmc.dir/convertpw4qmc.cpp.s"
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools && /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/convertpw4qmc.cpp -o CMakeFiles/convertpw4qmc.dir/convertpw4qmc.cpp.s

src/QMCTools/CMakeFiles/convertpw4qmc.dir/XmlRep.cpp.o: src/QMCTools/CMakeFiles/convertpw4qmc.dir/flags.make
src/QMCTools/CMakeFiles/convertpw4qmc.dir/XmlRep.cpp.o: ../src/QMCTools/XmlRep.cpp
src/QMCTools/CMakeFiles/convertpw4qmc.dir/XmlRep.cpp.o: src/QMCTools/CMakeFiles/convertpw4qmc.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/QMCTools/CMakeFiles/convertpw4qmc.dir/XmlRep.cpp.o"
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools && /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/QMCTools/CMakeFiles/convertpw4qmc.dir/XmlRep.cpp.o -MF CMakeFiles/convertpw4qmc.dir/XmlRep.cpp.o.d -o CMakeFiles/convertpw4qmc.dir/XmlRep.cpp.o -c /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/XmlRep.cpp

src/QMCTools/CMakeFiles/convertpw4qmc.dir/XmlRep.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/convertpw4qmc.dir/XmlRep.cpp.i"
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools && /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/XmlRep.cpp > CMakeFiles/convertpw4qmc.dir/XmlRep.cpp.i

src/QMCTools/CMakeFiles/convertpw4qmc.dir/XmlRep.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/convertpw4qmc.dir/XmlRep.cpp.s"
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools && /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/XmlRep.cpp -o CMakeFiles/convertpw4qmc.dir/XmlRep.cpp.s

src/QMCTools/CMakeFiles/convertpw4qmc.dir/WriteEshdf.cpp.o: src/QMCTools/CMakeFiles/convertpw4qmc.dir/flags.make
src/QMCTools/CMakeFiles/convertpw4qmc.dir/WriteEshdf.cpp.o: ../src/QMCTools/WriteEshdf.cpp
src/QMCTools/CMakeFiles/convertpw4qmc.dir/WriteEshdf.cpp.o: src/QMCTools/CMakeFiles/convertpw4qmc.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/QMCTools/CMakeFiles/convertpw4qmc.dir/WriteEshdf.cpp.o"
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools && /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/QMCTools/CMakeFiles/convertpw4qmc.dir/WriteEshdf.cpp.o -MF CMakeFiles/convertpw4qmc.dir/WriteEshdf.cpp.o.d -o CMakeFiles/convertpw4qmc.dir/WriteEshdf.cpp.o -c /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/WriteEshdf.cpp

src/QMCTools/CMakeFiles/convertpw4qmc.dir/WriteEshdf.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/convertpw4qmc.dir/WriteEshdf.cpp.i"
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools && /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/WriteEshdf.cpp > CMakeFiles/convertpw4qmc.dir/WriteEshdf.cpp.i

src/QMCTools/CMakeFiles/convertpw4qmc.dir/WriteEshdf.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/convertpw4qmc.dir/WriteEshdf.cpp.s"
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools && /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools/WriteEshdf.cpp -o CMakeFiles/convertpw4qmc.dir/WriteEshdf.cpp.s

# Object files for target convertpw4qmc
convertpw4qmc_OBJECTS = \
"CMakeFiles/convertpw4qmc.dir/convertpw4qmc.cpp.o" \
"CMakeFiles/convertpw4qmc.dir/XmlRep.cpp.o" \
"CMakeFiles/convertpw4qmc.dir/WriteEshdf.cpp.o"

# External object files for target convertpw4qmc
convertpw4qmc_EXTERNAL_OBJECTS =

bin/convertpw4qmc: src/QMCTools/CMakeFiles/convertpw4qmc.dir/convertpw4qmc.cpp.o
bin/convertpw4qmc: src/QMCTools/CMakeFiles/convertpw4qmc.dir/XmlRep.cpp.o
bin/convertpw4qmc: src/QMCTools/CMakeFiles/convertpw4qmc.dir/WriteEshdf.cpp.o
bin/convertpw4qmc: src/QMCTools/CMakeFiles/convertpw4qmc.dir/build.make
bin/convertpw4qmc: src/Utilities/libqmcutil.a
bin/convertpw4qmc: src/io/hdf/libqmcio_hdf.a
bin/convertpw4qmc: src/Message/libmessage.a
bin/convertpw4qmc: /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/hdf5-1.10.6-6wt7jwobpoa3iwvoyabbtoilhpsguglf/lib/libhdf5.so.103.2.0
bin/convertpw4qmc: src/io/OhmmsData/libqmcio_xml.a
bin/convertpw4qmc: src/Utilities/libcxx_helpers.a
bin/convertpw4qmc: src/Containers/libcontainers.a
bin/convertpw4qmc: src/Platforms/libplatform_runtime.a
bin/convertpw4qmc: src/Platforms/CPU/libplatform_cpu_runtime.a
bin/convertpw4qmc: src/Platforms/OMPTarget/libplatform_omptarget_runtime.a
bin/convertpw4qmc: src/Platforms/Host/libplatform_host_runtime.a
bin/convertpw4qmc: src/Utilities/libqmcrng.a
bin/convertpw4qmc: /usr/lib64/libxml2.so
bin/convertpw4qmc: /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/intel-oneapi-mkl-2021.4.0-uwfuzjcbz7p2ipcvv6wcrr3o5adh3ql2/mkl/2021.4.0/lib/intel64/libmkl_intel_lp64.so
bin/convertpw4qmc: /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/intel-oneapi-mkl-2021.4.0-uwfuzjcbz7p2ipcvv6wcrr3o5adh3ql2/mkl/2021.4.0/lib/intel64/libmkl_intel_thread.so
bin/convertpw4qmc: /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/intel-oneapi-mkl-2021.4.0-uwfuzjcbz7p2ipcvv6wcrr3o5adh3ql2/mkl/2021.4.0/lib/intel64/libmkl_core.so
bin/convertpw4qmc: /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/intel-oneapi-compilers-2021.1.2-lufoj3442adjwmyj2djuozq5aec3ofn2/compiler/2021.1.2/linux/compiler/lib/intel64_lin/libiomp5.so
bin/convertpw4qmc: src/QMCTools/CMakeFiles/convertpw4qmc.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable ../../bin/convertpw4qmc"
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/convertpw4qmc.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/QMCTools/CMakeFiles/convertpw4qmc.dir/build: bin/convertpw4qmc
.PHONY : src/QMCTools/CMakeFiles/convertpw4qmc.dir/build

src/QMCTools/CMakeFiles/convertpw4qmc.dir/clean:
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools && $(CMAKE_COMMAND) -P CMakeFiles/convertpw4qmc.dir/cmake_clean.cmake
.PHONY : src/QMCTools/CMakeFiles/convertpw4qmc.dir/clean

src/QMCTools/CMakeFiles/convertpw4qmc.dir/depend:
	cd /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/src/QMCTools /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee/src/QMCTools/CMakeFiles/convertpw4qmc.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/QMCTools/CMakeFiles/convertpw4qmc.dir/depend

