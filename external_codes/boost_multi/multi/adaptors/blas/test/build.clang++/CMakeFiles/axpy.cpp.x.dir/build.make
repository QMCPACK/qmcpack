# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Produce verbose output by default.
VERBOSE = 1

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/build.clang++

# Include any dependencies generated for this target.
include CMakeFiles/axpy.cpp.x.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/axpy.cpp.x.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/axpy.cpp.x.dir/flags.make

CMakeFiles/axpy.cpp.x.dir/axpy.cpp.o: CMakeFiles/axpy.cpp.x.dir/flags.make
CMakeFiles/axpy.cpp.x.dir/axpy.cpp.o: ../axpy.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/build.clang++/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/axpy.cpp.x.dir/axpy.cpp.o"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/axpy.cpp.x.dir/axpy.cpp.o -c /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/axpy.cpp

CMakeFiles/axpy.cpp.x.dir/axpy.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/axpy.cpp.x.dir/axpy.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/axpy.cpp > CMakeFiles/axpy.cpp.x.dir/axpy.cpp.i

CMakeFiles/axpy.cpp.x.dir/axpy.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/axpy.cpp.x.dir/axpy.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/axpy.cpp -o CMakeFiles/axpy.cpp.x.dir/axpy.cpp.s

# Object files for target axpy.cpp.x
axpy_cpp_x_OBJECTS = \
"CMakeFiles/axpy.cpp.x.dir/axpy.cpp.o"

# External object files for target axpy.cpp.x
axpy_cpp_x_EXTERNAL_OBJECTS =

axpy.cpp.x: CMakeFiles/axpy.cpp.x.dir/axpy.cpp.o
axpy.cpp.x: CMakeFiles/axpy.cpp.x.dir/build.make
axpy.cpp.x: /usr/lib/x86_64-linux-gnu/libopenblas.so
axpy.cpp.x: /usr/lib/x86_64-linux-gnu/libboost_unit_test_framework.so.1.71.0
axpy.cpp.x: CMakeFiles/axpy.cpp.x.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/build.clang++/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable axpy.cpp.x"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/axpy.cpp.x.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/axpy.cpp.x.dir/build: axpy.cpp.x

.PHONY : CMakeFiles/axpy.cpp.x.dir/build

CMakeFiles/axpy.cpp.x.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/axpy.cpp.x.dir/cmake_clean.cmake
.PHONY : CMakeFiles/axpy.cpp.x.dir/clean

CMakeFiles/axpy.cpp.x.dir/depend:
	cd /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/build.clang++ && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/build.clang++ /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/build.clang++ /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/build.clang++/CMakeFiles/axpy.cpp.x.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/axpy.cpp.x.dir/depend

