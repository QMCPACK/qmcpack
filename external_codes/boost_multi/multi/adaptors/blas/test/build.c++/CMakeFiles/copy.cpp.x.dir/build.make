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
CMAKE_BINARY_DIR = /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/build.c++

# Include any dependencies generated for this target.
include CMakeFiles/copy.cpp.x.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/copy.cpp.x.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/copy.cpp.x.dir/flags.make

CMakeFiles/copy.cpp.x.dir/copy.cpp.o: CMakeFiles/copy.cpp.x.dir/flags.make
CMakeFiles/copy.cpp.x.dir/copy.cpp.o: ../copy.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/build.c++/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/copy.cpp.x.dir/copy.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/copy.cpp.x.dir/copy.cpp.o -c /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/copy.cpp

CMakeFiles/copy.cpp.x.dir/copy.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/copy.cpp.x.dir/copy.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/copy.cpp > CMakeFiles/copy.cpp.x.dir/copy.cpp.i

CMakeFiles/copy.cpp.x.dir/copy.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/copy.cpp.x.dir/copy.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/copy.cpp -o CMakeFiles/copy.cpp.x.dir/copy.cpp.s

# Object files for target copy.cpp.x
copy_cpp_x_OBJECTS = \
"CMakeFiles/copy.cpp.x.dir/copy.cpp.o"

# External object files for target copy.cpp.x
copy_cpp_x_EXTERNAL_OBJECTS =

copy.cpp.x: CMakeFiles/copy.cpp.x.dir/copy.cpp.o
copy.cpp.x: CMakeFiles/copy.cpp.x.dir/build.make
copy.cpp.x: /usr/lib/x86_64-linux-gnu/libopenblas.so
copy.cpp.x: /usr/lib/x86_64-linux-gnu/libboost_unit_test_framework.so.1.71.0
copy.cpp.x: CMakeFiles/copy.cpp.x.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/build.c++/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable copy.cpp.x"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/copy.cpp.x.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/copy.cpp.x.dir/build: copy.cpp.x

.PHONY : CMakeFiles/copy.cpp.x.dir/build

CMakeFiles/copy.cpp.x.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/copy.cpp.x.dir/cmake_clean.cmake
.PHONY : CMakeFiles/copy.cpp.x.dir/clean

CMakeFiles/copy.cpp.x.dir/depend:
	cd /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/build.c++ && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/build.c++ /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/build.c++ /home/correaa/prj/inq/external_libs/multi/adaptors/blas/test/build.c++/CMakeFiles/copy.cpp.x.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/copy.cpp.x.dir/depend

