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
CMAKE_SOURCE_DIR = /home/correaa/prj/inq/external_libs/multi/adaptors/cuda

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/build

# Utility rule file for Nightly.

# Include the progress variables for this target.
include thrust/test/CMakeFiles/Nightly.dir/progress.make

thrust/test/CMakeFiles/Nightly:
	cd /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/build/thrust/test && /usr/bin/ctest -D Nightly

Nightly: thrust/test/CMakeFiles/Nightly
Nightly: thrust/test/CMakeFiles/Nightly.dir/build.make

.PHONY : Nightly

# Rule to build all files generated by this target.
thrust/test/CMakeFiles/Nightly.dir/build: Nightly

.PHONY : thrust/test/CMakeFiles/Nightly.dir/build

thrust/test/CMakeFiles/Nightly.dir/clean:
	cd /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/build/thrust/test && $(CMAKE_COMMAND) -P CMakeFiles/Nightly.dir/cmake_clean.cmake
.PHONY : thrust/test/CMakeFiles/Nightly.dir/clean

thrust/test/CMakeFiles/Nightly.dir/depend:
	cd /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/correaa/prj/inq/external_libs/multi/adaptors/cuda /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/thrust/test /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/build /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/build/thrust/test /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/build/thrust/test/CMakeFiles/Nightly.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : thrust/test/CMakeFiles/Nightly.dir/depend

