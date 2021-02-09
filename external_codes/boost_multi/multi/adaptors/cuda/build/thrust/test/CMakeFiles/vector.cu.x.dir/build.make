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

# Include any dependencies generated for this target.
include thrust/test/CMakeFiles/vector.cu.x.dir/depend.make

# Include the progress variables for this target.
include thrust/test/CMakeFiles/vector.cu.x.dir/progress.make

# Include the compile flags for this target's objects.
include thrust/test/CMakeFiles/vector.cu.x.dir/flags.make

thrust/test/CMakeFiles/vector.cu.x.dir/vector.cu.o: thrust/test/CMakeFiles/vector.cu.x.dir/flags.make
thrust/test/CMakeFiles/vector.cu.x.dir/vector.cu.o: ../thrust/test/vector.cu
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/correaa/prj/inq/external_libs/multi/adaptors/cuda/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CUDA object thrust/test/CMakeFiles/vector.cu.x.dir/vector.cu.o"
	cd /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/build/thrust/test && /usr/local/cuda/bin/nvcc  $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -x cu -c /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/thrust/test/vector.cu -o CMakeFiles/vector.cu.x.dir/vector.cu.o

thrust/test/CMakeFiles/vector.cu.x.dir/vector.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CUDA source to CMakeFiles/vector.cu.x.dir/vector.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

thrust/test/CMakeFiles/vector.cu.x.dir/vector.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CUDA source to assembly CMakeFiles/vector.cu.x.dir/vector.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

# Object files for target vector.cu.x
vector_cu_x_OBJECTS = \
"CMakeFiles/vector.cu.x.dir/vector.cu.o"

# External object files for target vector.cu.x
vector_cu_x_EXTERNAL_OBJECTS =

thrust/test/vector.cu.x: thrust/test/CMakeFiles/vector.cu.x.dir/vector.cu.o
thrust/test/vector.cu.x: thrust/test/CMakeFiles/vector.cu.x.dir/build.make
thrust/test/vector.cu.x: /usr/lib/x86_64-linux-gnu/libboost_unit_test_framework.so.1.71.0
thrust/test/vector.cu.x: thrust/test/CMakeFiles/vector.cu.x.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/correaa/prj/inq/external_libs/multi/adaptors/cuda/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CUDA executable vector.cu.x"
	cd /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/build/thrust/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/vector.cu.x.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
thrust/test/CMakeFiles/vector.cu.x.dir/build: thrust/test/vector.cu.x

.PHONY : thrust/test/CMakeFiles/vector.cu.x.dir/build

thrust/test/CMakeFiles/vector.cu.x.dir/clean:
	cd /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/build/thrust/test && $(CMAKE_COMMAND) -P CMakeFiles/vector.cu.x.dir/cmake_clean.cmake
.PHONY : thrust/test/CMakeFiles/vector.cu.x.dir/clean

thrust/test/CMakeFiles/vector.cu.x.dir/depend:
	cd /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/correaa/prj/inq/external_libs/multi/adaptors/cuda /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/thrust/test /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/build /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/build/thrust/test /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/build/thrust/test/CMakeFiles/vector.cu.x.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : thrust/test/CMakeFiles/vector.cu.x.dir/depend

