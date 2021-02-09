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
CMAKE_SOURCE_DIR = /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/thrust/test

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/thrust/test/build

# Include any dependencies generated for this target.
include CMakeFiles/array.cu.x.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/array.cu.x.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/array.cu.x.dir/flags.make

CMakeFiles/array.cu.x.dir/array.cu.o: CMakeFiles/array.cu.x.dir/flags.make
CMakeFiles/array.cu.x.dir/array.cu.o: ../array.cu
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/correaa/prj/inq/external_libs/multi/adaptors/cuda/thrust/test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CUDA object CMakeFiles/array.cu.x.dir/array.cu.o"
	/usr/local/cuda-11.1/bin/nvcc  $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -x cu -c /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/thrust/test/array.cu -o CMakeFiles/array.cu.x.dir/array.cu.o

CMakeFiles/array.cu.x.dir/array.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CUDA source to CMakeFiles/array.cu.x.dir/array.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/array.cu.x.dir/array.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CUDA source to assembly CMakeFiles/array.cu.x.dir/array.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

# Object files for target array.cu.x
array_cu_x_OBJECTS = \
"CMakeFiles/array.cu.x.dir/array.cu.o"

# External object files for target array.cu.x
array_cu_x_EXTERNAL_OBJECTS =

array.cu.x: CMakeFiles/array.cu.x.dir/array.cu.o
array.cu.x: CMakeFiles/array.cu.x.dir/build.make
array.cu.x: /usr/lib/x86_64-linux-gnu/libboost_unit_test_framework.so.1.71.0
array.cu.x: CMakeFiles/array.cu.x.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/correaa/prj/inq/external_libs/multi/adaptors/cuda/thrust/test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CUDA executable array.cu.x"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/array.cu.x.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/array.cu.x.dir/build: array.cu.x

.PHONY : CMakeFiles/array.cu.x.dir/build

CMakeFiles/array.cu.x.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/array.cu.x.dir/cmake_clean.cmake
.PHONY : CMakeFiles/array.cu.x.dir/clean

CMakeFiles/array.cu.x.dir/depend:
	cd /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/thrust/test/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/thrust/test /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/thrust/test /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/thrust/test/build /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/thrust/test/build /home/correaa/prj/inq/external_libs/multi/adaptors/cuda/thrust/test/build/CMakeFiles/array.cu.x.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/array.cu.x.dir/depend

