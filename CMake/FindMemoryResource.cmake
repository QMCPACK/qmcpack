# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:

FindMemoryResource
##############

This module supports the C++17 standard library's memory resource specification. Use the
:imp-target:`std::memory_resource` imported target to use this.

Because at the moment it is unlikely that a libc++ will support this.

Options
*******

The ``COMPONENTS`` argument to this module supports the following values:

.. find-component:: Experimental
    :name: fs.Experimental

    Allows the module to find the "experimental" memory_resource TS version of the
    memory_resource library. This is the library that should be used with the
    ``std::experimental::memory_resource`` namespace.

.. find-component:: Final
    :name: fs.Final

    Finds the final C++17 standard version of the memory_resource library.

If no components are provided, behaves as if the
:find-component:`fs.Final` component was specified.

If both :find-component:`fs.Experimental` and :find-component:`fs.Final` are
provided, first looks for ``Final``, and falls back to ``Experimental`` in case
of failure. If ``Final`` is found, :imp-target:`std::memory_resource` and all
:ref:`variables <fs.variables>` will refer to the ``Final`` version.


Imported Targets
****************

.. imp-target:: std::memory_resource

    The ``std::memory_resource`` imported target is defined when any requested
    version of the C++ memory_resource library has been found, whether it is
    *Experimental* or *Final*.

    If no version of the memory_resource library is available, this target will not
    be defined.

    .. note::
        This target has ``cxx_std_17`` as an ``INTERFACE``
        :ref:`compile language standard feature <req-lang-standards>`. Linking
        to this target will automatically enable C++17 if no later standard
        version is already required on the linking target.


.. _mr.variables:

Variables
*********

.. variable:: CXX_MEMORYRESOURCE_IS_EXPERIMENTAL

    Set to ``TRUE`` when the :find-component:`mr.Experimental` version of C++
    memory_resource library was found, otherwise ``FALSE``.

.. variable:: CXX_MEMORYRESOURCE_HAVE_MR

    Set to ``TRUE`` when a memory_resource header was found.

.. variable:: CXX_MEMORYRESOURCE_HEADER

    Set to either ``memory_resource`` or ``experimental/memory_resource`` depending on
    whether :find-component:`mr.Final` or :find-component:`mr.Experimental` was
    found.

.. variable:: CXX_MEMORYRESOURCE_NAMESPACE

    Set to either ``std::memory_resource`` or ``std::experimental::memory_resource``
    depending on whether :find-component:`mr.Final` or
    :find-component:`mr.Experimental` was found.


Examples
********

Using `find_package(MemoryResource)` with no component arguments:

.. code-block:: cmake

    find_package(MemoryResource REQUIRED)

    add_executable(my-program main.cpp)
    target_link_libraries(my-program PRIVATE std::memory_resource)


#]=======================================================================]


if(TARGET std::memory_resource)
    # This module has already been processed. Don't do it again.
    return()
endif()

cmake_minimum_required(VERSION 3.10)

include(CMakePushCheckState)
include(CheckIncludeFileCXX)

# If we're not cross-compiling, try to run test executables.
# Otherwise, assume that compile + link is a sufficient check.
if(CMAKE_CROSSCOMPILING)
    include(CheckCXXSourceCompiles)
    macro(_cmcm_check_cxx_source code var)
        check_cxx_source_compiles("${code}" ${var})
    endmacro()
else()
    include(CheckCXXSourceRuns)
    macro(_cmcm_check_cxx_source code var)
        check_cxx_source_runs("${code}" ${var})
    endmacro()
endif()

cmake_push_check_state()

#set(CMAKE_REQUIRED_QUIET ${memory_resource_FIND_QUIETLY})

# All of our tests required C++17 or later
set(CMAKE_CXX_STANDARD 17)

# Normalize and check the component list we were given
set(want_components ${memory_resource_FIND_COMPONENTS})
message ("memory_resource_FIND_COMPONENTS ${memory_resource_FIND_COMPONENTS}")
set(want_components Final Experimental)

message("want_component at ${want_components}")

# Warn on any unrecognized components
set(extra_components ${want_components})
list(REMOVE_ITEM extra_components Final Experimental)
foreach(component IN LISTS extra_components)
    message(WARNING "Extraneous find_package component for memory_resource: ${component}")
endforeach()

# Detect which of Experimental and Final we should look for
set(find_experimental TRUE)
set(find_final TRUE)
if(NOT "Final" IN_LIST want_components)
  set(find_final FALSE)
  message("Final not in list")
endif()
if(NOT "Experimental" IN_LIST want_components)
    set(find_experimental FALSE)
endif()

set(_CXX_MEMORYRESOURCE_HAVE_HEADER "" CACHE INTERNAL "")

if(find_final)
  message("searching for memory_resource")
  check_include_file_cxx(memory_resource CXX_MEMORYRESOURCE_HAVE_HEADER)
  message("_CXX_... ${CXX_MEMORYRESOURCE_HAVE_HEADER}" )
# workaround single check include file cxx seems throughly messed up on osx
  if(CXX_MEMORYRESOURCE_HAVE_HEADER)
        # We found the non-experimental header. Don't bother looking for the
        # experimental one.
        message("Found <memory_resource> header")
        set(find_experimental FALSE)
    endif()
else()
    set(CXX_MEMORYRESOURCE_HAVE_HEADER FALSE)
endif()

if(find_experimental)
  message("searching for experimental/memory_resource")
  check_include_file_cxx(experimental/memory_resource WHATEVER)
  message("_CXX_... ${_CXX_MEMORYRESOURCE_HAVE_EXPERIMENTAL_HEADER} ${WHATEVER}")
  if(WHATEVER)
    set(_CXX_MEMORYRESOURCE_HAVE_EXPERIMENTAL_HEADER 1)
  endif()
else()
    set(_CXX_MEMORYRESOURCE_HAVE_EXPERIMENTAL_HEADER FALSE)
endif()

if(CXX_MEMORYRESOURCE_HAVE_HEADER)
    set(_have_mr TRUE)
    set(_mr_header memory_resource)
    set(_mr_namespace std::pmr)
    set(_is_experimental FALSE)
elseif(_CXX_MEMORYRESOURCE_HAVE_EXPERIMENTAL_HEADER)
    set(_have_mr TRUE)
    set(_mr_header experimental/memory_resource)
    set(_mr_namespace std::experimental::pmr)
    set(_is_experimental TRUE)
else()
    set(_have_mr FALSE)
endif()

set(CXX_MEMORYRESOURCE_HAVE_MR ${_have_mr} CACHE BOOL "TRUE if we have the C++ memory_resource headers")
message("have C++ memory_resource header: ${CXX_MEMORYRESOURCE_HAVE_MR}")
set(CXX_MEMORYRESOURCE_HEADER ${_mr_header} CACHE STRING "The header that should be included to obtain the memory_resource APIs")
set(CXX_MEMORYRESOURCE_NAMESPACE ${_mr_namespace} CACHE STRING "The C++ namespace that contains the memory_resource APIs")
set(CXX_MEMORYRESOURCE_IS_EXPERIMENTAL ${_is_experimental} CACHE BOOL "TRUE if the C++ memory_resource library is the experimental version")

message ("CXX_MEMORYRESOURCE_HEADER: ${CXX_MEMORYRESOURCE_HEADER}")
message ("CXX_MEMORYRESOURCE_NAMESPACE: ${CXX_MEMORYRESOURCE_NAMESPACE}")

set(_found FALSE)

if(CXX_MEMORYRESOURCE_HAVE_MR)
    # We have some memory_resource library available. Do link checks
    string(CONFIGURE [[
        #include <cstdlib>
        #include <@CXX_MEMORYRESOURCE_HEADER@>

        int main() {
            auto mem_resource = @CXX_MEMORYRESOURCE_NAMESPACE@::get_default_resource();
            auto ptr = mem_resource->allocate(256,64);
            mem_resource->deallocate(ptr,256,64);
            return EXIT_SUCCESS;
        }
    ]] code @ONLY)

    # Check a simple memory_resource program without any linker flags
    _cmcm_check_cxx_source("${code}" CXX_MEMORYRESOURCE_NO_LINK_NEEDED)

    set(can_link ${CXX_MEMORYRESOURCE_NO_LINK_NEEDED})
    message("Can build without memory_resource lib flag: ${CXX_MEMORYRESOURCE_NO_LINK_NEEDED}")
endif()
    
    if(NOT CXX_MEMORYRESOURCE_NO_LINK_NEEDED)
        set(prev_libraries ${CMAKE_REQUIRED_LIBRARIES})
        # Add the libstdc++ flag
        set(CMAKE_REQUIRED_LIBRARIES ${prev_libraries} -lstdc++)
        message("Trying libstdc++")
        _cmcm_check_cxx_source("${code}" CXX_MEMORYRESOURCE_STDCPPFS_NEEDED)
        set(can_link ${CXX_MEMORYRESOURCE_STDCPPFS_NEEDED})
        if(NOT CXX_MEMORYRESOURCE_STDCPPFS_NEEDED)
            message("Trying libc++experimental flag")
            set(CMAKE_REQUIRED_LIBRARIES ${prev_libraries} -lc++experimental)
            _cmcm_check_cxx_source("${code}" CXX_MEMORYRESOURCE_CPPFS_NEEDED)
            set(can_link ${CXX_MEMORYRESOURCE_CPPFS_NEEDED})
        endif()
    endif()

#     if(can_link)
#         add_library(std::memory_resource INTERFACE IMPORTED)
#         set_property(TARGET std::memory_resource APPEND PROPERTY INTERFACE_COMPILE_FEATURES cxx_std_17)
#         set(_found TRUE)

#     endif()
# endif()

# cmake_pop_check_state()

# set(MemoryResource_FOUND ${_found} CACHE BOOL "TRUE if we can run a program using std::memory_resource" FORCE)

# if(MemoryResource_FIND_REQUIRED AND NOT MemoryResource_FOUND)
#     message(FATAL_ERROR "Cannot find std::memory_resource")
# endif()
