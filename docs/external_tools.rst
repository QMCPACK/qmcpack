.. _external-tools:

External Tools
==============

This chapter provides some information on using QMCPACK with external tools.

.. _Sanitizer-Libraries:

Sanitizer Libraries
-------------------

Using CMake, set one of these flags for using the clang sanitizer libraries with or without lldb.

::

   -DENABLE_SANITIZER  link with the GNU or Clang sanitizer library for asan, ubsan, tsan or msan (default=none)
   
In general: 

- address sanitizer (asan):  catches most pointer-based errors and memory leaks (via lsan) by default. 
- undefined behavior sanitizer (ubsan): low-overhead, catches undefined behavior accessing misaligned memory or signed or float to integer overflows.
- undefined behavior sanitizer (tsan): catches potential race conditions in threaded code.
- memory sanitizer (msan): catches using uninitialized memory errors, but is difficult to use without a full set of msan-instrumented libraries.

These set the basic flags required to build with either of these sanitizer libraries which are mutually exclusive. Depending on your system and linker, these may be incompatible with the "Release" build, so set ``-DCMAKE_BUILD_TYPE=Debug`` or ``-DCMAKE_BUILD_TYPE=RelWithDebInfo``. They are tested on GitHub Actions CI using deterministic tests ``ctest -L deterministic`` (currently ubsan). See the following links for additional information on use, run time, and build options of the sanitizers: https://clang.llvm.org/docs/AddressSanitizer.html & https://clang.llvm.org/docs/MemorySanitizer.html.


Intel VTune
-----------

Intel's VTune profiler has an API that allows program control over collection (pause/resume) and can add information to the profile data (e.g., delineating tasks).

VTune API
~~~~~~~~~

If the variable ``USE_VTUNE_API`` is set, QMCPACK will check that the
include file (``ittnotify.h``) and the library (``libittnotify.a``) can be found.
To provide CMake with the VTune search paths, add ``VTUNE_ROOT`` which contains ``include`` and ``lib64`` sub-directories.

An example of options to be passed to CMake:

::

  -DUSE_VTUNE_API=ON \
  -DVTUNE_ROOT=/opt/intel/vtune_amplifier_xe

Timers as Tasks
~~~~~~~~~~~~~~~

To aid in connecting the timers in the code to the profile data, the start/stop of
timers will be recorded as a task if ``USE_VTUNE_TASKS`` is set.

In addition to compiling with ``USE_VTUNE_TASKS``, an option needs to be set at run time to collect the task API data.
In the graphical user interface (GUI), select the checkbox labeled "Analyze user tasks" when setting up the analysis type.
For the command line, set the ``enable-user-tasks`` knob to ``true``. For example,

::

  amplxe-cl -collect hotspots -knob enable-user-tasks=true ...

Collection with the timers set at "fine" can generate too much task data in the profile.
Collection with the timers at "medium" collects a more reasonable amount of task data.

NVIDIA Tools Extensions
-----------------------

NVIDIA's Tools Extensions (NVTX) API enables programmers to annotate their source code when used with the NVIDIA profilers.

NVTX API
~~~~~~~~

If the variable ``USE_NVTX_API`` is set, QMCPACK will add the library (``libnvToolsExt.so``) to the QMCPACK target. To add NVTX annotations
to a function, it is necessary to include the ``nvToolsExt.h`` header file and then make the appropriate calls into the NVTX API. For more information
about the NVTX API, see https://docs.nvidia.com/cuda/profiler-users-guide/index.html#nvtx. Any additional calls to the NVTX API should be guarded by
the ``USE_NVTX_API`` compiler define.

Scitools Understand
-------------------

Scitools Understand (https://scitools.com/) is a tool for static
code analysis. The easiest configuration route is to use the JSON output
from CMake, which the Understand project importer can read directly:

#. Configure QMCPACK by running CMake with ``CMAKE_EXPORT_COMPILE_COMMANDS=ON``, for example:

   ::

      cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++
      -DQMC_MPI=0 -DCMAKE_EXPORT_COMPILE_COMMANDS=ON ../qmcpack/

#. Run Understand and create a new C++ project. At the import files
   and settings dialog, import the ``compile_commands.json`` created by
   CMake in the build directory.
