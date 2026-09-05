.. _integration_tests:

Integration Tests
=================

Unlike unit tests requiring only a specific part of QMCPACK being built for testing, integration tests require the qmcpack
executable. In this category, tests are made based on realistic simulations although the amount of statistics collected depends on
sub-categories:

* Deterministic integration tests must be 100% reliable, quick to run, and always pass. They usually run one or a few walkers for
  a very few steps in a few seconds. They are used to rapidly identify changes as part of the continuous integration testing, to
  verify installations, and for development work.
* Short integration tests mostly run 16 walkers in a few hundred steps within a minutes. These are usually stochastic and should
  pass with very high reliability.
* Long integration tests mostly run 16 walkers in a few thousand steps within 10 minutes. These are usually stochastic and should
  pass with very high reliability.

To keep overall testing costs down, electron counts are usually kept small while still being large enough to comprehensively test
the code e.g. 3-10. The complete test set except for the long tests has to be able to be run on a laptop or modest workstation in
a reasonable amount of time.

Integration test organization
-----------------------------

Integration tests are placed under directories such as ``tests/heg``, ``tests/solids`` and ``tests/molecules`` from the top
directory and one sub-directory for each simulation system. Each test source directory contains input XML files, orbital h5 files,
pseudo-potential files and reference data (qmc_ref). These files may be shared by a few tests to minimize duplicated files. When
cmake is invoked in the build directory, one directory per test is created and necessary files correspond to a given test are
softlinked. It serves as a working directory when that test is being executed. To minimize the number file operation and make the
cmake execution fast, there is limitation on file names used by tests. The filenames are given below and implemented in the
COPY_DIRECTORY_USING_SYMLINK_LIMITED function in Cmake/macros.cmake.

::

  qmc-ref/qmc_ref for reference data folder.
  *.opt.xml/*.ncpp.xml/*.BFD.xml/*.ccECP.xml for pseudo-potential files.
  *.py/*.sh for result checking helper scripts.
  *.wfj.xml/*.wfnoj.xml/*.wfs.xml for standalone wavefunction input files.
  *.structure.xml/*.ptcl.xml for standalone structure/particleset input files.

How to add a integration test
-----------------------------

#. Generate reference data using a very long (many blocks >=2000) and possibly wide run (many nodes). This reduces both the
   error bar and the error bar of the error bar (10x samples than long test, 100x samples than short test). A folder named qmc-ref
   containing input.xml, scalar.dat and output file is required with the commit. The number of blocks should be about 200 to avoid
   large text files (a simple way to obtain these files is to repeat the reference run with 10x fewer blocks and 10x more steps).
#. Generate the short/long run input files. Use the reference error bar to appropriately estimate the error bar for the long and
   short tests. These error bars are sqrt(10+1) and sqrt(100+1) times larger than the very long reference. 10x grade is not a hard
   requirement but ref >= 10 long, long >= 10 short are required.
#. Short tests must be less than 20 sec VMC, 1 min OPT/DMC on a 16core Xeon processor. Long tests are preferably in the 5-10min
   range. For systems containing more than just a few electrons submitting only a long test may be appropriate.
#. Deterministic tests require a different approach: use of a fixed seed value, and for example, 3 blocks of 2 steps and a single
   walker. The intent of these tests is to exercise the code paths but keep the run short enough that the numerical deviations do
   not build up. Different reference data may be needed for mixed precision vs full precision runs.

Suggested procedure to add a test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#. Study some of the existing tests and their CMakeLists.txt configuration file to see what is required and the typical system
   sizes and run lengths used.
#. Perform a run ~30s in length on a 16 core machine (200 blocks min) using the CPU version of the code with 16 MPI and 1 thread
   per MPI. Decide if the resulting error bar is meaningful for a test. If so, short and long tests should be created. If not,
   possibly only a long test is appropriate.
#. Perform a reference run by increasing steps and blocks by 10x each (2000 blocks) and obtain reference mean and error bars.
   Long and short test error bars are then sqrt(100+1) and sqrt(10+1) of the reference.
#. Generate reference scalar data by redoing the reference run with 200 blocks and 100x steps. These data are should be committed
   in a qmc-ref directory with the test.
#. Create short (1x blocks, 1x steps) and long (1x blocks, 10x steps) input files (200 blocks each). Make one set of input files
   for CPU runs (walkers=1) and another for GPU runs (walkers=16).
#. Create CMakeLists.txt by following the example in other tests. CPU runs should include at least a 4 MPI, 4 thread test since
   this tests OpenMP, MPI, and any possible interactions between them. A GPU test should have 1 MPI and 16 threads.
#. Create a README file with information describing the tests and the reference data.
#. Check that the tests run properly with ctest on your local machine.
#. Submit a pull request with the final tests.

