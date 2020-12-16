.. _integration_tests:

Integration Tests
=================

Unlike unit tests requiring only a specific part of QMCPACK being built for testing, integration tests require the qmcpack executable.
In this category, tests are made based on realistic simulations although the amount of statistics collected depends on sub-categories:

- Deterministic integration tests run one or a few walkers in a few steps within 10 seconds.
- Short integration tests mostly run 16 walkers in a few hundred steps within a minutes.
- Long integration tests mostly run 16 walkers in a few thousand steps within 10 minutes.

Integration test organization
-----------------------------

Integration tests are placed under ``tests/heg``, ``tests/solids`` and ``tests/molecules`` from the top directory and one sub-directory for each simulation system.
Each test source directory contains input XML files, orbital h5 files, pseudo-potential files and reference data (qmc_ref).
These files may be shared by a few tests to minimize duplicated files.
When cmake is invoked in the build directory, one directory per test is created and necessary files correspond to a given test are softlinked.
It serves as a working directory when that test is being executed. To minimize the number file operation and make the cmake execution fast,
there is limitation on file names used by tests.

::

  qmc-ref/qmc_ref for reference data folder.
  *.opt.xml/*.ncpp.xml/*.BFD.xml for pseudo-potential files.
  *.py/*.sh for result checking helper scripts.
  *.wfj.xml/*.wfnoj.xml/*.wfs.xml for standalone wavefunction input files.
  *.structure.xml/*.ptcl.xml for standalone structure/particleset input files.

How to add a integration test
-----------------------------

#. Have a very long run (many blocks >=2000) and possibly wide (many nodes) to generate reference data. This reduces both the error bar and the error bar of the error bar (10x samples than long test, 100x samples than short test). A folder named qmc-ref containing input.xml, scalar.dat and output file is required with the commit. The number of blocks should be about 200 to avoid large text files (a simple way to obtain these files is to redo the reference run with 10x fewer blocks and 10x more steps).
#. Generate the short/long run input files. Use the reference error bar to appropriately estimate the error bar for the long and short tests. These error bars are sqrt(10+1) and sqrt(100+1) times larger than the very long reference. 10x grade is not a hard requirement but ref >= 10 long, long >= 10 short are required.
#. Short tests must be less than 20 sec VMC, 1 min OPT/DMC on a 16 core Xeon processor. Long tests are preferably in the 5-10min range. For systems containing more than just a few electrons submitting only a long test may be appropriate.

A convenient way to proceed is as follows:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#. Perform a run ~30s in length on a 16 core machine (200 blocks min) using the CPU version of the code with 16 MPI and 1 thread per MPI. Decide if the resulting error bar is meaningful for a test. If so, short and long tests should be created. If not, possibly only a long test is appropriate.
#. Perform a reference run by increasing steps and blocks by 10x each (2000 blocks) and obtain reference mean and error bars. Long and short test error bars are then sqrt(100+1) and sqrt(10+1) of the reference.
#. Generate reference scalar data by redoing the reference run with 200 blocks and 100x steps. These data are should be committed in a qmc-ref directory with the test.
#. Create short (1x blocks, 1x steps) and long (1x blocks, 10x steps) input files (200 blocks each). Make one set of input files for CPU runs (walkers=1) and another for GPU runs (walkers=16).
#. Create CMakeLists.txt by following the example in other tests. CPU runs should include at least a 4 MPI, 4 thread test since this tests OpenMP, MPI, and any possible interactions between them. A GPU test should have 1 MPI and 16 threads.
#. Create a README file with information describing the tests and the reference data.
#. Check that the tests run properly with ctest on your local machine.
#. Submit a pull request with the final tests.

