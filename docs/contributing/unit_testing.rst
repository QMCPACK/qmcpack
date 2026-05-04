.. _unit-testing:

Unit Testing
============

Unit testing is a standard software engineering practice to aid in ensuring a quality product. A good suite of unit tests provides confidence in refactoring and changing code, furnishes some documentation on how classes and functions are used, and can drive a more decoupled design.

If unit tests do not already exist for a section of code, you are encouraged to add them when modifying that section of code.  New code additions should also include unit tests.
When possible, fixes for specific bugs should also include a unit test that would have caught the bug.

Unit testing framework
----------------------

The Catch framework is used for unit testing.
See the project site for a tutorial and documentation: https://github.com/philsquared/Catch.

Catch consists solely of header files. It is distributed as a single include file about 400 KB in size.  In QMCPACK, it is stored in ``external_codes/catch``.

Unit test organization
----------------------

The source for the unit tests is located in the ``tests`` directory under each directory in ``src`` (e.g., ``src/QMCWavefunctions/tests``).
All of the tests in each ``tests`` directory get compiled into an executable.
After building the project, the individual unit test executables can be found in ``build/tests/bin``.
For example, the tests in ``src/QMCWavefunctions/tests`` are compiled into ``build/tests/bin/test_wavefunction``.

All the unit test executables are collected under ctest with the ``unit`` label.
When checking the whole code, it is useful to run through CMake (``cmake -L unit``).
When working on an individual directory, it is useful to run the individual executable.
One can work from one of the `tests` directories beneath the build directory for a faster build and test cycle.


Some of the tests reference input files. The unit test CMake setup places those input files in particular locations under the ``tests`` directory (e.g., ``tests/xml_test``).  The individual test needs to be run from that directory to find the expected input files.

Command line options are available on the unit test executables.  Some of the more useful ones are

-  ``-h`` List command line options.

-  ``-l`` List all the tests in the executable.

-  ``--turn-on-printout`` See all the QMCPACK printouts which are suppressed by default in unit tests.

A test name can be given on the command line to execute just that test.  This is useful when iterating
on a particular test or when running in the debugger.   Test names often contain spaces, so most command line environments require enclosing the test name in single or double quotes.

Example
-------

The first example is one test from ``src/Numerics/tests/test_grid_functor.cpp``.

.. code-block::
  :caption: Unit test example using Catch.
  :name: Listing 75

  TEST_CASE("double_1d_grid_functor", "[numerics]")
  {
    LinearGrid<double> grid;
    OneDimGridFunctor<double> f(&grid);

    grid.set(0.0, 1.0, 3);

    REQUIRE(grid.size() == 3);
    REQUIRE(grid.rmin() == 0.0);
    REQUIRE(grid.rmax() == 1.0);
    REQUIRE(grid.dh() == Approx(0.5));
    REQUIRE(grid.dr(1) == Approx(0.5));
  }

The test function declaration is
``TEST_CASE("double_1d_grid_functor","[numerics]")``.
The first argument is the test name, and it must be unique in the test suite.
The second argument is an optional list of tags.  Each tag is a name surrounded by brackets (``"[tag1][tag2]"``).  It can also be the empty string.

The ``REQUIRE`` macro accepts expressions with C++ comparison operators and records an error if the value of the expression is false.

Floating point numbers may have small differences due to roundoff, etc.   The ``Approx`` class adds some tolerance to the comparison.  Place it on either side of the comparison (e.g., ``Approx(a) == 0.3`` or ``a = Approx(0.3)``).   To adjust the tolerance, use the ``epsilon`` and ``scale`` methods to ``Approx`` (``REQUIRE(Approx(a).epsilon(0.001) = 0.3);``.

Expected output
~~~~~~~~~~~~~~~

When running the test executables individually, the output of a run with no failures should look like

::

  ===============================================================================
  All tests passed (26 assertions in 4 test cases)

A test with failures will look like

::

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  test_numerics is a Catch v1.4.0 host application.
  Run with -? for options

  -------------------------------------------------------------------------------
  double_1d_grid_functor
  -------------------------------------------------------------------------------
  /home/user/qmcpack/src/Numerics/tests/test_grid_functor.cpp:29
  ...............................................................................

  /home/user/qmcpack/src/Numerics/tests/test_grid_functor.cpp:39: FAILED:
    REQUIRE( grid.dh() == Approx(0.6) )
  with expansion:
    0.5 == Approx( 0.6 )

  ===============================================================================
  test cases:  4 |  3 passed | 1 failed
  assertions: 25 | 24 passed | 1 failed

Adding tests
------------

Three scenarios are covered here: adding a new test in an existing file, adding a new test file, and adding a new ``tests`` directory.

Adding a test to existing file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Copy an existing test or from the example shown here.  Be sure to change the test name.

Adding a test file
~~~~~~~~~~~~~~~~~~

When adding a new test file,
create a file in the test directory, or copy from an existing file.  Add the file name to the ``ADD_EXECUTABLE`` in the ``CMakeLists.txt`` file in that directory.
The pattern for the test file name is ``test_<ClassName>.cpp``.  Many older tests do not follow this pattern, but new tests should.

One (and only one) file must define the ``main`` function for the test executable by defining ``CATCH_CONFIG_MAIN`` before including the Catch header.  If more than one file defines this value, there will be linking errors about multiply defined values.

Some of the tests need to shut down MPI properly to avoid extraneous error messages. Those tests include ``Message/catch_mpi_main.hpp`` instead of defining ``CATCH_CONFIG_MAIN``.

Adding a test directory
~~~~~~~~~~~~~~~~~~~~~~~

Copy the ``CMakeLists.txt`` file from an existing ``tests`` directory.
Change the ``SRC_DIR`` name and the  files in the ``ADD_EXECUTABLES`` line.  The libraries to link in ``TARGET_LINK_LIBRARIES`` may need to be updated.

Add the new test directory to ``src/CMakeLists.txt`` in the ``BUILD_UNIT_TESTS`` section near the end.

Testing with random numbers
---------------------------

Many algorithms and parts of the code depend on random numbers, which makes validating the results difficult.
One solution is to verify that certain properties hold for any random number.
This approach is valuable at some levels of testing, but is unsatisfying at the unit test level.

The ``Utilities`` directory contains a "fake" random number generator that can be used for deterministic tests of these parts of the code.
Currently it outputs a single, fixed value every time it is called, but it could be expanded to produce more varied, but still deterministic, sequences.
See ``src/QMCDrivers/test_vmc.cpp`` for an example of using the fake random number generator.

Setting up objects
------------------
One of the more difficult parts of writing tests is constructing the object, and prerequisite objects.
There are three routes to building an object:

1. Construct the object directly.
2. Use an XML fragment and use the XML parsing paths to construct an object
3. For updated classes, construct the Input object and use that in the construction path.


Constructing the object directly can be the most difficult in terms of building all the prerequisite objects.
Building an object from an XML fragment has an advantage of being similar to input files.

Building from XML
-----------------

Use C++ raw string literals (strings delimited with ``R"(`` and ``)"``) to use XML fragments in the code.
The ``Libxml2Document`` class has a ``parseFromString`` function to parse XML input for testing.

The following code fragment to read the xml is common

.. code-block::

  const char* xml_str = R"(<tmp><wavefunction></wavefunction></tmp>)";
  Libxml2Document doc;
  bool okay = doc.parseFromString(xml_str);
  REQUIRE(okay);

After parsing, the ``Libxml2Document`` class has a ``getRoot`` function to the the root XML node.
The QMCPACK parsing functions often expect the tags they are parsing to be a child of the node
that is passed to the function.
For this case, put an additional tag as a parent of the target elements (The reason for ``<tmp>...</tmp>`` in the example above.

The ``Libxml2Document`` class can also read XML from a file with the ``parse`` function
Reading from a file can make the test code smaller, at the expense of maintaining an extra file.
