# Change Log

Notable changes to QMCPACK are documented in this file.

## [3.8.0] - 2019-07-23

### Notes

This release includes Quantum Espresso v6.4.1 support, new examples
for adding wavefunctions and Jastrow functions, and many updates to the
AFQMC code functionality. Additionally, all the updated scripts and
functionality utilized during the [2019 QMCPACK
workshop](https://github.com/QMCPACK/qmcpack_workshop_2019) are
provided; this link also includes several new tutorials. A large
number of feature refinements, bugfixes, testing improvements and
source code cleanup have also been performed.

* Quantum Espresso v6.4.1 support [\#1732](https://github.com/QMCPACK/qmcpack/pull/1732).
* New tutorial for adding a simple wavefunction (He) [\#1621](https://github.com/QMCPACK/qmcpack/pull/1621).
* New tutorial and capability for adding Jastrow functors from symbolic expressions written in Python [\#1557](https://github.com/QMCPACK/qmcpack/pull/1557).
* [Updated compiler and library support policy](https://github.com/QMCPACK/qmcpack#prerequisites), and matching testing. We aim to support open source compilers and
  libraries within two years of release. Use of older software is discouraged and untested. Support for closed source compilers over the same period may require use of an exact version.
* Many updates to AFQMC code to support more compilers and libraries.
* Newly expanded deterministic test set should now pass on all platforms and be usable as an.
  installation check. Recommend to run "ctest -L deterministic" after building QMCPACK.
* AFQMC code now only reads HDF5 format data to improve I/O performance and storage utilization.
* K-point AFQMC code usable in production (e.g. bug fix [\#1524](https://github.com/QMCPACK/qmcpack/pull/1524)).
* Updated AFQMC workflow scripts for interfacing with PySCF.
* Faster initial cusp correction calculation for all-electron calculations, e.g. [\#1643](https://github.com/QMCPACK/qmcpack/pull/1643).
* Improved stability of cusp correction calculation [\#1594](https://github.com/QMCPACK/qmcpack/pull/1594).
* New short-ranged e-n Jastrow [\#1680](https://github.com/QMCPACK/qmcpack/pull/1680).
* Substantially faster 1-body reduced density matrix (1DRDM) estimator [\#1672](https://github.com/QMCPACK/qmcpack/pull/1672).
* Performance tests added for LCAO code and Gaussian basis sets [\#1639](https://github.com/QMCPACK/qmcpack/pull/1639).
* Reduced configuration output by default. Use -DQMC_VERBOSE_CONFIGURATION=1 on CMake line for greater detail.
* Partial support for forces in LCAO basis e.g. [\#1559](https://github.com/QMCPACK/qmcpack/pull/1559). See details given at
  2019 QMCPACK Workshop and in manual.
* Improved human-readable Jastrow output [\#1525](https://github.com/QMCPACK/qmcpack/pull/1525).
* Improved MPI implementation. QMCPACK is now compatible with OpenMPI v4.
* Majority of the manual has been professionally edited.

### Known bugs

See [list of open bugs](https://github.com/QMCPACK/qmcpack/issues?q=is%3Aissue+is%3Aopen+label%3Abug).

* There is a general problem with MVAPICH2 involving aligned memory allocations that will cause
  QMCPACK to crash if MVAPICH is compiled using defaults. See  [\#1703](https://github.com/QMCPACK/qmcpack/issues/1703) for details and workaround.

### NEXUS

* Examples added of PySCF molecular and solid-state workflows [\#1552](https://github.com/QMCPACK/qmcpack/pull/1552).
* Update support for Quantum Package 2.0 [\#1538](https://github.com/QMCPACK/qmcpack/pull/1538).
* Support for additional machines including SuperMUC-NG [\#1665](https://github.com/QMCPACK/qmcpack/pull/1665).
* Support for ghost atoms [\#1653](https://github.com/QMCPACK/qmcpack/pull/1653).
* Update outdated cubic specifier to alat for QE [\#1642](https://github.com/QMCPACK/qmcpack/pull/1642).
* K-point grids are symmetrized with spglib [\#1544](https://github.com/QMCPACK/qmcpack/pull/1544).
* Bundling of jobs at NERSC [\#1748](https://github.com/QMCPACK/qmcpack/pull/1748).

## [3.7.0] - 2019-03-29

### Notes

This release includes GPU support for the AFQMC implementation,
Quantum Espresso v6.4 support, and in the real-space code makes the
structure-of-arrays (SoA) code path the default. A large number of
feature refinements, bugfixes, testing improvements and source code
cleanup have been performed.

* The improved structures of arrays (SoA) build is now the
  default. This is generally significantly faster and uses less memory
  than the AoS build due to better algorithms, but does not yet have
  the full range of functionality. The older AoS build can be selected
  with -DENABLE_SOA=0.
* AFQMC code fully supports GPU acceleration via NVIDIA CUDA. Use -DENABLE_CUDA=1.
* Quantum Espresso v6.4 is supported.  [\#1457](https://github.com/QMCPACK/qmcpack/pull/1457)
* Better error handling e.g.  [\#1423](https://github.com/QMCPACK/qmcpack/issues/1423)
* Workarounds for MPI support on Summit.  [\#1479](https://github.com/QMCPACK/qmcpack/pull/1479)
* ppconvert should be more reliable.  [\#891](https://github.com/QMCPACK/qmcpack/issues/891)
* Delayed update implementation on GPUs.  [\#1279](https://github.com/QMCPACK/qmcpack/pull/1279) 
* Continued improvements to the testing system and test coverage. While still under
  development, a new set of deterministic tests is intended to rapidly
  and reliably test the code, with good coverage. Tests pass for real
  and complex, but not yet mixed-precision or GPU builds.
* Source code has been formatted with clang-format for consistency throughout. 

### Known Bugs

See [list of open bugs](https://github.com/QMCPACK/qmcpack/issues?q=is%3Aissue+is%3Aopen+label%3Abug).

* Theres is a bug that could result in an incorrect local
  electron-ion pseudopotential energy with CUDA v9.1 and Kepler GPUs. This is still being
  investigated. [\#1440](https://github.com/QMCPACK/qmcpack/issues/1440)

* QMCPACK will not build with OpenMPI v4 due to use of deprecated
  functions. This will be addressed when the new MPI wrappers are
  fully adopted. Older OpenMPI libraries are fully capable.

### NEXUS

* A collection of training material is at https://github.com/QMCPACK/nexus_training
* Improved generation of QMCPACK inputs. [\#1471](https://github.com/QMCPACK/qmcpack/pull/1471)
* Improved Gaussian Process optimization. [\#1498](https://github.com/QMCPACK/qmcpack/pull/1498)
* Updated Cori support. [\#1463](https://github.com/QMCPACK/qmcpack/pull/1463)
* Supercell tiling is more robust. [\#1432](https://github.com/QMCPACK/qmcpack/pull/1432)
* Summit support. [\#1394](https://github.com/QMCPACK/qmcpack/pull/1394)
* Improved handling of excited state calculations. [\#1365](https://github.com/QMCPACK/qmcpack/pull/1365)
* Fixed CHGCAR conversion. [\#1351](https://github.com/QMCPACK/qmcpack/pull/1351)

## [3.6.0] - 2018-12-19

### Notes

This release includes a completely new AFQMC implementation,
significant performance improvements for large runs, greater
functionality in the structure-of-arrays (SoA) code path, support for
larger spline data on multiple GPUs, and support for new machines and
compilers. The manual has been improved, bugs have been fixed, and
source code cleanup continued.

A C++14 and C99 capable compiler, Boost 1.61.0, and CMake 3.6 or
greater are now required.

* Completely updated AFQMC implementation including reduced scaling separable density
  fitting https://arxiv.org/abs/1810.00284 Documentation and examples
  will be added in v3.7.0. Contact the developers for use instructions
  in the interim. [\#1245](https://github.com/QMCPACK/qmcpack/pull/1245)

* Implementation of delayed updates for CPU. Substantial speedups for
  runs with 100s of electrons, with increasing gains at larger
  electron counts. See manual for details.  [\#1170](https://github.com/QMCPACK/qmcpack/pull/1170)

* Initial support for nested OpenMP to further reduce time-to-solution for large problems. [\#1082](https://github.com/QMCPACK/qmcpack/pull/1082)

* Support for splitting/distributing spline orbital data across multiple GPUs on a
  single node. [\#1101](https://github.com/QMCPACK/qmcpack/pull/1101)

* Cusp correction for all electron calculations is implemented in the
  SoA version. [\#1172](https://github.com/QMCPACK/qmcpack/pull/1172)

* Backflow is implemented in the SoA version. [\#1225](https://github.com/QMCPACK/qmcpack/pull/1225)

* K-points with real coefficients are supported in periodic LCAO. [\#1006](https://github.com/QMCPACK/qmcpack/pull/1006)

* Initial support for Summit at OLCF. Revisions may be needed in
  January 2019 as the software stack is updated. This will be
  addressed in a new version as required.

* Initial support for PGI compiler.

* Build instructions for ARM-based systems.  [\#1148](https://github.com/QMCPACK/qmcpack/pull/1148)

* Setup scripts are python 2 and 3 compatible. [\#1261](https://github.com/QMCPACK/qmcpack/pull/1261)

* QMCPACK and NEXUS can now be installed by "make install" after
  configuring CMake with CMAKE_PREFIX_PATH. [\#1020](https://github.com/QMCPACK/qmcpack/issues/1020)

* Significantly reworked test labeling and categorization system. [\#1155](https://github.com/QMCPACK/qmcpack/pull/1155)

* Partial transition to a new MPI wrapper implementation for greater compatibility.

* Utilities have been renamed for clarity and to avoid name collisions
  with other applications. getSupercell is renamed
  qmc-get-supercell. extract-eshdf-kvectors is renamed
  qmc-extract-eshdf-kvectors.

### Known bugs

Several potentially significant bugs are outstanding and will be addressed in the
next release. See [list of open bugs](https://github.com/QMCPACK/qmcpack/issues?q=is%3Aissue+is%3Aopen+label%3Abug).

* LCAO (Gaussian basis) molecular calculations are incorrect with
  certain diffusion functions. The reason for this bug is currently
  unclear. [\#1145](https://github.com/QMCPACK/qmcpack/issues/1145)

* On NVIDIA Volta GPUs some runs show inconsistencies with the CPU
  version. Standard carbon diamond and LiH tests pass with good agreement
  with the CPU implementation. [\#1054](https://github.com/QMCPACK/qmcpack/issues/1054)

* QMCPACK will not build with OpenMPI v4.0.0 due to use of deprecated
  functions. This will be addressed in the next version as the new MPI
  wrappers are fully adopted. Older OpenMPI libraries are fully capable.


### NEXUS

* Interface to and support for PySCF. [\#1220](https://github.com/QMCPACK/qmcpack/pull/1220)
* Interface to and support for Quantum Package (QP). [\#1093](https://github.com/QMCPACK/qmcpack/pull/1093)
* Support for excited state calculations. [\#1200](https://github.com/QMCPACK/qmcpack/pull/1200)
* qfit is renamed qmc-fit.
* ntest, sim, redo are renamed nxs-test, nxs-sim, nxs-redo.
* Many smaller improvements.

## [3.5.0] - 2018-08-02

### Notes

This release includes support for the latest Quantum Espresso version 6.3,
an initial implementation of periodic Gaussian support via PySCF, and a new
version of the hybrid or "APW" representation of orbitals. Many minor
bugs have been fixed, configuration and documentation improved.

Note that the PDF manuals are no longer included with the source. Versions
are available online via https://qmcpack.org . The PDFs can be built
using manual/build_manual.sh and
nexus/documentation/user_guide_source/build_nexus_user_guide.sh

Attention developers: This version contains substantially fewer source
lines than previous versions due to clean out of old code and unused
execution paths. Refactoring to improve the internal structure of
QMCPACK is ongoing. Track the develop branch and follow discussion on
GitHub closely to avoid difficult merges.

* Support for Quantum Espresso 6.3 and 6.2.1. Check documentation to
  ensure compiled with required HDF5 support.
* Support for periodic gaussians and PySCF generated
  wavefunctions. Initial version is limited to Gamma-point.
* Improved hybrid representation of single particle orbitals
  (APW-like) for significantly reduced memory usage and possible
  accuracy increase compared to conventional spline representation. https://arxiv.org/abs/1805.07406
* Norms of orbitals are checked inside QMCPACK to catch conversion errors.
* Added verbosity setting to QMCPACK output.
* CUDA can now be enabled with SoA builds.
* Many improvements to QMCPACK manual, including all new features, CIPSI, 3-body
  jastrow factor description, spack package, and enabling HTML generation.
* CMake configuration improvements, particularly around MKL handling.
* Extensive cleanup of unused source files and unused code paths
  removed, reducing the number of source lines by over 30 percent.

### Known bugs

* Weight of first block of DMC density is incorrect in CPU
  code. DMC densities in CUDA GPU code are incorrect for all
  blocks. [\#934](https://github.com/QMCPACK/qmcpack/issues/934) and [\#925](https://github.com/QMCPACK/qmcpack/issues/925)
* Runs with only a single electron may crash. [\#945](https://github.com/QMCPACK/qmcpack/issues/945)

### NEXUS

* Support for GAMESS HDF5 workflows.
* Nexus accepts command line inputs.
* Nexus testing via ntest executable.
* Added GAMESS-NEXUS examples for RHF, CISD, and CASSCF wavefunction.
* Added support for -nojastrow workflows.
* Added support for Stampede supercomputer.
* Added script to build NEXUS user guide.
* Various bugfixes including to GAMESS input parsing.

## [3.4.0] - 2018-01-29

### Notes

This release includes size-consistent t-moves, and improvements to
load balancing and memory usage that will be visible in large
runs. Significant revisions have been made to the gaussian
wavefunction reader and a PySCF interface is in progress.  A bug
affecting non-git installs (from release tarballs) is fixed. Feedback
is particularly welcome on the new features.

* Size consistent t-moves implemented (Casula 2010 algorithm).
  Enabled via nonlocalmoves parameter, see manual.
* Bugfix: For non-git builds, build process failed on some systems
  due to git-rev.h handling.
* Optimized load balancing in DMC. Command line option async_swap
  removed. Parameter use_nonblocking now disables non-blocking MPI
  load balancing. Non-blocking MPI is now enabled by default.
* Improved memory handling and usage in SoA code, increases
  performance.
* Improved stability of GPU matrix inversion for large runs.
* Ongoing improvements to output to improve readability.
* Initial interface to PySCF for real space QMC trial wavefunctions.
* Enabled use of HDF5 files for Gaussian based wavefunctions
  with SoA implementation.
* Added Appendix to manual listing all known QMCPACK publications.
  This will be updated on an ongoing basis. Please advise of any
  missing publications.
* Optimized momentum distribution estimator. Supported by SoA and
  1,2,3-body Jastrow functions.
* Support for labeled timers in Intel VTune based profiling.

### NEXUS

* Minor bugfixes and improvements.

### Known limitations

* PySCF interface is preliminary. convert4qmc is updated, but manual
  entries are not yet provided. This will be improved in later
  versions. The interface is currently only for isolated molecular
  systems. A full periodic implementation is in progress.

* Documentation, examples and tutorials are not yet consistent with
  the updated converter convert4qmc.

## [3.3.0] - 2017-12-18

### Notes

This release includes new methods, converter updates, and many
optimizations, feature improvements, and bug fixes. It is a
recommended update for all users.

### QMCPACK updates

* Support for finite difference linear response (FDLR) method and
  wavefunctions, developed and contributed by Nick Blunt and Eric
  Neuscamman, see Journal of Chemical Physics 147, 194101 (2017),
  https://doi.org/10.1063/1.4998197 and
  https://arxiv.org/abs/1707.09439 .
* Major update to convert4qmc, conversion from GAMESS and other
  gaussian basis set codes. HDF5 output is now supported for large
  wavefunctions with -hdf5 option. Significantly improved example
  inputs \*.qmc.in.xml.
* Gaussian based trial wavefunctions now supported by structure of
  arrays implementation (ENABLE_SOA=1). A full reimplementation that
  will also support gaussians in periodic boundary conditions, e.g.
  from pyscf, is in progress.
* Initialization of multideterminant wavefunctions improved for faster
  startup and lower memory usage. In practice this significantly
  raises the usable maximum number of determinants.
* Maximum CPU time setting (maxcpusecs): QMC drivers will not
  start a new block if there is not enough estimated time remaining to
  complete the next block and gracefully shut down.
* Homogeneous electron gas wavefunction support and tests.
* New command line verbosity command line flag -verbosity. Output of
  QMCPACK will be overhauled over the next few releases to
  support low, high, and debug options, and also to significantly improve
  readability and utility.
* Bugfix: Umrigar drift diffusion term is now consistent with the
  Umrigar small time step error algorithm with complex wavefunctions.
* Bugfix: Momentum distribution is now correctly weighted and also
  correctly signed for twist averaging.
* Renamed performance tests with atom and electron count.
* Removed support for "buffering" of non-local pseudopotential
  wavefunction components during optimization (useBuffer setting) to
  reduce memory usage and for simplicity.
* doxygen documentation for developer-level documentation of the code
  and file structure. Produced via make in qmcpack/doxygen. HTML
  currently published at http://docs.qmcpack.org/doxygen/doxy/
* Many minor bug fixes and improved tests.

### NEXUS

* Improved postprocessing support for Quantum Espresso.
* Various minor bug fixes.

### Known issues and limitations

* Documentation, examples and tutorials are not yet consistent with
 the updated converter convert4qmc.
* Core functionality is largely compatible with ENABLE_SOA but
  some specialized wavefunctions and observables are not.
* Use of GNU compilers with glibc 2.23 builds will crash due to a bug
  in libmvec of glibc. The glibc version can be verified by
  "ldd --version".

## [3.2.0] - 2017-09-21

### Notes

This release provides a significant speed increase for
many calculations. A C++11 compiler is now required. It is a
recommended update.

### QMCPACK updates

* Major speedup for calculations using spline wavefunctions via
  initial implementation of "Structure of Arrays" data layout and
  improved algorithms. Enabled via -DENABLE_SOA=1. Benefits all CPU
  architectures. Many runs are doubled in speed. Not yet available for
  Gaussian-basis sets or for all observables and QMC methods. See
  writeup in manual for guidance.
* A compiler supporting C++11 is now required.
* DMC respects MaxCPUSecs parameter and will gracefully shut down and
  not start a new block if there is not sufficient estimated time to
  complete it.
* Checkpointing code rewritten for robustness and performance at scale.
  Parallel as well as serial HDF5 supported and autodetected.
* Improved beta-release of AFQMC code and documentation.
* Backflow documentation and optimization tips added.
* Correlated sampling VMC drivers reactivated.
* Added carbon graphite performance test similar to CORAL benchmark.
* Improvements to CMake and CTest usage.
* Build instructions for NERSC, ALCF, and OLCF machines updated.
* Latest manual PDF now available at http://docs.qmcpack.org

### NEXUS

* Significantly improved manual entry for "qmca" analysis tool, the
  main recommended tool for statistical analysis of QMCPACK data.
* Added time step fitting tool "qfit" for timestep extrapolation. Uses
  jack-knife statistical technique.
* Improved density file postprocessing.
* Support for Makov-Payne corrections.

## [3.1.1] - 2017-08-01

### Notes

This is a bugfix release and recommended update.

### QMCPACK updates

* Added numerical tolerance to check of jastrow cutoff and Wigner Seitz
  radius.
* CMake correctly configures when MPI is not present.
* Improved support for test coverage measurements.
* Added unit tests for some estimators.

### NEXUS

* IPython compatible exit handling (from Duy Le)


## [3.1.0] - 2017-06-21

### Notes

This release incorporates an improved DMC equilibration scheme,
numerous bugfixes, small improvements, and significantly improved
testing. It is a recommended update.

### QMCPACK updates

* Improved population control during DMC equilibration. Reduces variance on larger runs.
* Bugfix: Real valued wavefunction GPU code gave incorrect result for some non-gamma twists that could be made real, e.g. X point. Complex code (QMC_COMPLEX=1) was always correct.
* All particle move VMC and DMC algorithms enabled, tests added.
* Reptation Monte Carlo (RMC) enabled, tests added.
* Significantly improved AFQMC implementation.
* Added NiO based VMC and DMC performance tests and description in manual. Wavefunction files accessed via QMC_DATA.
* Added DMC tests with locality and t-moves approximations.
* Added AFQMC tests.
* Added test of real space QMC restart capabilities.
* Added tests for several estimators.
* Added unit test for DMC walker propagation, effective core potentials, and OhmmsPETE.
* To avoid filesystem limitations, QMC_SYMLINK_TEST_FILES can be set to symlink (1) or copy test files (0).
* Fixed mixed precision Ceperley force evaluation.
* Many updated tests to improve statistical reliability. Removed flux estimator from short tests because they were not reliable enough.
* Tests that rely on non-standard python modules that are not available are skipped.
* Error trap jastrow factors with cutoff radii larger than Wigner Seitz radius.
* Bugfix: Prevent users from adding correlation terms on non-existing electron pairs, e.g. up-down correlation terms when only up-spin particles are present.
* Support for measuring test coverage and performing coverage runs with cmake and ctest.
* Support for GCC7 and IBM XL (non Blue Gene) compiler.
* Support selecting GPU microarchitecture via -DCUDA_ARCH=sm_35(default).
* SummitDev IBM Minsky build recipe (Power8 + NVIDIA Pascal P100 GPUs).
* Significantly updated optimizer description in manual, including excited state optimization.
* Added description of using Intel MKL with non-Intel compilers in manual.
* Added description of MPIEXEC and MPIEXEC_NUMPROCS_FLAG to manual for systems where MPI runner is non-standard.
* Updated labs with correct pseudopotentials, basis set files.
* Many updated error messages and warnings.

### Known problems

* AFQMC without MKL will fail, e.g. short-afqmc-N2_vdz-4-1 test fails.

### NEXUS updates

* Improved selection algorithm to obtain optimally tiled supercells.
* Support for parallel pw2qmcpack workflows.
* Support for HPC resources at the Leibniz Supercomputing Center.
* Better consistency checks for the Structure class.
* Bugfix: forbid job bundling for simulations that depend on each other.
* Bugfix: correctly select low spin polarization in primitive and tiled (net_spin="low" option).


## [3.0.0] - 2017-01-30

### Notes

We are adopting [Semantic Versioning](http://semver.org) with this
release. It is the first to be made from the git repository on GitHub,
and the first named release since 2016-06-02 and subversion
revision 6964.

A potentially severe bug is fixed for periodic wavefunctions in this version,
in addition to many usability improvements and bugfixes. All users are
strongly recommended to upgrade.

NEXUS updates are listed after QMCPACK updates.

### QMCPACK updates

* IMPORTANT BUGFIX: Real-valued wavefunction code would occasionally make a numerically
  unstable choice for constructing real-valued periodic wavefunctions, leading to
  large variances and poor energies. Algorithm for constructing
  wavefunctions improved.
* Fully parallel pw2qmcpack.x for QE 5.3, enables conversion of large
  wavefunctions and use of same parallel setup as pw.x runs.
* Full testing of Quantum Espresso workflows (pw.x -> pw2qmcpack.x ->
  qmcpack). Specify directory containing QE binaries via QE_BIN during configuration.
* Added open boundary conditions tests using QE wavefunctions,
  as might be used for molecular work. Requires QE_BIN and computes
  trial wavefunction on the fly.
* Added DMC, optimizer and additional system tests.
* Added unit tests using the Catch framework.
* Plane wave wavefunctions can be evaluated in plane waves, use "pw"
  as determinantset type. Slow, but useful for checking spline accuracy. Tests added.
* Complex implementation on GPUs, supports arbitrary twists and
  complex phase wavefunctions as per CPU code.
* Flux estimator correct for complex wavefunctions.
* Mixed precision CPU implementation, activated via -DQMC_MIXED_PRECISION=1.
* Double precision GPU implementation, complementing existing
  mixed precision implementation, activated via -DQMC_MIXED_PRECISION=0.
* GAMESS CI converter improved.
* C++11 detection and support.
* Initial release of new optimizer, requires C++11 (contact Eric Neuscamman).
* Initial release of orbital-based AFQMC code, requires C+11 and MKL (contact Miguel Morales).
* Fine grained timers implemented, activated via -DENABLE_TIMERS=1.
* Improved Intel math and vector math library support. MKL and MKL VML more easily
  supported with GCC as well as Intel compilers.
* Many code updates to eliminate CLANG warnings.
* Configure scripts, printed headers, manual updated for git. Git
  version printed during configure and on standard output.
* Source files headers updated to consistently show UIUC/NCSA open source
  license and list development history.
* Numerous manual updates.
* Updated QMCPACK tutorial laboratories.
* Many small bug fixes, improvements and optimizations.

### NEXUS updates

* General
  *  Nexus output now tracks time instead of poll number.
  *  Reported memory use now includes child processes.
* Workflow generator
  *  Major new capability to generate simple to complex workflows involving QE, VASP, and QMCPACK.
  *  Aim is to allow single notebook/worksheet describing all simulation workflows needed in a project.
  *  Users can succinctly create any subchain of the workflow: relax->scf->nscf->orbital_conv->qmc.
  *  Additional elements can be added to workflow chains over time as needed.
  *  Scans of structural parameters and input parameters at any level of the chain are possible.
  *  No programming constructs are required (for/if, etc).
  *  Directory substructure is automatically generated in the case of scans.
  *  Native support for visualizing workflows via pydot is provided.
  *  Documentation for this feature is pending.
* Quantum Espresso workflows
  *  Support for vdW functional input.
  *  Fixes to SCF->NSCF workflows for QE 5.3.0+.
  *  Support for automatic restarts of SCF runs.
  *  Native support for workflows involving post-processing tools
    * pp.x, dos.x, bands.x, projwfc.x, cppp.x, pw_export.x supported.
    * Postprocessing and summary of Lowdin charge data from projwfc.x.
* QMCPACK workflows
  *  Fixes for QE/VASP structural relaxation -> QMCPACK workflows.
  *  Fixed job bundling of twist averaged runs.
  *  Support for partitioned sposet input.
* Supercomputing environments
  *  Native support for several supercomputing environments located at Sandia Nat. Labs.
* Atomic structure manipulation
  *  Ability to find optimal supercells, similar to getSupercell tool.
  *  Robustness fixes to tiling operations.
*  Tools
  *  qmca
    *  Fix for twist averaging with user-provided weights.
  *  qmcfit
    * New command line tool for jack-knife fitting of QMCPACK data.
    * Timestep extrapolation currently supported.
    * General binding/equation of state fitting pending.
