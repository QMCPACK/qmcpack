# Change Log

Notable changes to QMCPACK will be documented in this file. 

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
 
