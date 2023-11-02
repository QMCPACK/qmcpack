# This file is configured by CMake automatically as DartConfiguration.tcl
# If you choose not to use CMake, this file may be hand configured, by
# filling in the required variables.


# Configuration directories and files
SourceDirectory: /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix
BuildDirectory: /gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix/build-complex-soa-fp-cee

# Where to place the cost data store
CostDataFile: 

# Site is something like machine.domain, i.e. pragmatic.crd
Site: cee-compute036

# Build name is osname-revision-compiler, i.e. Linux-2.4.2-2smp-c++
BuildName: Linux-mpic++

# Subprojects
LabelsForSubprojects: 

# Submission information
SubmitURL: http://cdash.qmcpack.org/CDash/submit.php?project=QMCPACK
SubmitInactivityTimeout: 

# Dashboard start time
NightlyStartTime: 00:00:00 EDT

# Commands for the build/test/submit cycle
ConfigureCommand: "/projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/cmake-3.23.1-v3pkd6okl3lzh723twjgmlfqs54zqrdv/bin/cmake" "/gpfs/rclay/qmcpack_worktree/complex_orb_opt_fix"
MakeCommand: /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/cmake-3.23.1-v3pkd6okl3lzh723twjgmlfqs54zqrdv/bin/cmake --build . --config "${CTEST_CONFIGURATION_TYPE}"
DefaultCTestConfigurationType: Release

# version control
UpdateVersionOnly: 

# CVS options
# Default is "-d -P -A"
CVSCommand: 
CVSUpdateOptions: 

# Subversion options
SVNCommand: 
SVNOptions: 
SVNUpdateOptions: 

# Git options
GITCommand: /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/git-2.35.2-tf3qm22shnekjfqgkzk435xf7tjvuj3g/bin/git
GITInitSubmodules: 
GITUpdateOptions: 
GITUpdateCustom: 

# Perforce options
P4Command: 
P4Client: 
P4Options: 
P4UpdateOptions: 
P4UpdateCustom: 

# Generic update command
UpdateCommand: /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-10.3.0/git-2.35.2-tf3qm22shnekjfqgkzk435xf7tjvuj3g/bin/git
UpdateOptions: 
UpdateType: git

# Compiler info
Compiler: /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/intel-2021.1.2/openmpi-4.1.2-awpydu2boorwxkfkfgtkgilngfsokjy4/bin/mpic++
CompilerVersion: 2021.1.0.20201208

# Dynamic analysis (MemCheck)
PurifyCommand: 
ValgrindCommand: 
ValgrindCommandOptions: 
DrMemoryCommand: 
DrMemoryCommandOptions: 
CudaSanitizerCommand: 
CudaSanitizerCommandOptions: 
MemoryCheckType: 
MemoryCheckSanitizerOptions: 
MemoryCheckCommand: /usr/bin/valgrind
MemoryCheckCommandOptions: 
MemoryCheckSuppressionFile: 

# Coverage
CoverageCommand: /projects/cde/v3/cee/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/gcc-10.3.0-bxx4ufppf7l3xlbvdinnusa4e6ykqok5/bin/gcov
CoverageExtraFlags: -l

# Testing options
# TimeOut is the amount of time in seconds to wait for processes
# to complete during testing.  After TimeOut seconds, the
# process will be summarily terminated.
# Currently set to 25 minutes
TimeOut: 3600

# During parallel testing CTest will not start a new test if doing
# so would cause the system load to exceed this value.
TestLoad: 

UseLaunchers: 
CurlOptions: 
# warning, if you add new options here that have to do with submit,
# you have to update cmCTestSubmitCommand.cxx

# For CTest submissions that timeout, these options
# specify behavior for retrying the submission
CTestSubmitRetryDelay: 5
CTestSubmitRetryCount: 3
