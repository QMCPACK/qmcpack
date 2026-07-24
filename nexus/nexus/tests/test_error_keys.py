from io import BytesIO, StringIO

import pytest

from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.ERROR_KEYS)


from ..error_keys import find_error_keys


def test_input_forms_and_return_lines(tmp_path):
    """Check supported input forms and optional matching-line collection."""
    text = (
        'QMCPACK execution begins\n'
        'APP_ABORT: invalid particleset reference\n'
        'MPI_ABORT was invoked on rank 0\n'
        )

    assert find_error_keys(text, qmcpack=True)
    assert find_error_keys(StringIO(text), qmcpack=True)

    output_file = tmp_path / 'app.err'
    output_file.write_text(text)
    assert find_error_keys(output_file, qmcpack=True)
    assert find_error_keys(str(output_file), qmcpack=True)

    found, lines = find_error_keys(
        text,
        mpi=True,
        qmcpack=True,
        return_lines=True,
        )
    assert found
    assert lines == [
        'APP_ABORT: invalid particleset reference',
        'MPI_ABORT was invoked on rank 0',
        ]

    with pytest.raises(TypeError):
        find_error_keys(BytesIO(text.encode()), qmcpack=True)


def test_selectors():
    """Check individual, batch, combined, and all-error selectors."""
    # Every group is opt-in, including operating-system diagnostics.
    assert find_error_keys('Segmentation fault (core dumped)', shell=True)
    assert not find_error_keys('calculation completed normally', shell=True)

    # Batch selectors enable their documented constituent sets.
    assert find_error_keys('UCX ERROR transport endpoint failed', hpc=True)
    assert find_error_keys(
        'terminate called after throwing an instance of std::runtime_error',
        code=True,
        )
    assert find_error_keys(
        'HDF5-DIAG: Error detected in HDF5',
        code_library=True,
        )
    assert find_error_keys(
        'scipy.sparse.linalg.ArpackNoConvergence: no convergence',
        python_module=True,
        )

    # Individual groups remain isolated unless selected together.
    assert not find_error_keys('APP_ABORT: invalid input', lapack=True)
    assert find_error_keys(
        'APP_ABORT: invalid input',
        lapack=True,
        qmcpack=True,
        )

    # all_errors includes application-specific sets.
    assert find_error_keys(
        'RMG Error: domain decomposition failed',
        all_errors=True,
        )

    with pytest.raises(ValueError):
        find_error_keys('anything')


def test_operating_system_catches():
    assert find_error_keys(
        'run.sh: line 8: 4217 Segmentation fault (core dumped)',
        shell=True,
        )
    assert find_error_keys('EDAC MC0: Hardware Error', shell=True)
    assert find_error_keys(
        'rank 3 terminated by SIGSEGV',
        linux_signals=True,
        )
    assert find_error_keys(
        'bash: pw.x: No such file or directory',
        posix=True,
        )
    assert find_error_keys(
        'fatal error: No space left on device',
        posix=True,
        )
    assert find_error_keys('I/O failed with errno ENOSPC', posix=True)


def test_hpc_environment_catches():
    assert find_error_keys(
        '[1712] UCX ERROR endpoint timed out',
        infiniband=True,
        )
    assert find_error_keys(
        'LNet peer unreachable: fatal transport error',
        lustre=True,
        )
    assert find_error_keys('GPFS: [ERROR] disk unavailable', gpfs=True)
    assert find_error_keys(
        'slurmstepd: error: Detected 1 oom-kill event',
        slurm=True,
        )
    assert find_error_keys(
        'JobState=OUT_OF_MEMORY State=OUT_OF_MEMORY',
        slurm=True,
        )
    assert find_error_keys(
        'resources_used.walltime=01:00:00 exit_status=137',
        pbs=True,
        )
    assert find_error_keys(
        'MPI_Abort was invoked on rank 2 in communicator MPI_COMM_WORLD',
        mpi=True,
        )
    assert find_error_keys(
        'libgomp: Thread creation failed: Resource unavailable',
        openmp=True,
        )


def test_compiled_code_catches():
    assert find_error_keys(
        "libstdc++.so: version `GLIBCXX_3.4.30' not found",
        linking=True,
        )
    assert find_error_keys(
        'forrtl: severe (174): SIGSEGV, segmentation fault occurred',
        fortran=True,
        )
    assert find_error_keys(
        'ERROR: AddressSanitizer: heap-use-after-free',
        cplusplus=True,
        )
    assert find_error_keys(
        'kernel launch returned cudaErrorIllegalAddress',
        cuda=True,
        )
    assert find_error_keys(
        'NCCL WARN Error: failed to extend /dev/shm/nccl-a1',
        cuda=True,
        )
    assert find_error_keys(
        'kernel returned hipErrorLaunchFailure',
        hip=True,
        )


def test_python_runtime_catches():
    assert find_error_keys('Traceback (most recent call last):', python=True)
    assert find_error_keys(
        'pyscf.lib.exceptions.SomeError: failed operation',
        python=True,
        )


def test_compiled_library_catches():
    assert find_error_keys(
        'Intel MKL FATAL ERROR: Cannot load libmkl_avx2.so',
        blas=True,
        )
    assert find_error_keys(
        'LAPACK computational failure: eigensolver failed',
        lapack=True,
        )
    assert find_error_keys(
        'HDF5-DIAG: Error detected in HDF5 (1.14.0) thread 0:',
        hdf5=True,
        )
    assert find_error_keys(
        'input.xml: parser error : Premature end of data',
        libxml2=True,
        )


def test_python_module_catches():
    assert find_error_keys(
        'numpy.linalg.LinAlgError: Singular matrix',
        numpy=True,
        )
    assert find_error_keys(
        'scipy.spatial.QhullError: Initial simplex is flat',
        scipy=True,
        )
    assert find_error_keys(
        'OSError: unable to open file (file signature not found)',
        h5py=True,
        )


def test_operating_system_near_misses():
    assert not find_error_keys(
        'The documentation discusses a segmentation fault example.',
        shell=True,
        )
    assert not find_error_keys(
        'Installed checkpoint handlers for SIGTERM and SIGUSR1',
        linux_signals=True,
        )
    assert not find_error_keys(
        'EINTR and EAGAIN are retried by the I/O loop',
        posix=True,
        )
    assert not find_error_keys(
        'Optional cache file does not exist; continuing normally',
        posix=True,
        )


def test_hpc_environment_near_misses():
    assert not find_error_keys('FI_PROVIDER=FI_VERBS;FI_PSM2', infiniband=True)
    assert not find_error_keys('LNet: initialized 8 peer interfaces', lustre=True)
    assert not find_error_keys('State=COMPLETED ExitCode=0:0', slurm=True)
    assert not find_error_keys(
        'resources_used.walltime=00:02:14 exit_status=0',
        pbs=True,
        )
    assert not find_error_keys(
        'MPI_ERR_LASTCODE is the upper bound for implementation codes',
        mpi=True,
        )


def test_compiled_code_near_misses():
    assert not find_error_keys(
        'cudaEventQuery returned cudaErrorNotReady',
        cuda=True,
        )
    assert not find_error_keys(
        'ncclCommGetAsyncError returned ncclInProgress',
        cuda=True,
        )
    assert not find_error_keys(
        'hipEventQuery returned hipErrorNotReady',
        hip=True,
        )


def test_python_runtime_near_misses():
    assert not find_error_keys(
        'except ValueError: use the default input',
        python=True,
        )
    assert not find_error_keys(
        'The RuntimeError class is translated internally',
        python=True,
        )


def test_compiled_library_near_misses():
    assert not find_error_keys(
        'fftw_execute(plan); imported wisdom unavailable, replanning',
        fftw=True,
        )
    assert not find_error_keys(
        'SCF iteration 4 failed to converge; reducing mixing',
        lapack=True,
        )
    assert not find_error_keys(
        'capacity exceeded the estimate, resizing buffer',
        hdf5=True,
        )
    assert not find_error_keys(
        'Validation failed: optional metadata ignored',
        libxml2=True,
        )


def test_python_module_near_misses():
    assert not find_error_keys(
        'FileExistsError is caught before creating the HDF5 file',
        h5py=True,
        )


def test_pwscf_output():
    assert find_error_keys(
        '''
     iteration # 100     ecut= 50.00 Ry
     estimated scf accuracy < 1.0E-3 Ry
     convergence NOT achieved after 100 iterations: stopping
''',
        pwscf=True,
        )
    assert find_error_keys(
        '''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Error in routine cdiaghg (865):
     problems computing cholesky
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
''',
        pwscf=True,
        )
    assert not find_error_keys(
        '''
     iteration # 8
     convergence has been achieved in 8 iterations
     JOB DONE.
''',
        pwscf=True,
        )


def test_pyscf_output():
    assert find_error_keys(
        '''
cycle= 50 E= -75.983948  delta_E= -1.2e-7  |g|= 0.00031
SCF not converged.
SCF energy = -75.983948
''',
        pyscf=True,
        )
    assert not find_error_keys(
        '''
converged SCF energy = -75.983948
<S^2> = 0  2S+1 = 1
''',
        pyscf=True,
        )


def test_quantum_package_output():
    assert find_error_keys(
        '''
Summary at N_det = 1000000
qp run: Error: Selection failed before the requested PT2 threshold
''',
        quantum_package=True,
        )
    assert not find_error_keys(
        '''
Selection completed
Davidson diagonalization converged
N_det = 42560
''',
        quantum_package=True,
        )


def test_rmg_output():
    assert find_error_keys(
        '''
RMGDFT: starting domain decomposition
RMG Error: domain decomposition failed for the requested processor grid
''',
        rmg=True,
        )
    assert not find_error_keys(
        '''
RMGDFT: SCF convergence achieved
Total energy = -15.231901 Ha
''',
        rmg=True,
        )


def test_qmcpack_output():
    assert find_error_keys(
        '''
QMCPACK 4.1.0
Reading particlesets from input.xml
APP_ABORT: inconsistent input settings in determinantset
''',
        qmcpack=True,
        )
    assert not find_error_keys(
        '''
QMCPACK execution completed
Total Execution time = 12.4 sec
''',
        qmcpack=True,
        )


def test_vasp_output():
    assert find_error_keys(
        '''
 -----------------------------------------------------------------------------
|     VERY BAD NEWS! internal error in subroutine IBZKPT:                    |
|     Reciprocal lattice and k-lattice belong to different class of lattices |
 -----------------------------------------------------------------------------
''',
        vasp=True,
        )
    assert find_error_keys(
        '''
  25 F= -.78260631E+03 E0= -.78254989E+03 d E =-.117534E-05
  ZBRENT: fatal error in bracketing
  please rerun with smaller EDIFF, or copy CONTCAR to POSCAR and continue
''',
        vasp=True,
        )
    assert find_error_keys(
        '''
DAV:   4    -0.124848E+03   -0.34210E-03
EDDDAV: Call to ZHEGV failed. Returncode = 6 3 8
''',
        vasp=True,
        )
    assert not find_error_keys(
        '''
WARNING: small aliasing (wrap around) errors must be expected
reached required accuracy - stopping structural energy minimisation
General timing and accounting informations for this job:
''',
        vasp=True,
        )


def test_gamess_output():
    assert find_error_keys(
        '''
 SCF IS UNCONVERGED, TOO MANY ITERATIONS
 FINAL UHF ENERGY IS 0.0000000000 AFTER 30 ITERATIONS
 EXECUTION OF GAMESS TERMINATED -ABNORMALLY- AT Thu Nov 14 13:08:57 2019
''',
        gamess=True,
        )
    assert find_error_keys(
        '''
 ***** ERROR: MEMORY REQUEST EXCEEDS AVAILABLE MEMORY
 PROCESS NO. 0 WORDS REQUIRED= 300899439 AVAILABLE= 250000000
 EXECUTION OF GAMESS TERMINATED -ABNORMALLY-
''',
        gamess=True,
        )
    assert find_error_keys(
        '''
 DDI Process 0: error code 911
 ddikick.x: application process 0 quit unexpectedly.
 ddikick.x: Execution terminated due to error(s).
''',
        gamess=True,
        )
    assert not find_error_keys(
        '''
          FINAL RHF ENERGY IS      -75.9839481234
 EXECUTION OF GAMESS TERMINATED NORMALLY
 TOTAL WALL CLOCK TIME= 1.2 SECONDS, CPU UTILIZATION IS 99.4%
''',
        gamess=True,
        )
