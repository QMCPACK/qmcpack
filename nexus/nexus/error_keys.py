"""Error diagnostics commonly printed by scientific applications.

The tuples ending in ``_errors`` contain readable examples of diagnostics.
The tuples ending in ``_error_patterns`` contain regular expressions used to
match variable portions or to add context that reduces false positives.
"""

import os
import re
from functools import cache


# Operating-system errors

shell_errors = (
    'Segmentation fault',
    'Floating point exception',
    'Illegal instruction',
    'Bus error',
    'Bad system call',
    'Aborted',
    'Killed',
    'Terminated',
    'Stack overflow',
    'Out of memory',
    'Cannot allocate memory',
    # Kernel OOM messages can identify a different process on the node.
    # 'oom-kill',
    # 'invoked oom-killer',
    # 'Killed process 123',
    'run.sh: line 8: 4217 Killed',
    # Machine-check and EDAC diagnostics can report corrected hardware errors.
    # 'Machine Check Exception',
    # 'EDAC Hardware Error',
    'stack smashing detected',
    'general protection fault',
    )

shell_error_patterns = (
    r'^.*\b(?:segmentation fault|floating point exception|illegal instruction|bus error|bad system call)(?:\s+\(core dumped\))?\s*$',
    r'^\s*(?:aborted|killed)(?:\s+\(core dumped\))?\s*$',
    r'^.*:\s+line\s+\d+:\s+\d+\s+(?:aborted|killed)(?:\s+\(core dumped\))?(?:\s+.+)?$',
    r'^\s*terminated\s*$',
    r'^\s*(?:out of memory|cannot allocate memory)\s*$',
    r'\b(?:stack overflow|stack smashing detected|general protection fault)\b',
    # Machine-check and EDAC diagnostics can report corrected hardware errors.
    # r'\bMachine Check Exception\b',
    # r'\b(?:MCE|EDAC)[^\n]{0,80}\bHardware Error\b',
    )

# Signal names are matched only when termination context is present.  Several
# other POSIX signals are routinely used for job control and checkpointing.
linux_exit_signals = (
    'SIGHUP',
    'SIGILL',
    'SIGABRT',
    'SIGFPE',
    'SIGKILL',
    'SIGSEGV',
    'SIGPIPE',
    'SIGTERM',
    'SIGBUS',
    'SIGSYS',
    'SIGTRAP',
    'SIGXCPU',
    'SIGXFSZ',
    'terminated with signal 11',
    'exited on signal 11',
    )

linux_signal_error_patterns = (
    r'\b(?:terminated|killed|exited|aborted|died|received signal)\b[^\n]*\bSIG(?:HUP|ILL|ABRT|FPE|KILL|SEGV|PIPE|TERM|BUS|SYS|TRAP|XCPU|XFSZ)\b',
    r'\bSIG(?:HUP|ILL|ABRT|FPE|KILL|SEGV|PIPE|TERM|BUS|SYS|TRAP|XCPU|XFSZ)\b[^\n]*\b(?:terminated|killed|exited|aborted|died)\b',
    r'\bterminated with signal\s+\d+\b',
    r'\bexited on signal\s+\d+\b',
    )

posix_errors = (
    # Most errno messages can result from handled probes or retryable I/O.
    # Missing executables and explicitly fatal errors provide run-failure
    # context rather than relying on the errno message alone.
    'bash: pw.x: No such file or directory',
    'bash: ./pw.x: Permission denied',
    'fatal error: No such file or directory',
    'fatal error: Permission denied',
    'fatal error: Not a directory',
    'fatal error: Is a directory',
    'fatal error: No space left on device',
    'fatal error: Too many open files',
    'fatal error: Cannot allocate memory',
    'fatal error: Connection refused',
    'fatal error: Connection timed out',
    'fatal error: Network is unreachable',
    'fatal error: Address already in use',
    'fatal error: Broken pipe',
    'fatal error: errno ENOSPC',
    )

posix_errno_keys = (
    'ENOENT',
    'EACCES',
    'EISDIR',
    'ENOTDIR',
    'ENOSPC',
    'EMFILE',
    'ENOMEM',
    'ECONNREFUSED',
    'ETIMEDOUT',
    'ENETUNREACH',
    'EADDRINUSE',
    'EPIPE',
    'EIO',
    'ENXIO',
    'EBADF',
    'EBUSY',
    'ENODEV',
    'EROFS',
    'EDQUOT',
    'ECONNRESET',
    'EHOSTUNREACH',
    'ENOTCONN',
    )

posix_error_patterns = (
    r'^\s*(?:bash|sh|zsh|ksh):(?:\s+line\s+\d+:)?\s+[^:\n]+:\s+(?:no such file or directory|permission denied)\s*$',
    r'\bfatal(?:\s+error)?[^\n]*\b(?:no such file or directory|permission denied|not a directory|is a directory|no space left on device|too many open files|cannot allocate memory|connection refused|connection timed out|network is unreachable|address already in use|broken pipe)\b',
    r'\bfatal[^\n]*\b(?:errno\s+)?(?:ENOENT|EACCES|EISDIR|ENOTDIR|ENOSPC|EMFILE|ENOMEM|ECONNREFUSED|ETIMEDOUT|ENETUNREACH|EADDRINUSE|EPIPE|EIO|ENXIO|EBADF|EBUSY|ENODEV|EROFS|EDQUOT|ECONNRESET|EHOSTUNREACH|ENOTCONN)\b',
    )


# HPC environment errors

infiniband_errors = (
    # Provider initialization and endpoint errors can be retried or cause a
    # fallback to another transport; none alone establishes run failure.
    # 'UCX ERROR',
    # 'libibverbs: malformed packet',
    # 'ucp_ep_create failed',
    # 'uct_ep_connect_to_ep unreachable',
    # 'ucs_init timed out',
    # 'ibv_create_qp failed',
    # 'libfabric transport error',
    # 'ofi_endpoint unreachable',
    )

infiniband_error_patterns = (
    # See infiniband_errors: these require separate termination context.
    )

lustre_errors = (
    # Filesystem and network-layer errors can concern another client or be
    # recovered by retry/failover without invalidating this simulation.
    # 'LustreError:',
    # 'LBUG:',
    # 'LNetError transport failed',
    # 'Lustre client evicted',
    )

lustre_error_patterns = (
    # See lustre_errors: these require separate termination context.
    )

gpfs_errors = (
    # GPFS can retry, renew tokens, or fail over disks; these messages can also
    # describe node-wide events unrelated to the process being inspected.
    # 'GPFS: [ERROR]',
    # 'GPFS: [FATAL]',
    # 'GPFS deadlock detected',
    # 'GPFS disk unavailable',
    # 'GPFS unmounted abnormally',
    # 'GPFS token expired',
    )

gpfs_error_patterns = (
    # See gpfs_errors: these require separate termination context.
    )

slurm_errors = (
    # Generic launcher errors can describe non-fatal setup/cleanup issues.
    # 'slurmstepd: error:',
    # 'srun: error:',
    'srun: Force term; sending SIGKILL',
    'DUE TO TIME LIMIT',
    'Exceeded job memory limit',
    'State=FAILED',
    'State=TIMEOUT',
    'State=NODE_FAIL',
    'State=OUT_OF_MEMORY',
    'State=BOOT_FAIL',
    'State=DEADLINE',
    'State=CANCELLED',
    'State=PREEMPTED',
    'JOB CANCELLED',
    'STEP FAILED',
    'srun: launch failed',
    'slurmstepd: error: Detected 1 oom-kill event',
    )

slurm_error_patterns = (
    r'\bsrun:\s*Force term;\s*sending SIGKILL\b',
    r'\bDUE TO TIME LIMIT\b',
    r'\bExceeded job memory limit\b',
    r'\bState=(?:FAILED|TIMEOUT|NODE_FAIL|OUT_OF_MEMORY|BOOT_FAIL|DEADLINE|CANCELLED|PREEMPTED)\b',
    r'\b(?:JOB|STEP)[^\n]*\b(?:CANCELLED|FAILED|OUT_OF_MEMORY|TIMEOUT|NODE_FAIL)\b',
    r'\b(?:srun|slurmstepd):[^\n]*\blaunch failed\b',
    r'\bslurmstepd:[^\n]*\boom-kill\b',
    )

pbs_errors = (
    'PBS: job killed:',
    'ob_init: Unable to read server database',
    'cannot send job to mom',
    'qsub: Bad UID for job execution',
    'exit_status=1',
    )

pbs_error_patterns = (
    r'\bPBS:\s*job killed:',
    r'\bob_init:\s*Unable to read server database\b',
    r'\bcannot send job to mom\b',
    r'\bqsub:\s*Bad UID for job execution\b',
    r'\bexit_status\s*=\s*(?!0\b)-?\d+\b',
    )

mpi_errors = (
    # MPI error-class names are return values that applications may handle.
    # MPI_Abort is also an API name; termination context is required below.
    # 'MPI_ABORT',
    # 'MPI_ERR_BUFFER',
    # 'MPI_ERR_COUNT',
    # 'MPI_ERR_TYPE',
    # 'MPI_ERR_TAG',
    # 'MPI_ERR_COMM',
    # 'MPI_ERR_RANK',
    # 'MPI_ERR_REQUEST',
    # 'MPI_ERR_ROOT',
    # 'MPI_ERR_GROUP',
    # 'MPI_ERR_OP',
    # 'MPI_ERR_TOPOLOGY',
    # 'MPI_ERR_DIMS',
    # 'MPI_ERR_ARG',
    # 'MPI_ERR_UNKNOWN',
    # 'MPI_ERR_TRUNCATE',
    # 'MPI_ERR_OTHER',
    # 'MPI_ERR_INTERN',
    # 'MPI_ERR_IN_STATUS',
    # 'MPI_ERR_PENDING',
    # 'MPI_ERR_ACCESS',
    # 'MPI_ERR_AMODE',
    # 'MPI_ERR_ASSERT',
    # 'MPI_ERR_BAD_FILE',
    # 'MPI_ERR_BASE',
    # 'MPI_ERR_CONVERSION',
    # 'MPI_ERR_DISP',
    # 'MPI_ERR_DUP_DATAREP',
    # 'MPI_ERR_FILE_EXISTS',
    # 'MPI_ERR_FILE_IN_USE',
    # 'MPI_ERR_FILE',
    # 'MPI_ERR_INFO_KEY',
    # 'MPI_ERR_INFO_NOKEY',
    # 'MPI_ERR_INFO_VALUE',
    # 'MPI_ERR_INFO',
    # 'MPI_ERR_IO',
    # 'MPI_ERR_KEYVAL',
    # 'MPI_ERR_LOCKTYPE',
    # 'MPI_ERR_NAME',
    # 'MPI_ERR_NO_MEM',
    # 'MPI_ERR_NOT_SAME',
    # 'MPI_ERR_NO_SPACE',
    # 'MPI_ERR_NO_SUCH_FILE',
    # 'MPI_ERR_PORT',
    # 'MPI_ERR_PROC_ABORTED',
    # 'MPI_ERR_QUOTA',
    # 'MPI_ERR_READ_ONLY',
    # 'MPI_ERR_RMA_CONFLICT',
    # 'MPI_ERR_RMA_SYNC',
    # 'MPI_ERR_SERVICE',
    # 'MPI_ERR_SIZE',
    # 'MPI_ERR_SPAWN',
    # 'MPI_ERR_UNSUPPORTED_DATAREP',
    # 'MPI_ERR_UNSUPPORTED_OPERATION',
    # 'MPI_ERR_WIN',
    # 'MPI_T_ERR_MEMORY',
    # 'MPI_T_ERR_NOT_INITIALIZED',
    # 'MPI_T_ERR_CANNOT_INIT',
    # 'MPI_T_ERR_INVALID_INDEX',
    # 'MPI_T_ERR_INVALID_ITEM',
    # 'MPI_T_ERR_INVALID_HANDLE',
    # 'MPI_T_ERR_OUT_OF_HANDLES',
    # 'MPI_T_ERR_OUT_OF_SESSIONS',
    # 'MPI_T_ERR_INVALID_SESSION',
    # 'MPI_T_ERR_CVAR_SET_NOT_NOW',
    # 'MPI_T_ERR_CVAR_SET_NEVER',
    # 'MPI_T_ERR_PVAR_NO_STARTSTOP',
    # 'MPI_T_ERR_PVAR_NO_WRITE',
    # 'MPI_T_ERR_PVAR_NO_ATOMIC',
    # 'MPI_ERR_RMA_RANGE',
    # 'MPI_ERR_RMA_ATTACH',
    # 'MPI_ERR_RMA_FLAVOR',
    # 'MPI_ERR_RMA_SHARED',
    # 'MPI_T_ERR_INVALID',
    # 'MPI_T_ERR_INVALID_NAME',
    # 'MPI_ERR_SESSION',
    'mpirun: kill job',
    'orterun: kill job',
    'prterun: kill job',
    'mpirun noticed that process rank 1 exited on signal 11',
    # ORTE/PRTE can log handled internal errors during component selection.
    # 'ORTE_ERROR_LOG',
    # 'PRTE_ERROR_LOG',
    'the first job to fail is listed below',
    'job aborted:',
    'MPI_Abort was invoked',
    'one or more processes exited with non-zero status',
    'process returned a non-zero exit code',
    'Primary job terminated normally, but',
    'mpirun aborted',
    'mpiexec failed to launch the application',
    'orterun terminated',
    'prterun exited on signal 11',
    'exited on signal 11',
    'terminated with signal 11',
    # A callback name and generic cleanup text do not establish failure.
    # 'mpiexec_callback_proc',
    # 'cleaning up processes',
    'execvp error',
    )

mpi_error_patterns = (
    r'\b(?:mpirun|orterun|prterun):\s*kill job\b',
    r'\bmpirun noticed that process rank\s+\d+[^\n]*\b(?:non-zero|signal|terminated|aborted|died)\b',
    r'\bthe first job to fail is listed below\b',
    r'\bjob aborted:',
    r'\bMPI_Abort was invoked\b',
    r'\bone or more processes exited with non-zero status\b',
    r'\bprocess returned a non-zero exit code\b',
    r'\bPrimary job terminated normally, but\b',
    r'\b(?:mpirun|mpiexec|orterun|prterun)\b[^\n]*\b(?:aborted|failed|non-zero|signal|terminated)\b',
    r'\b(?:exited on|terminated with) signal(?:\s+\d+)?\b',
    r'\bexecvp error\b',
    )

openmp_errors = (
    'OMP: Error',
    'libgomp: Thread creation failed',
    'libgomp: Out of memory',
    'libiomp5: error',
    )

openmp_error_patterns = (
    r'\bOMP:\s*Error\b',
    r'\blibgomp:\s*(?:Thread creation failed|Out of memory)\b',
    r'\blibiomp5:\s*error\b',
    )


# Compiled-code and language-runtime errors

linking_errors = (
    'error while loading shared libraries',
    # These fragments can be emitted while probing an optional plugin.
    # 'cannot open shared object file',
    # 'undefined symbol:',
    # 'wrong ELF class',
    'symbol lookup error',
    'relocation error',
    'GLIBCXX_3.4.30 not found',
    'CXXABI_1.3.13 not found',
    "version 'GLIBC_2.34' not found",
    )

linking_error_patterns = (
    r'\b(?:error while loading shared libraries|symbol lookup error|relocation error)',
    r'\b(?:GLIBCXX|CXXABI)_[0-9.]+\b[^\n]*\bnot found\b',
    r'\bversion\s+[\'`][^\'`]+[\'`]\s+not found\b',
    )

fortran_runtime_errors = (
    'Fortran runtime error:',
    'ERROR STOP',
    'forrtl: severe (174):',
    'Coarray ERROR STOP',
    # A coarray status value may be inspected and handled by the application.
    # 'Stat_Stopped_Image',
    )

fortran_error_patterns = (
    r'\bFortran runtime error:',
    r'\bERROR STOP\b',
    r'\bforrtl:\s*severe\s*\(\d+\):',
    r'\bCoarray\s+ERROR STOP\b',
    )

cpp_errors = (
    'terminate called after throwing an instance of',
    'terminating with uncaught exception of type',
    'terminate called without an active exception',
    'terminate called recursively',
    'Assertion failed',
    'double free or corruption',
    'corrupted size vs. prev_size',
    'free(): invalid pointer',
    'free(): double free detected',
    'malloc(): memory corruption',
    'munmap_chunk(): invalid pointer',
    'pure virtual method called',
    'AddressSanitizer:DEADLYSIGNAL',
    'ERROR: AddressSanitizer',
    'SUMMARY: AddressSanitizer',
    # These sanitizers can report issues without terminating the calculation.
    # 'ERROR: LeakSanitizer',
    # 'WARNING: ThreadSanitizer',
    # 'UndefinedBehaviorSanitizer',
    )

cpp_error_patterns = (
    r'\bterminate called after throwing an instance of\b',
    r'\bterminating with uncaught exception of type\b',
    r'\bterminate called (?:without an active exception|recursively)\b',
    r'\bAssertion failed\b',
    r'(?:\bdouble free or corruption\b|\bcorrupted size vs\. prev_size\b|free\(\): (?:invalid pointer|double free detected)|malloc\(\): memory corruption|munmap_chunk\(\): invalid pointer|\bpure virtual method called\b)',
    r'\b(?:AddressSanitizer:DEADLYSIGNAL|ERROR: AddressSanitizer|SUMMARY: AddressSanitizer)\b',
    # A what() line alone does not establish that an exception was uncaught.
    # r'^\s*what\(\):\s+.+$',
    )

cuda_errors = (
    'CUDA error:',
    # CUDA/NCCL return values can be checked and handled, and connection
    # failures can trigger transport fallback.  Require an application's
    # explicit "CUDA error:" diagnostic rather than a bare status token.
    # 'cudaErrorMemoryAllocation',
    # 'cudaErrorInitializationError',
    # 'cudaErrorLaunchFailure',
    # 'cudaErrorLaunchTimeout',
    # 'cudaErrorLaunchOutOfResources',
    # 'cudaErrorIllegalAddress',
    # 'cudaErrorNoKernelImageForDevice',
    # 'cudaErrorInsufficientDriver',
    # 'cudaErrorSystemDriverMismatch',
    # 'cudaErrorECCUncorrectable',
    # 'cudaErrorUnknown',
    # 'ncclUnhandledCudaError',
    # 'ncclSystemError',
    # 'ncclInternalError',
    # 'ncclInvalidArgument',
    # 'ncclInvalidUsage',
    # 'ncclRemoteError',
    # 'NCCL call to connect failed',
    # 'UCX call to connect failed',
    # 'CUDA call to connect failed',
    # 'socket call to connect failed',
    # 'transport call to connect failed',
    # An Xid is a driver event, not proof that this process failed.  NCCL WARN
    # includes warnings as well as errors; concrete NCCL failures are matched
    # below.
    # 'NVRM: Xid',
    # 'NCCL WARN',
    )

cuda_error_patterns = (
    r'\bCUDA error:',
    )

hip_errors = (
    'HIP error:',
    # HIP return values can be handled by a fallback path.
    # 'hipErrorMemoryAllocation',
    # 'hipErrorInitializationError',
    # 'hipErrorLaunchFailure',
    # 'hipErrorLaunchTimeOut',
    # 'hipErrorLaunchOutOfResources',
    # 'hipErrorIllegalAddress',
    # 'hipErrorNoBinaryForGpu',
    # 'hipErrorInsufficientDriver',
    # 'hipErrorECCNotCorrectable',
    # 'hipErrorUnknown',
    # This also occurs in status labels such as "ECC Error Count: 0".
    # 'ECC Error',
    # Kernel GPU messages can concern a different process on the node.
    # 'amdgpu: Page Fault',
    # 'amdgpu GPU fault',
    # 'amdgpu ring timeout',
    # 'amdgpu GPU reset',
    # 'amdgpu uncorrectable error',
    # 'kfd GPU fault',
    # 'kfd page fault',
    # 'kfd ring timeout',
    # 'kfd GPU reset',
    # 'kfd uncorrectable error',
    )

hip_error_patterns = (
    r'\bHIP error:',
    )


# Python-runtime errors.  Exception names are documented separately for
# examples, while matching requires traceback/final-exception structure.
python_exception_names = (
    'IndexError',
    'KeyError',
    'ValueError',
    'TypeError',
    'AttributeError',
    'NameError',
    'RuntimeError',
    'FileNotFoundError',
    'PermissionError',
    'ModuleNotFoundError',
    'ImportError',
    'ZeroDivisionError',
    'RecursionError',
    'SyntaxError',
    'IndentationError',
    'MemoryError',
    'AssertionError',
    'OSError',
    'EOFError',
    'OverflowError',
    'FloatingPointError',
    'TimeoutError',
    )

python_errors = (
    'Traceback (most recent call last):',
    'ExceptionGroup Traceback',
    'Fatal Python error',
    )

python_error_patterns = (
    r'^\s*Traceback \(most recent call last\):',
    r'^\s*ExceptionGroup Traceback\b',
    r'^\s*Fatal Python error\b',
    # Exception lines can be printed by handlers; traceback/fatal markers above
    # provide termination context.
    # r'^\s*(?:[\w.]+\.)?[A-Za-z_]\w*(?:Error|Exception)\s*:\s*.*$',
    )


# Compiled scientific-library errors

blas_errors = (
    'Intel MKL FATAL ERROR:',
    # BLAS argument errors and non-fatal library errors return to the caller,
    # which can select a fallback or otherwise handle them.
    # 'Intel MKL ERROR:',
    # 'OpenBLAS Error:',
    # 'on entry to DGEMM parameter number 1 had an illegal value',
    )

blas_error_patterns = (
    r'\bIntel MKL FATAL ERROR:',
    )

lapack_errors = (
    # LAPACK reports status to its caller; these conditions can be handled by
    # fallback algorithms and do not establish failure of the simulation.
    # 'LAPACK error:',
    # 'LAPACK native error:',
    # 'LAPACK computational failure:',
    # These numerical conditions can be handled by fallback algorithms.
    # 'matrix is exactly singular',
    # 'matrix is singular',
    # 'is not positive definite',
    # 'decomposition constraint violation',
    )

lapack_error_patterns = (
    # Numerical solver conditions can be handled by a fallback algorithm.
    # r'\b(?:LAPACK|[sdcz][a-z0-9_]{3,})[^\n]{0,100}\b(?:matrix is singular|is not positive definite|failed to converge|computational failure)\b',
    )

# Failure to import FFTW wisdom is recoverable and fftw_execute is merely an
# API name.  Keep the group present for future confirmed fatal diagnostics.
fftw_errors = ()
fftw_error_patterns = ()

hdf5_errors = (
    # Applications routinely probe optional files and objects and recover from
    # the resulting HDF5 error stack.
    # 'unable to open file',
    # 'unable to create file',
    # 'unable to open group',
    # 'unable to open dataset',
    # A failed write or invalid selection can concern optional checkpoint or
    # metadata output, and HDF5 returns these errors to the caller.
    # 'parallel write failed',
    # This is only an HDF5 error-stack classification, not a terminal outcome.
    # 'major: Parallel HDF5',
    # 'data space selection exceeds dataset dimensions',
    )

hdf5_error_patterns = (
    # HDF5 prints an error stack for failed probes even when the caller recovers.
    # r'HDF5-DIAG:\s*Error\s*detected',
    # r'\b(?:major|minor):\s*(?:file accessibility|unable to open file|unable to create file|write failed|read failed|object not found|bad value)\b',
    )

libxml2_errors = (
    # Parsing and validation failures can concern optional XML content and are
    # returned to the caller; termination must be established elsewhere.
    # 'parser error :',
    # 'This element is not expected',
    # 'Schemas validity error',
    # External entities can be optional and failure to load them is recoverable.
    # 'I/O error : Permission denied to access system file',
    # 'failed to load external entity',
    # 'Opening and ending tag mismatch',
    # 'Premature end of data',
    # 'XML validation failed',
    # 'XML element is not expected',
    )

libxml2_error_patterns = (
    # See libxml2_errors: these require separate termination context.
    )


# Python-module errors

numpy_errors = (
    # NumPy exceptions can be caught and handled by the calling application.
    # 'LinAlgError: calculation failed',
    # 'AxisError: calculation failed',
    # 'DTypePromotionError: calculation failed',
    # 'TooHardError: calculation failed',
    # '_ArrayMemoryError: calculation failed',
    )

numpy_error_patterns = (
    # A surrounding uncaught traceback must establish run failure.
    )

scipy_errors = (
    # SciPy exceptions and partial-convergence results can be handled.
    # 'ArpackError: calculation failed',
    # 'ArpackNoConvergence: calculation failed',
    # 'NoConvergence: calculation failed',
    # 'QhullError: calculation failed',
    # These are exception messages that callers can catch and recover from.
    # 'ARPACK error',
    # 'ARPACK iteration did not converge',
    # 'SuperLU factorization failed',
    # 'Factor is exactly singular',
    )

scipy_error_patterns = (
    # A surrounding uncaught traceback must establish run failure.
    )

h5py_errors = (
    # h5py exceptions are routinely caught during optional file/object probes.
    # 'CheckWriteEligibilityError: write is not permitted',
    # 'OSError: unable to open file',
    # 'RuntimeError: unable to create file',
    # 'ValueError: unable to read dataset',
    # 'OSError: unable to write dataset',
    # 'OSError: file signature not found',
    # "RuntimeError: object doesn't exist",
    # 'OSError: bad object header',
    # 'ValueError: address overflow',
    # 'OSError: no write intent',
    # These fragments can come from caught exceptions during optional probes.
    # 'file signature not found',
    # "object doesn't exist",
    # 'bad object header',
    # 'address overflow',
    # 'no write intent',
    )

h5py_error_patterns = (
    # A surrounding uncaught traceback must establish run failure.
    )


# Simulation-code errors

pwscf_errors = (
    'Error in routine',
    'Error in routine cdiaghg (1):',
    'bfgs failed',
    'bfgs failed: convergence not achieved',
    'convergence NOT achieved',
    'convergence NOT achieved after 100 iterations',
    'problems computing cholesky',
    'too many bands are not converged',
    )

pwscf_error_patterns = (
    r'\bError in routine\b',
    r'\bbfgs failed\b',
    r'\bconvergence\s+NOT\s+achieved\b',
    r'\bproblems computing cholesky\b',
    r'\btoo many bands are not converged\b',
    r'\bError in routine\s+[a-z0-9_]+\s*\(\d+\):',
    r'\bconvergence\s+NOT\s+achieved\s+after\s+\d+\s+iterations\b',
    r'\bbfgs failed\b[^\n]*\bconvergence not achieved\b',
    )

pyscf_errors = (
    'LibxcError: functional is not available',
    'SCF not converged',
    'CASSCF not converged',
    'UCASSCF not converged',
    'CCSD not converged',
    'Newton not converged',
    )

pyscf_error_patterns = (
    r'\b(?:SCF|CASSCF|UCASSCF|CCSD|Newton)[^\n]*\bnot converged\b',
    r'^\s*(?:pyscf[\w.]*\.)?LibxcError\s*:',
    )

quantum_package_errors = (
    'EZFIO error:',
    'FATAL ERROR:',
    'irp_error',
    'IRP_FATAL',
    'qp run: Error',
    'Too many determinants',
    'Selection failed',
    'Davidson not converged',
    'CIPSI not converged',
    'SCF not converged',
    'selection not converged',
    )

quantum_package_error_patterns = (
    r'(?:\bEZFIO error:|\bFATAL ERROR:|\birp_error\b|\bIRP_FATAL\b|\bqp run:\s*Error\b|\bToo many determinants\b|\bSelection failed\b)',
    r'\b(?:Davidson|CIPSI|SCF|selection)[^\n]*\bnot converged\b',
    )

rmg_errors = (
    'FATAL ERROR:',
    'CRITICAL:',
    'RMG Error:',
    'RMG Fatal:',
    'RMG Critical:',
    'RMGDFT Error:',
    'RMGDFT Fatal:',
    'RMGDFT Critical:',
    'Fatal RMG error',
    'Critical RMG error',
    'Fatal RMGDFT error',
    'Critical RMGDFT error',
    'SCF failed to converge',
    'SCF not converged',
    'multigrid failed',
    'Davidson breakdown',
    'subspace not converged',
    'domain decomposition failed',
    'grid decomposition failed',
    )

rmg_error_patterns = (
    r'^\s*(?:FATAL ERROR|CRITICAL):',
    r'\bRMG(?:DFT)?\s*(?:Error|Fatal|Critical)\s*:',
    r'\b(?:Fatal|Critical)\s+RMG(?:DFT)?\s+error\b',
    r'\bSCF[^\n]*\b(?:failed to converge|not converged)\b',
    r'\b(?:multigrid|Davidson|subspace)[^\n]*\b(?:failed|breakdown|not converged)\b',
    r'\b(?:domain decomposition|grid decomposition)[^\n]*\bfailed\b',
    )

qmcpack_errors = (
    'APP_ABORT',
    'Fatal Error',
    'Aborting at',
    'inconsistent input settings',
    'UniformCommunicateError',
    'barrier_and_abort',
    'Communicate::abort',
    )

qmcpack_error_patterns = (
    r'\bAPP_ABORT\b',
    r'\bUniformCommunicateError\b',
    r'\b(?:barrier_and_abort|Communicate::abort)\b',
    r'\bFatal Error\b',
    r'\bAborting at\b',
    r'\binconsistent input settings\b',
    )

vasp_errors = (
    'VERY BAD NEWS! internal error in subroutine',
    'ZBRENT: fatal error in bracketing',
    # VASP can continue from this warning and subsequently converge.
    # 'BRMIX: very serious problems',
    'EDDDAV: Call to ZHEGV failed',
    'EDDDAV: Call to ZHEEV failed',
    'EDDRMM: Call to ZHEGV failed',
    'EDDRMM: Call to ZHEEV failed',
    'LAPACK: Routine ZPOTRF failed',
    'ERROR FEXCP:',
    'ERROR: the triple product of the basis vectors',
    'ERROR: there must be 1 or 3 items on line 2 of POSCAR',
    )

vasp_error_patterns = (
    r'^\s*(?:\|\s*)?(?:VERY BAD NEWS!\s*)?(?:internal\s+)?error in subroutine\b',
    r'^\s*ZBRENT:\s*fatal\s+(?:error|internal)[^\n]*\bbracket',
    # BRMIX can be transient and followed by a converged, valid calculation.
    # r'^\s*BRMIX:\s*very serious problems\b',
    r'^\s*(?:EDDDAV|EDDRMM):[^\n]*(?:ZHEGV|ZHEEV)[^\n]*failed\b',
    r'^\s*LAPACK:[^\n]*\bfailed\b',
    r'^\s*ERROR FEXCP:',
    r'^\s*ERROR:\s*the triple product of the basis vectors\b',
    r'^\s*ERROR:\s*there must be 1 or 3 items on line 2 of POSCAR\b',
    )

gamess_errors = (
    'EXECUTION OF GAMESS TERMINATED -ABNORMALLY-',
    'SCF IS UNCONVERGED, TOO MANY ITERATIONS',
    'SCF DID NOT CONVERGE',
    'MEMORY REQUEST EXCEEDS AVAILABLE MEMORY',
    'WORDS OF MEMORY UNAVAILABLE',
    '1024 WORDS OF MEMORY UNAVAILABLE',
    'INPUT HAS AT LEAST ONE SPELLING OR LOGIC MISTAKE',
    'THIS JOB CANNOT CONTINUE',
    'ddikick.x: Fatal error detected',
    'ddikick.x: application process quit unexpectedly',
    'ddikick.x: application process 0 quit unexpectedly',
    'ddikick.x: Execution terminated due to error(s)',
    'DDI Process 0: error code 1',
    '*** ERROR TERMINATION ***',
    )

gamess_error_patterns = (
    r'\bEXECUTION OF GAMESS TERMINATED\s+-?ABNORMALLY-?(?!\w)',
    r'\bddikick\.x:\s*application process(?:\s+\d+)?\s+quit unexpectedly\b',
    r'\bDDI Process\s+\d+:\s*error code\s+(?!0\b)\d+\b',
    r'\bSCF\s+(?:IS UNCONVERGED,\s+TOO MANY ITERATIONS|DID NOT CONVERGE)\b',
    r'\bMEMORY REQUEST EXCEEDS AVAILABLE MEMORY\b',
    r'\b(?:\d+\s+)?WORDS OF MEMORY UNAVAILABLE\b',
    r'\bINPUT HAS AT LEAST ONE SPELLING OR LOGIC MISTAKE\b',
    r'\bTHIS JOB CANNOT CONTINUE\b',
    r'\bddikick\.x:\s*Fatal error detected\b',
    r'\bddikick\.x:\s*Execution terminated due to error\(s\)\.?',
    r'\*{3}\s*ERROR TERMINATION\s*\*{3}',
    )


# Human-readable keys and precise regex patterns are kept separately so the
# former can be used to document and test what the latter are intended to find.
_error_keys = {
    'shell'           : (),
    'linux_signals'   : (),
    'posix'           : (),
    'infiniband'      : infiniband_errors,
    'lustre'          : lustre_errors,
    'gpfs'            : gpfs_errors,
    'slurm'           : slurm_errors,
    'pbs'             : pbs_errors,
    'mpi'             : mpi_errors,
    'openmp'          : openmp_errors,
    'linking'         : linking_errors,
    'fortran'         : fortran_runtime_errors,
    'cpp'             : cpp_errors,
    'cuda'            : cuda_errors,
    'hip'             : hip_errors,
    'python'          : python_errors,
    'blas'            : blas_errors,
    'lapack'          : lapack_errors,
    'fftw'            : fftw_errors,
    'hdf5'            : hdf5_errors,
    'libxml2'         : libxml2_errors,
    'numpy'           : (),
    # Exception class names from scipy_errors are matched only as structured
    # exception lines by scipy_error_patterns.  Message fragments alone do not
    # establish that the exception was left unhandled.
    'scipy'           : (
        # 'ARPACK error',
        # 'ARPACK iteration did not converge',
        # 'SuperLU factorization failed',
        # 'Factor is exactly singular',
        ),
    # CheckWriteEligibilityError is anchored as an exception line by
    # h5py_error_patterns.  HDF5 message fragments can result from caught
    # exceptions during optional probes.
    'h5py'            : (
        # 'file signature not found',
        # "object doesn't exist",
        # 'bad object header',
        # 'address overflow',
        # 'no write intent',
        ),
    'pwscf'           : pwscf_errors,
    # LibxcError is matched as a structured exception line by
    # pyscf_error_patterns.
    'pyscf'           : ('SCF not converged',),
    'quantum_package' : quantum_package_errors,
    'rmg'             : rmg_errors,
    'qmcpack'         : qmcpack_errors,
    'vasp'            : vasp_errors,
    'gamess'          : gamess_errors,
    }

_error_patterns = {
    'shell'           : shell_error_patterns,
    'linux_signals'   : linux_signal_error_patterns,
    'posix'           : posix_error_patterns,
    'infiniband'      : infiniband_error_patterns,
    'lustre'          : lustre_error_patterns,
    'gpfs'            : gpfs_error_patterns,
    'slurm'           : slurm_error_patterns,
    'pbs'             : pbs_error_patterns,
    'mpi'             : mpi_error_patterns,
    'openmp'          : openmp_error_patterns,
    'linking'         : linking_error_patterns,
    'fortran'         : fortran_error_patterns,
    'cpp'             : cpp_error_patterns,
    'cuda'            : cuda_error_patterns,
    'hip'             : hip_error_patterns,
    'python'          : python_error_patterns,
    'blas'            : blas_error_patterns,
    'lapack'          : lapack_error_patterns,
    'fftw'            : fftw_error_patterns,
    'hdf5'            : hdf5_error_patterns,
    'libxml2'         : libxml2_error_patterns,
    'numpy'           : numpy_error_patterns,
    'scipy'           : scipy_error_patterns,
    'h5py'            : h5py_error_patterns,
    'pwscf'           : pwscf_error_patterns,
    'pyscf'           : pyscf_error_patterns,
    'quantum_package' : quantum_package_error_patterns,
    'rmg'             : rmg_error_patterns,
    'qmcpack'         : qmcpack_error_patterns,
    'vasp'            : vasp_error_patterns,
    'gamess'          : gamess_error_patterns,
    }

_error_set_names = tuple(_error_keys)


def _literal_error_pattern(error_key):
    """Escape a readable key while allowing flexible whitespace."""
    pattern = r'\s+'.join(re.escape(part) for part in error_key.split())
    if error_key and (error_key[0].isalnum() or error_key[0] == '_'):
        pattern = r'(?<!\w)' + pattern
    if error_key and (error_key[-1].isalnum() or error_key[-1] == '_'):
        pattern += r'(?!\w)'
    return pattern


@cache
def _combined_error_pattern(enabled_sets):
    patterns = []
    seen = set()
    for set_name in enabled_sets:
        set_patterns = [
            _literal_error_pattern(key) for key in _error_keys[set_name]
            ]
        set_patterns.extend(_error_patterns[set_name])
        for pattern in set_patterns:
            if pattern not in seen:
                seen.add(pattern)
                patterns.append(pattern)
    if not patterns:
        return None
    patterns.sort(key=len, reverse=True)
    return re.compile(
        r'(?:{})'.format('|'.join(patterns)),
        re.IGNORECASE | re.MULTILINE,
        )


def _read_error_text(source):
    if hasattr(source, 'read'):
        text = source.read()
    elif isinstance(source, os.PathLike):
        with open(source, 'r', errors='replace') as infile:
            text = infile.read()
    elif isinstance(source, str):
        try:
            is_file = '\n' not in source and os.path.isfile(source)
        except OSError:
            is_file = False
        if is_file:
            with open(source, 'r', errors='replace') as infile:
                text = infile.read()
        else:
            text = source
    else:
        raise TypeError(
            'source must be text, a path-like object, or an open text file'
            )
    if not isinstance(text, str):
        raise TypeError('source must provide text rather than binary data')
    return text


def find_error_keys(
        source,
        # select error batches
        *,
        all_errors         = False,
        operating_system   = False,
        hpc                = False,
        code               = False,
        code_library       = False,
        python_module      = False,
        # operating system errors
        shell              = False,
        linux_signals      = False,
        posix              = False,
        # hpc errors
        infiniband         = False,
        lustre             = False,
        gpfs               = False,
        slurm              = False,
        pbs                = False,
        mpi                = False,
        openmp             = False,
        # code errors
        linking            = False,
        fortran            = False,
        cpp                = False,
        cuda               = False,
        hip                = False,
        # python code errors
        python             = False,
        # code library errors
        blas               = False,
        lapack             = False,
        fftw               = False,
        hdf5               = False,
        libxml2            = False,
        # python module errors
        numpy              = False,
        scipy              = False,
        h5py               = False,
        # simulation code errors
        pwscf              = False,
        pyscf              = False,
        quantum_package    = False,
        rmg                = False,
        qmcpack            = False,
        vasp               = False,
        gamess             = False,
        # return lines found or not
        return_lines       = False,
        ):
    """Find likely failure diagnostics in scientific-application output.

    Parameters
    ----------
    source : str, os.PathLike, or text file
        Text to search, a path to a text file, or an open text stream.  A
        string naming an existing file is interpreted as a path; all other
        strings are interpreted as text.  Callers inspecting a simulation
        must search its standard output and standard error separately.
    all_errors : bool, optional
        Enable every individual error set.
    operating_system : bool, optional
        Enable ``shell``, ``linux_signals``, and ``posix``.
    hpc : bool, optional
        Enable ``infiniband``, ``lustre``, ``gpfs``, ``slurm``, ``pbs``,
        ``mpi``, and ``openmp``.
    code : bool, optional
        Enable ``linking``, ``fortran``, ``cpp``, ``cuda``, and ``hip``.
    code_library : bool, optional
        Enable ``blas``, ``lapack``, ``fftw``, ``hdf5``, and ``libxml2``.
    python_module : bool, optional
        Enable ``numpy``, ``scipy``, and ``h5py``.  These sets intentionally
        exclude standalone exception lines; also enable ``python`` to detect
        uncaught failures through their traceback.
    shell, linux_signals, posix : bool, optional
        Select individual operating-system error sets.
    infiniband, lustre, gpfs, slurm, pbs, mpi, openmp : bool, optional
        Select individual HPC environment error sets.
    linking, fortran, cpp, cuda, hip : bool, optional
        Select individual compiled-code and runtime error sets.
    python : bool, optional
        Select uncaught Python and interpreter errors.
    blas, lapack, fftw, hdf5, libxml2 : bool, optional
        Select individual compiled-library error sets.
    numpy : bool, optional
        Select NumPy-specific errors.  A NumPy exception can be handled, so
        no standalone exception is currently definitive; enable ``python``
        to detect an uncaught exception through its traceback.
    scipy : bool, optional
        Select SciPy-specific errors.  A SciPy exception can be handled, so
        no standalone exception is currently definitive; enable ``python``
        to detect an uncaught exception through its traceback.
    h5py : bool, optional
        Select h5py-specific errors.  An h5py exception can be handled, so no
        standalone exception is currently definitive; enable ``python`` to
        detect an uncaught exception through its traceback.
    pwscf, quantum_package, rmg, qmcpack, vasp, gamess : bool, optional
        Select individual simulation-code error sets.
    pyscf : bool, optional
        Select PySCF errors and uncaught Python/interpreter errors.  Enabling
        this selector also enables ``python``.
    return_lines : bool, optional
        If ``True``, return the lines containing matching diagnostics in
        addition to the status.

    Returns
    -------
    found : bool
        Whether at least one enabled failure expression was found.
    lines : list of str
        Matching lines, returned only when ``return_lines=True``.  Repeated
        matching lines are retained.

    Notes
    -----
    Selectors are additive and all default to ``False``.  If no set is
    selected, a :class:`ValueError` is raised.

    Enabled keys and regexes are combined into one cached, case-insensitive
    expression.  The patterns favor diagnostics that make failure to produce
    intended simulation output likely; completion and output validity are
    expected to be assessed separately.

    Examples
    --------
    >>> find_error_keys("run complete", qmcpack=True)
    False
    >>> find_error_keys(
    ...     "QMCPACK fatal path: APP_ABORT invalid input",
    ...     qmcpack=True,
    ...     )
    True
    >>> output = '''
    ... step 1
    ... LAPACK computational failure: ZHEEV did not converge
    ... '''
    >>> find_error_keys(output, lapack=True, return_lines=True)
    (True, ['LAPACK computational failure: ZHEEV did not converge'])
    """
    flags = {
        'shell'           : shell,
        'linux_signals'   : linux_signals,
        'posix'           : posix,
        'infiniband'      : infiniband,
        'lustre'          : lustre,
        'gpfs'            : gpfs,
        'slurm'           : slurm,
        'pbs'             : pbs,
        'mpi'             : mpi,
        'openmp'          : openmp,
        'linking'         : linking,
        'fortran'         : fortran,
        'cpp'             : cpp,
        'cuda'            : cuda,
        'hip'             : hip,
        'python'          : python,
        'blas'            : blas,
        'lapack'          : lapack,
        'fftw'            : fftw,
        'hdf5'            : hdf5,
        'libxml2'         : libxml2,
        'numpy'           : numpy,
        'scipy'           : scipy,
        'h5py'            : h5py,
        'pwscf'           : pwscf,
        'pyscf'           : pyscf,
        'quantum_package' : quantum_package,
        'rmg'             : rmg,
        'qmcpack'         : qmcpack,
        'vasp'            : vasp,
        'gamess'          : gamess,
        }

    if all_errors:
        for name in flags:
            flags[name] = True
    else:
        if operating_system:
            for name in ('shell', 'linux_signals', 'posix'):
                flags[name] = True
        if hpc:
            for name in (
                    'infiniband', 'lustre', 'gpfs', 'slurm', 'pbs', 'mpi',
                    'openmp',
                    ):
                flags[name] = True
        if code:
            for name in ('linking', 'fortran', 'cpp', 'cuda', 'hip'):
                flags[name] = True
        if code_library:
            for name in ('blas', 'lapack', 'fftw', 'hdf5', 'libxml2'):
                flags[name] = True
        if python_module:
            for name in ('numpy', 'scipy', 'h5py'):
                flags[name] = True

    # PySCF is a Python simulation code.  A traceback or fatal interpreter
    # diagnostic is therefore a PySCF run failure even when its exception type
    # is not one of the PySCF-specific diagnostics above.
    if flags['pyscf']:
        flags['python'] = True

    enabled_sets = tuple(name for name in _error_set_names if flags[name])
    if len(enabled_sets)==0:
        raise ValueError(
            'at least one error set must be requested by find_error_keys'
            )

    pattern = _combined_error_pattern(enabled_sets)
    if pattern is None:
        if not return_lines:
            return False
        else:
            return False,list()

    text = _read_error_text(source)
    if not return_lines:
        return pattern.search(text) is not None

    lines = [line for line in text.splitlines() if pattern.search(line) is not None]
    errors_found = len(lines)>0
    return errors_found,lines
#end def find_error_keys
