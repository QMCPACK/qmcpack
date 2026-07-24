"""Error diagnostics commonly printed by scientific applications.

The tuples ending in ``_errors`` contain readable examples of diagnostics.
The tuples ending in ``_error_patterns`` contain regular expressions used to
match variable portions or to add context that reduces false positives.
"""

import os
import re
from functools import lru_cache


# Operating-system errors

shell_errors = (
    'Segmentation fault',
    'Floating point exception',
    'Illegal instruction',
    'Bus error',
    'Bad system call',
    'Stack overflow',
    'Out of memory',
    'invoked oom-killer',
    'Machine Check Exception',
    'EDAC Hardware Error',
    'stack smashing detected',
    'general protection fault',
    )

shell_error_patterns = (
    r'^.*\b(?:segmentation fault|floating point exception|'
    r'illegal instruction|bus error|bad system call)'
    r'(?:\s+\(core dumped\))?\s*$',
    r'^.*\b(?:aborted|killed)(?:\s+\(core dumped\))?\s*$',
    r'^.*\bterminated\s*$',
    r'\b(?:out of memory|cannot allocate memory|oom-kill(?:er)?|'
    r'invoked oom-killer|killed process\s+\d+)\b',
    r'\b(?:stack smashing detected|general protection fault|'
    r'Machine Check Exception)\b',
    r'\b(?:MCE|EDAC)[^\n]{0,80}\bHardware Error\b',
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
    )

linux_signal_error_patterns = (
    r'\b(?:terminated|killed|exited|aborted|died|received signal)\b'
    r'[^\n]{0,80}\bSIG(?:HUP|ILL|ABRT|FPE|KILL|SEGV|PIPE|TERM|BUS|SYS|'
    r'TRAP|XCPU|XFSZ)\b',
    r'\bSIG(?:HUP|ILL|ABRT|FPE|KILL|SEGV|PIPE|TERM|BUS|SYS|TRAP|XCPU|'
    r'XFSZ)\b[^\n]{0,80}\b(?:terminated|killed|exited|aborted|died)\b',
    r'\bterminated with signal\s+\d+\b',
    r'\bexited on signal\s+\d+\b',
    )

posix_errors = (
    'No such file or directory',
    'Permission denied',
    'Not a directory',
    'Is a directory',
    'No space left on device',
    'Too many open files',
    'Cannot allocate memory',
    'Connection refused',
    'Connection timed out',
    'Network is unreachable',
    'Address already in use',
    'Broken pipe',
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
    r'^(?:.*:\s*)?(?:no such file or directory|permission denied|'
    r'not a directory|is a directory|no space left on device|'
    r'too many open files|cannot allocate memory|connection refused|'
    r'connection timed out|network is unreachable|address already in use|'
    r'broken pipe)\s*$',
    r'\b(?:fatal|error|exception|failed|cannot|unable)[^\n]{0,100}'
    r'\b(?:no such file or directory|permission denied|not a directory|'
    r'is a directory|no space left on device|too many open files|'
    r'cannot allocate memory|connection refused|connection timed out|'
    r'network is unreachable|address already in use|broken pipe)\b',
    r'\b(?:errno|error|failed|failure|fatal)[^\n]{0,40}'
    r'\b(?:ENOENT|EACCES|EISDIR|ENOTDIR|ENOSPC|EMFILE|ENOMEM|'
    r'ECONNREFUSED|ETIMEDOUT|ENETUNREACH|EADDRINUSE|EPIPE|EIO|ENXIO|'
    r'EBADF|EBUSY|ENODEV|EROFS|EDQUOT|ECONNRESET|EHOSTUNREACH|ENOTCONN)\b',
    )


# HPC environment errors

infiniband_errors = (
    'UCX ERROR',
    'libibverbs: malformed packet',
    )

infiniband_error_patterns = (
    r'\bUCX\s+ERROR\b',
    r'\b(?:ucp|uct|ucs)_[a-z0-9_]+\b[^\n]{0,80}'
    r'\b(?:failed|error|unreachable|timed out)\b',
    r'\bibv_[a-z0-9_]+\b[^\n]{0,80}\b(?:failed|error)\b',
    r'\b(?:libfabric|ofi_[a-z0-9_]+)\b[^\n]{0,80}'
    r'\b(?:error|failed|unreachable)\b',
    )

lustre_errors = (
    'LustreError:',
    'LBUG:',
    )

lustre_error_patterns = (
    r'\bLNet(?:Error)?\b[^\n]{0,80}\b(?:error|failed|fatal|timeout|'
    r'unreachable)\b',
    r'\bLustre\b[^\n]{0,80}\b(?:error|failed|fatal|evicted)\b',
    )

gpfs_errors = (
    'GPFS: [ERROR]',
    'GPFS: [FATAL]',
    )

gpfs_error_patterns = (
    r'\bGPFS\b[^\n]{0,80}\b(?:deadlock detected|disk unavailable|'
    r'unmounted abnormally|token expired)\b',
    )

slurm_errors = (
    'slurmstepd: error:',
    'srun: error:',
    'srun: Force term; sending SIGKILL',
    'DUE TO TIME LIMIT',
    'Exceeded job memory limit',
    )

slurm_error_patterns = (
    r'\bState=(?:FAILED|TIMEOUT|NODE_FAIL|OUT_OF_MEMORY|BOOT_FAIL|'
    r'DEADLINE|CANCELLED|PREEMPTED)\b',
    r'\b(?:JOB|STEP)[^\n]{0,80}\b(?:CANCELLED|FAILED|OUT_OF_MEMORY|'
    r'TIMEOUT|NODE_FAIL)\b',
    r'\b(?:launch failed|oom-kill)\b',
    )

pbs_errors = (
    'PBS: job killed:',
    'ob_init: Unable to read server database',
    'cannot send job to mom',
    'qsub: Bad UID for job execution',
    )

pbs_error_patterns = (
    r'\bexit_status\s*=\s*(?!0\b)-?\d+\b',
    )

mpi_errors = (
    'MPI_ABORT',
    'MPI_ERR_BUFFER',
    'MPI_ERR_COUNT',
    'MPI_ERR_TYPE',
    'MPI_ERR_TAG',
    'MPI_ERR_COMM',
    'MPI_ERR_RANK',
    'MPI_ERR_REQUEST',
    'MPI_ERR_ROOT',
    'MPI_ERR_GROUP',
    'MPI_ERR_OP',
    'MPI_ERR_TOPOLOGY',
    'MPI_ERR_DIMS',
    'MPI_ERR_ARG',
    'MPI_ERR_UNKNOWN',
    'MPI_ERR_TRUNCATE',
    'MPI_ERR_OTHER',
    'MPI_ERR_INTERN',
    'MPI_ERR_IN_STATUS',
    'MPI_ERR_PENDING',
    'MPI_ERR_ACCESS',
    'MPI_ERR_AMODE',
    'MPI_ERR_ASSERT',
    'MPI_ERR_BAD_FILE',
    'MPI_ERR_BASE',
    'MPI_ERR_CONVERSION',
    'MPI_ERR_DISP',
    'MPI_ERR_DUP_DATAREP',
    'MPI_ERR_FILE_EXISTS',
    'MPI_ERR_FILE_IN_USE',
    'MPI_ERR_FILE',
    'MPI_ERR_INFO_KEY',
    'MPI_ERR_INFO_NOKEY',
    'MPI_ERR_INFO_VALUE',
    'MPI_ERR_INFO',
    'MPI_ERR_IO',
    'MPI_ERR_KEYVAL',
    'MPI_ERR_LOCKTYPE',
    'MPI_ERR_NAME',
    'MPI_ERR_NO_MEM',
    'MPI_ERR_NOT_SAME',
    'MPI_ERR_NO_SPACE',
    'MPI_ERR_NO_SUCH_FILE',
    'MPI_ERR_PORT',
    'MPI_ERR_PROC_ABORTED',
    'MPI_ERR_QUOTA',
    'MPI_ERR_READ_ONLY',
    'MPI_ERR_RMA_CONFLICT',
    'MPI_ERR_RMA_SYNC',
    'MPI_ERR_SERVICE',
    'MPI_ERR_SIZE',
    'MPI_ERR_SPAWN',
    'MPI_ERR_UNSUPPORTED_DATAREP',
    'MPI_ERR_UNSUPPORTED_OPERATION',
    'MPI_ERR_WIN',
    'MPI_T_ERR_MEMORY',
    'MPI_T_ERR_NOT_INITIALIZED',
    'MPI_T_ERR_CANNOT_INIT',
    'MPI_T_ERR_INVALID_INDEX',
    'MPI_T_ERR_INVALID_ITEM',
    'MPI_T_ERR_INVALID_HANDLE',
    'MPI_T_ERR_OUT_OF_HANDLES',
    'MPI_T_ERR_OUT_OF_SESSIONS',
    'MPI_T_ERR_INVALID_SESSION',
    'MPI_T_ERR_CVAR_SET_NOT_NOW',
    'MPI_T_ERR_CVAR_SET_NEVER',
    'MPI_T_ERR_PVAR_NO_STARTSTOP',
    'MPI_T_ERR_PVAR_NO_WRITE',
    'MPI_T_ERR_PVAR_NO_ATOMIC',
    'MPI_ERR_RMA_RANGE',
    'MPI_ERR_RMA_ATTACH',
    'MPI_ERR_RMA_FLAVOR',
    'MPI_ERR_RMA_SHARED',
    'MPI_T_ERR_INVALID',
    'MPI_T_ERR_INVALID_NAME',
    'MPI_ERR_SESSION',
    'mpirun: kill job',
    'mpirun noticed that process rank',
    'ORTE_ERROR_LOG',
    'PRTE_ERROR_LOG',
    'the first job to fail is listed below',
    'job aborted:',
    'mpiexec_callback_proc',
    'cleaning up processes',
    'execvp error',
    )

mpi_error_patterns = (
    r'\bMPI_Abort was invoked\b',
    r'\bone or more processes exited with non-zero status\b',
    r'\bprocess returned a non-zero exit code\b',
    r'\bPrimary job terminated normally, but\b',
    r'\b(?:mpirun|mpiexec|orterun|prterun)\b[^\n]{0,120}'
    r'\b(?:aborted|failed|non-zero|signal|terminated)\b',
    r'\b(?:exited on|terminated with) signal(?:\s+\d+)?\b',
    )

openmp_errors = (
    'OMP: Error',
    'libgomp: Thread creation failed',
    'libgomp: Out of memory',
    'libiomp5: error',
    )

openmp_error_patterns = ()


# Compiled-code and language-runtime errors

linking_errors = (
    'error while loading shared libraries',
    'cannot open shared object file',
    'undefined symbol:',
    'wrong ELF class',
    'symbol lookup error',
    'relocation error',
    )

linking_error_patterns = (
    r'\b(?:GLIBCXX|CXXABI)_[0-9.]+\b[^\n]{0,40}\bnot found\b',
    r'\bversion\s+[\'`][^\'`]+[\'`]\s+not found\b',
    )

fortran_runtime_errors = (
    'Fortran runtime error:',
    'ERROR STOP',
    'Stat_Stopped_Image',
    )

fortran_error_patterns = (
    r'\bforrtl:\s*severe\s*\(\d+\):',
    r'\bCoarray\s+ERROR STOP\b',
    )

cplusplus_errors = (
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
    'ERROR: LeakSanitizer',
    'WARNING: ThreadSanitizer',
    'UndefinedBehaviorSanitizer',
    )

cplusplus_error_patterns = (
    r'^\s*what\(\):\s+.+$',
    )

cuda_errors = (
    'CUDA error:',
    'NVRM: Xid',
    'NCCL WARN',
    )

cuda_error_patterns = (
    r'\bcudaError(?:MemoryAllocation|InitializationError|LaunchFailure|'
    r'LaunchTimeout|LaunchOutOfResources|IllegalAddress|'
    r'NoKernelImageForDevice|InsufficientDriver|SystemDriverMismatch|'
    r'ECCUncorrectable|Unknown)\b',
    r'\bnccl(?:UnhandledCudaError|SystemError|InternalError|'
    r'InvalidArgument|InvalidUsage|RemoteError)\b',
    r'\b(?:NCCL|UCX|CUDA|socket|transport)[^\n]{0,80}'
    r'\bcall to connect failed\b',
    )

hip_errors = (
    'HIP error:',
    'ECC Error',
    'amdgpu: Page Fault',
    )

hip_error_patterns = (
    r'\bhipError(?:MemoryAllocation|InitializationError|LaunchFailure|'
    r'LaunchTimeOut|LaunchOutOfResources|IllegalAddress|NoBinaryForGpu|'
    r'InsufficientDriver|ECCNotCorrectable|Unknown)\b',
    r'\b(?:amdgpu|kfd)[^\n]{0,100}\b(?:GPU fault|page fault|ring timeout|'
    r'GPU reset|uncorrectable)\b',
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
    r'^\s*(?:[\w.]+\.)?[A-Za-z_]\w*(?:Error|Exception)\s*:\s*.*$',
    )


# Compiled scientific-library errors

blas_errors = (
    'Intel MKL ERROR:',
    'Intel MKL FATAL ERROR:',
    'OpenBLAS Error:',
    )

blas_error_patterns = (
    r'\bon entry to\s+[a-z0-9_]+\s+parameter(?: number)?\s+\d+'
    r'\s+had an illegal value\b',
    )

lapack_errors = (
    'LAPACK error:',
    'LAPACK native error:',
    'LAPACK computational failure:',
    'matrix is exactly singular',
    'matrix is singular',
    'is not positive definite',
    'decomposition constraint violation',
    )

lapack_error_patterns = (
    r'\b(?:LAPACK|[sdcz][a-z0-9_]{3,})[^\n]{0,100}'
    r'\b(?:matrix is singular|is not positive definite|failed to converge|'
    r'computational failure)\b',
    )

# Failure to import FFTW wisdom is recoverable and fftw_execute is merely an
# API name.  Keep the group present for future confirmed fatal diagnostics.
fftw_errors = ()
fftw_error_patterns = ()

hdf5_errors = (
    'unable to open file',
    'unable to create file',
    'unable to open group',
    'unable to open dataset',
    'parallel write failed',
    'major: Parallel HDF5',
    'data space selection exceeds dataset dimensions',
    )

hdf5_error_patterns = (
    r'HDF5-DIAG:\s*Error\s*detected',
    r'\b(?:major|minor):\s*(?:file accessibility|unable to open file|'
    r'unable to create file|write failed|read failed|object not found|'
    r'bad value)\b',
    )

libxml2_errors = (
    'parser error :',
    'This element is not expected',
    'Schemas validity error',
    'I/O error : Permission denied to access system file',
    'failed to load external entity',
    'Opening and ending tag mismatch',
    'Premature end of data',
    )

libxml2_error_patterns = (
    r'\b(?:parser|schemas?|xml)[^\n]{0,80}'
    r'\b(?:error|validation failed|not expected|failed to load)\b',
    )


# Python-module errors

numpy_errors = (
    'LinAlgError',
    'AxisError',
    'DTypePromotionError',
    'TooHardError',
    '_ArrayMemoryError',
    )

numpy_error_patterns = (
    r'^\s*(?:numpy[\w.]*\.)?(?:LinAlgError|AxisError|'
    r'DTypePromotionError|TooHardError|_ArrayMemoryError)\s*:',
    )

scipy_errors = (
    'ArpackError',
    'ArpackNoConvergence',
    'NoConvergence',
    'QhullError',
    'ARPACK error',
    'ARPACK iteration did not converge',
    'SuperLU factorization failed',
    'Factor is exactly singular',
    )

scipy_error_patterns = (
    r'^\s*(?:scipy[\w.]*\.)?(?:ArpackError|ArpackNoConvergence|'
    r'NoConvergence|QhullError)\s*:',
    )

h5py_errors = (
    'CheckWriteEligibilityError',
    'file signature not found',
    "object doesn't exist",
    'bad object header',
    'address overflow',
    'no write intent',
    )

h5py_error_patterns = (
    r'^\s*(?:h5py[\w.]*\.)?CheckWriteEligibilityError\s*:',
    r'^\s*(?:OSError|RuntimeError|ValueError):[^\n]*'
    r'(?:unable to (?:open|create|read|write)|file signature not found|'
    r'object doesn\'t exist|bad object header|address overflow|'
    r'no write intent)',
    )


# Simulation-code errors

pwscf_errors = (
    'Error in routine',
    'bfgs failed',
    'convergence NOT achieved',
    'problems computing cholesky',
    'too many bands are not converged',
    )

pwscf_error_patterns = (
    r'\bError in routine\s+[a-z0-9_]+\s*\(\d+\):',
    r'\bconvergence\s+NOT\s+achieved\s+after\s+\d+\s+iterations\b',
    r'\bbfgs failed\b[^\n]*\bconvergence not achieved\b',
    )

pyscf_errors = (
    'LibxcError',
    'SCF not converged',
    )

pyscf_error_patterns = (
    r'\b(?:SCF|CASSCF|UCASSCF|CCSD|Newton)[^\n]{0,40}'
    r'\bnot converged\b',
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
    )

quantum_package_error_patterns = (
    r'\b(?:Davidson|CIPSI|SCF|selection)[^\n]{0,60}'
    r'\bnot converged\b',
    )

rmg_errors = (
    'FATAL ERROR:',
    'CRITICAL:',
    'RMG Error:',
    )

rmg_error_patterns = (
    r'\bRMG(?:DFT)?\s*(?:Error|Fatal|Critical)\s*:',
    r'\b(?:Fatal|Critical)\s+RMG(?:DFT)?\s+error\b',
    r'\bSCF[^\n]{0,60}\b(?:failed to converge|not converged)\b',
    r'\b(?:multigrid|Davidson|subspace)[^\n]{0,60}'
    r'\b(?:failed|breakdown|not converged)\b',
    r'\b(?:domain decomposition|grid decomposition)[^\n]{0,60}'
    r'\bfailed\b',
    )

qmcpack_errors = (
    'APP_ABORT',
    'Fatal Error',
    'Aborting at',
    'inconsistent input settings',
    'UniformCommunicateError',
    'barrier_and_abort',
    )

qmcpack_error_patterns = (
    r'\bAPP_ABORT\b',
    r'\bUniformCommunicateError\b',
    r'\b(?:barrier_and_abort|Communicate::abort)\b',
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
    'cplusplus'       : cplusplus_errors,
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
    # exception lines by scipy_error_patterns.  These library messages are
    # sufficiently specific to search for literally.
    'scipy'           : (
        'ARPACK error',
        'ARPACK iteration did not converge',
        'SuperLU factorization failed',
        'Factor is exactly singular',
        ),
    # CheckWriteEligibilityError is anchored as an exception line by
    # h5py_error_patterns.  These HDF5-specific message fragments are retained
    # as safe literal matches.
    'h5py'            : (
        'file signature not found',
        "object doesn't exist",
        'bad object header',
        'address overflow',
        'no write intent',
        ),
    'pwscf'           : pwscf_errors,
    # LibxcError is matched as a structured exception line by
    # pyscf_error_patterns; PySCF's explicit convergence message is safe to
    # search for literally.
    'pyscf'           : ('SCF not converged',),
    'quantum_package' : quantum_package_errors,
    'rmg'             : rmg_errors,
    'qmcpack'         : qmcpack_errors,
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
    'cplusplus'       : cplusplus_error_patterns,
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


@lru_cache(maxsize=None)
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
        cplusplus          = False,
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
        # return lines found or not
        return_lines       = False,
        ):
    """Find likely failure diagnostics in scientific-application output.

    Parameters
    ----------
    source : str, os.PathLike, or text file
        Text to search, a path to a text file, or an open text stream.  A
        string naming an existing file is interpreted as a path; all other
        strings are interpreted as text.
    all_errors : bool, optional
        Enable every individual error set.
    operating_system : bool, optional
        Enable ``shell``, ``linux_signals``, and ``posix``.
    hpc : bool, optional
        Enable ``infiniband``, ``lustre``, ``gpfs``, ``slurm``, ``pbs``,
        ``mpi``, and ``openmp``.
    code : bool, optional
        Enable ``linking``, ``fortran``, ``cplusplus``, ``cuda``, and ``hip``.
    code_library : bool, optional
        Enable ``blas``, ``lapack``, ``fftw``, ``hdf5``, and ``libxml2``.
    python_module : bool, optional
        Enable ``numpy``, ``scipy``, and ``h5py``.
    shell, linux_signals, posix : bool, optional
        Select individual operating-system error sets.
    infiniband, lustre, gpfs, slurm, pbs, mpi, openmp : bool, optional
        Select individual HPC environment error sets.
    linking, fortran, cplusplus, cuda, hip : bool, optional
        Select individual compiled-code and runtime error sets.
    python : bool, optional
        Select uncaught Python and interpreter errors.
    blas, lapack, fftw, hdf5, libxml2 : bool, optional
        Select individual compiled-library error sets.
    numpy, scipy, h5py : bool, optional
        Select individual Python-module error sets.
    pwscf, pyscf, quantum_package, rmg, qmcpack : bool, optional
        Select individual simulation-code error sets.
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
    >>> find_error_keys(
    ...     "step 1\\nLAPACK computational failure: ZHEEV did not converge",
    ...     lapack=True,
    ...     return_lines=True,
    ...     )
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
        'cplusplus'       : cplusplus,
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
            for name in ('linking', 'fortran', 'cplusplus', 'cuda', 'hip'):
                flags[name] = True
        if code_library:
            for name in ('blas', 'lapack', 'fftw', 'hdf5', 'libxml2'):
                flags[name] = True
        if python_module:
            for name in ('numpy', 'scipy', 'h5py'):
                flags[name] = True

    enabled_sets = tuple(name for name in _error_set_names if flags[name])
    if not enabled_sets:
        raise ValueError(
            'at least one error set must be requested by find_error_keys'
            )

    pattern = _combined_error_pattern(enabled_sets)
    if pattern is None:
        return (False, []) if return_lines else False

    text = _read_error_text(source)
    if not return_lines:
        return pattern.search(text) is not None

    lines = [line for line in text.splitlines() if pattern.search(line)]
    return bool(lines), lines
#end def find_error_keys
