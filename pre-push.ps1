<#
    Windows equivalent of ./pre-push (the bash script) for boost-multi.

    Builds and tests the project with whatever compiler toolchains are actually
    installed on this machine, mirroring the spirit of the Linux/macOS script's
    multi-compiler matrix (several build dirs, warnings-as-errors, a sanitizer
    build) using the MSVC/LLVM tools that ship with Visual Studio.

    Detected on this machine: cl.exe (MSVC) + bundled Ninja, via
    "Visual Studio Build Tools 2022".
    ccache is optional; the script skips using it as a compiler launcher with
    a warning instead of failing when it is missing. To add it later:

      winget install ccache   # or: scoop install ccache / choco install ccache

    nvcc (CUDA Toolkit) is also optional: when found, a 4th variant builds
    with -DENABLE_CUDA=1 using cl.exe as the CUDA host compiler, exercising
    the .cu sources under include/boost/multi/adaptors/cuda/cublas/test.
    Skipped with a warning when nvcc isn't installed. To add it (pinned to
    the version this script is tested against):

      winget install --id Nvidia.CUDA --version 12.9

    clang-cl (LLVM's MSVC-compatible driver) is also optional: when found, a
    5th variant builds with clang-cl.exe standing in for cl.exe -- same
    command-line flags/ABI/STL, but a different frontend, so it catches a
    different set of warnings. Skipped with a warning when it isn't installed.
    To add it, either install the "C++ Clang tools for Windows" component via
    the Visual Studio Installer, or install standalone LLVM:

      winget install LLVM.LLVM

    Usage:
      pwsh -File .\pre-push.ps1                    # build + test everything
      pwsh -File .\pre-push.ps1 multi_array_ref     # limit to one target/ctest filter
#>

[CmdletBinding()]
param(
    [Parameter(Position = 0)]
    [string]$Target
)

$ErrorActionPreference = 'Stop'

$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $repoRoot

function Write-Section([string]$Message) {
    Write-Host "`n=== $Message ===" -ForegroundColor Cyan
}

function Find-VsInstallation {
    $vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
    if (-not (Test-Path $vswhere)) { return $null }
    $installPath = & $vswhere -latest -products * `
        -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 `
        -property installationPath
    if (-not $installPath) { return $null }
    return $installPath.Trim()
}

function Find-CudaCompiler {
    $nvcc = Get-Command nvcc.exe -ErrorAction SilentlyContinue
    if ($nvcc) { return $nvcc.Source }
    if ($env:CUDA_PATH) {
        $candidate = Join-Path $env:CUDA_PATH 'bin\nvcc.exe'
        if (Test-Path $candidate) { return $candidate }
    }
    # Fall back to scanning the default install root directly (picking the
    # newest version by directory name) in case nvcc isn't on PATH and
    # CUDA_PATH wasn't set for this process (e.g. installed after this shell
    # started).
    $cudaRoot = Join-Path ${env:ProgramFiles} 'NVIDIA GPU Computing Toolkit\CUDA'
    if (Test-Path $cudaRoot) {
        $latest = Get-ChildItem $cudaRoot -Directory -Filter 'v*' -ErrorAction SilentlyContinue |
            Sort-Object Name -Descending | Select-Object -First 1
        if ($latest) {
            $candidate = Join-Path $latest.FullName 'bin\nvcc.exe'
            if (Test-Path $candidate) { return $candidate }
        }
    }
    return $null
}

function Find-ClangCl([string]$VsPath) {
    # Prefer the VS-bundled clang-cl (matches the cl.exe/STL version this
    # script otherwise pins) over a possibly-mismatched standalone install.
    $bundled = Join-Path $VsPath 'VC\Tools\Llvm\x64\bin\clang-cl.exe'
    if (Test-Path $bundled) { return $bundled }
    $onPath = Get-Command clang-cl.exe -ErrorAction SilentlyContinue
    if ($onPath) { return $onPath.Source }
    return $null
}

function Import-VsDevEnvironment([string]$VsPath) {
    $vsDevCmd = Join-Path $VsPath 'Common7\Tools\VsDevCmd.bat'
    if (-not (Test-Path $vsDevCmd)) { throw "VsDevCmd.bat not found at $vsDevCmd" }
    # Use %ComSpec% rather than a bare 'cmd.exe': the latter depends on the
    # invoking shell's PATH containing System32, which isn't always true
    # (trimmed/restricted PATH, unusual launch context, etc.), whereas
    # ComSpec is set by Windows itself to cmd.exe's full path unconditionally.
    $comSpec = if ($env:ComSpec) { $env:ComSpec } else { Join-Path $env:SystemRoot 'System32\cmd.exe' }
    $envDump = & $comSpec /c "`"$vsDevCmd`" -arch=x64 -host_arch=x64 -no_logo && set"
    foreach ($line in $envDump) {
        if ($line -match '^([A-Za-z_][A-Za-z0-9_]*)=(.*)$') {
            [System.Environment]::SetEnvironmentVariable($Matches[1], $Matches[2], 'Process')
        }
    }
}

# ---- locate toolchain ----------------------------------------------------
$vsPath = Find-VsInstallation
if (-not $vsPath) {
    Write-Error "No Visual Studio C++ toolset found via vswhere. Install the 'Desktop development with C++' workload."
    exit 1
}
Write-Section "Using Visual Studio at $vsPath"
Import-VsDevEnvironment $vsPath

$ninja = Join-Path $vsPath 'Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja\ninja.exe'

# Full path to the cl.exe that Import-VsDevEnvironment just put on PATH (from
# *this* $vsPath), so it can be pinned explicitly below instead of leaving
# CMake to cache whatever cl.exe/clang-cl.exe it first resolved. Otherwise, if
# a build dir was configured before a newer/older VS install appeared (e.g.
# after installing a preview alongside 2022), CMakeCache.txt keeps pointing at
# the old compiler binary while the environment's INCLUDE/LIB (re-imported
# fresh every run, above) point at the *new* install's MSVC STL headers --
# silent version-mismatch UB at best, a hard STL1000 static_assert at worst.
$clExe = (Get-Command cl.exe -ErrorAction SilentlyContinue).Source

$nvccExe = Find-CudaCompiler
if ($nvccExe) {
    Write-Host "nvcc detected at $nvccExe"
    # nvcc locates its sibling tools (ptxas.exe, cudafe++.exe, fatbinary.exe,
    # ...) via PATH rather than relative to its own full path, so pinning
    # $nvccExe alone isn't enough -- without this, the build fails with
    # "'ptxas' is not recognized as an internal or external command".
    $env:PATH = "$(Split-Path $nvccExe);$env:PATH"
} else {
    Write-Warning 'nvcc not found (checked PATH, CUDA_PATH, and the default install root); skipping the CUDA variant. Install it with: winget install --id Nvidia.CUDA --version 12.9'
}

$clangClExe = Find-ClangCl $vsPath
if ($clangClExe) {
    Write-Host "clang-cl detected at $clangClExe"
} else {
    Write-Warning 'clang-cl not found (checked the VS "C++ Clang tools for Windows" component and PATH); skipping the clang-cl variant. Install it via the Visual Studio Installer, or: winget install LLVM.LLVM'
}

if (Test-Path $ninja) {
    $env:PATH = "$(Split-Path $ninja);$env:PATH"
    $env:CMAKE_GENERATOR = 'Ninja'
    $env:NINJA_STATUS    = '[%f/%t||%r] '  # matches .gitlab-ci-correaa.yml's NINJA_STATUS
    Write-Host "Generator: Ninja ($ninja)"
} else {
    Write-Warning 'Bundled Ninja not found; falling back to the default Visual Studio generator (multi-config).'
}

$env:CMAKE_COLOR_DIAGNOSTICS = 'ON'

# ---- vcpkg (for find_package(Boost) to succeed) -----------------------------
$vcpkgToolchainArgs = @()
$vcpkgRoot = 'C:\vcpkg'
$vcpkgToolchain = Join-Path $vcpkgRoot 'scripts\buildsystems\vcpkg.cmake'
if (Test-Path $vcpkgToolchain) {
    $vcpkgToolchainArgs = @("-DCMAKE_TOOLCHAIN_FILE=$vcpkgToolchain")
    Write-Host "vcpkg detected at $vcpkgRoot; passing its toolchain file so find_package(Boost) can succeed."
} else {
    Write-Warning "vcpkg not found at $vcpkgRoot; Boost will likely not be found, so tests will not be built (see CMakeLists.txt's find_package(Boost) warning)."
}

$ccache = Get-Command ccache.exe -ErrorAction SilentlyContinue
if ($ccache) {
    $env:CMAKE_CXX_COMPILER_LAUNCHER = 'ccache'
    Write-Host 'ccache detected, using as compiler launcher.'
} else {
    Write-Warning 'ccache not found; builds will not be cached. (winget install ccache)'
}

# ---- parallelism -----------------------------------------------------------
$physCores = (Get-CimInstance Win32_Processor | Measure-Object -Property NumberOfCores -Sum).Sum
if (-not $physCores) { $physCores = [Environment]::ProcessorCount }
$env:CMAKE_BUILD_PARALLEL_LEVEL = "$physCores"
$env:CTEST_PARALLEL_LEVEL       = "$physCores"
$env:CTEST_OUTPUT_ON_FAILURE    = '1'

# ---- target / filter arg ----------------------------------------------------
$buildTargetArgs = @()
$ctestFilterArgs = @()
if ($Target) {
    Write-Host "Target/filter: $Target"
    $buildTargetArgs = @('--target', $Target)
    $ctestFilterArgs = @('-R', $Target)
} else {
    Write-Host 'No target arg; building/testing everything.'
}

$script:failed = @()
$script:noTests = @()

function Invoke-Variant {
    param(
        [string]$Name,
        [string]$BuildDir,
        [string[]]$ConfigureArgs,
        [string]$Config = 'Debug'
    )
    Write-Section $Name
    try {
        # Always reconfigure (cheap no-op when nothing changed) instead of only
        # configuring on first creation of $BuildDir: skipping it left stale
        # CMakeCache.txt entries (e.g. a compiler path from a VS install that's
        # since been superseded) silently un-synced with the freshly re-imported
        # environment above, see $clExe comment.
        # -Wno-unused-cli: on a reconfigure of an *existing* cache,
        # CMAKE_TOOLCHAIN_FILE is legitimately re-derived from the cache rather
        # than re-consumed (it only matters before the first project()), so
        # CMake would otherwise flag our always-passed -DCMAKE_TOOLCHAIN_FILE
        # as "unused" every run. (--no-warn-unused-cli is the deprecated spelling.)
        & cmake -Wno-unused-cli -S . -B $BuildDir @ConfigureArgs @vcpkgToolchainArgs
        if ($LASTEXITCODE -ne 0) { throw 'configure failed' }

        & cmake --build $BuildDir --config $Config @buildTargetArgs
        if ($LASTEXITCODE -ne 0) {
            # A clean re-run only needs to recompile the object(s) that hit a
            # transient failure, same rationale as the ctest retry below.
            Write-Warning "$Name`: build failed, retrying once."
            & cmake --build $BuildDir --config $Config @buildTargetArgs
            if ($LASTEXITCODE -ne 0) { throw 'build failed (twice)' }
        }

        $ctestOutput = & ctest --test-dir $BuildDir -C $Config @ctestFilterArgs 2>&1
        $ctestOutput | ForEach-Object { Write-Host $_ }
        if ($LASTEXITCODE -ne 0) {
            & ctest --test-dir $BuildDir -C $Config --rerun-failed @ctestFilterArgs
            if ($LASTEXITCODE -ne 0) { throw 'tests failed' }
        }
        if ($ctestOutput -match 'No tests were found') {
            Write-Warning "$Name`: no tests were actually compiled/run (likely CMakeLists.txt's find_package(Boost) failed to locate Boost -- see the 'Cannot find Boost' warning above). This variant only verified the headers compile, nothing was tested."
            $script:noTests += $Name
        }
    } catch {
        Write-Warning "$Name FAILED: $_"
        $script:failed += $Name
    }
}

# Pin cl.exe by full path (recomputed above from the current $vsPath) on every
# configure -- otherwise CMake just keeps whatever compiler path got cached
# the first time this build dir was configured, even after a newer VS install
# changes which cl.exe/headers Import-VsDevEnvironment puts on PATH.
$msvcCompilerArgs = @()
if ($clExe) {
    $msvcCompilerArgs = @("-DCMAKE_C_COMPILER=$clExe", "-DCMAKE_CXX_COMPILER=$clExe")
} else {
    Write-Warning 'cl.exe not found on PATH after importing the VS dev environment; letting CMake auto-detect the compiler (may pick up a stale cached one on reconfigure).'
}

# ---- variant 1: MSVC debug, warnings as errors ------------------------------
Invoke-Variant -Name 'MSVC (cl.exe) Debug -WX' -BuildDir '.build.msvc' -Config 'Debug' -ConfigureArgs (@(
    '-DCMAKE_BUILD_TYPE=Debug',
    '-DCMAKE_COMPILE_WARNING_AS_ERROR=ON'
) + $msvcCompilerArgs)

# ---- variant 2: MSVC + AddressSanitizer -------------------------------------
# RelWithDebInfo (not Debug) because CMake's default Debug flags add /RTC1,
# which MSVC refuses to combine with /fsanitize=address.
# /Zi + linker /DEBUG: without debug info, cl.exe emits warning C5072 ("ASAN
# enabled without debug information emission") and ASAN error reports fall
# back to raw addresses instead of symbolized file/line stack traces.
Invoke-Variant -Name 'MSVC (cl.exe) AddressSanitizer' -BuildDir '.build.msvc.asan' -Config 'RelWithDebInfo' -ConfigureArgs (@(
    '-DCMAKE_BUILD_TYPE=RelWithDebInfo',
    '-DCMAKE_COMPILE_WARNING_AS_ERROR=ON',
    '-DCMAKE_CXX_FLAGS=/fsanitize=address /Zi'
) + $msvcCompilerArgs)

# ---- variant 3: MSVC release, C++23 ------------------------------------------
Invoke-Variant -Name 'MSVC (cl.exe) Release C++23' -BuildDir '.build.msvc.c++23' -Config 'Release' -ConfigureArgs (@(
    '-DCMAKE_BUILD_TYPE=Release',
    '-DCMAKE_COMPILE_WARNING_AS_ERROR=ON',
    '-DCMAKE_CXX_STANDARD=23'
) + $msvcCompilerArgs)

# ---- variant 4: MSVC + nvcc CUDA --------------------------------------------
# Host compiler is cl.exe (via $msvcCompilerArgs and -DCMAKE_CUDA_HOST_COMPILER
# below), mirroring the Linux script's nvcc variants (pre-push:68) which use
# g++ as host compiler there. Skipped entirely when nvcc isn't installed.
# The cublas test subdirectory (include/boost/multi/adaptors/cuda/cublas/test)
# also requires a host BLAS (e.g. `vcpkg install openblas`); if that isn't
# found, its CMakeLists.txt warns and returns without adding any CUDA test
# executables, so this variant would compile but run no tests -- the same
# "no tests found" note the Boost-less variants above already surface.
if ($nvccExe) {
    $cudaHostCompilerArgs = @()
    if ($clExe) { $cudaHostCompilerArgs = @("-DCMAKE_CUDA_HOST_COMPILER=$clExe") }
    Invoke-Variant -Name 'MSVC (cl.exe) + nvcc CUDA' -BuildDir '.build.msvc.cuda' -Config 'Release' -ConfigureArgs (@(
        '-DCMAKE_BUILD_TYPE=Release',
        # Explicit ON (not just "unset"): this variant reuses .build.msvc.cuda
        # across runs, and CMake cache BOOLs persist across reconfigures, so
        # this pins the setting regardless of what an earlier ad-hoc configure
        # left cached, matching variant 1's intent of catching real warnings.
        # Requires every CUDA-Toolkit/Windows-SDK header warning under /Wall
        # to be suppressed per-directory (see the -Xcompiler=/wdNNNN lists in
        # fftw/test, thrust/test, cublas/test, etc.) so only warnings from our
        # own code turn into build failures here.
        '-DCMAKE_COMPILE_WARNING_AS_ERROR=ON',
        '-DENABLE_CUDA=1',
        "-DCMAKE_CUDA_COMPILER=$nvccExe",
        '-DCMAKE_CUDA_ARCHITECTURES=native',
        '-DCMAKE_CUDA_FLAGS=--extended-lambda -Wno-deprecated-gpu-targets'
    ) + $msvcCompilerArgs + $cudaHostCompilerArgs)
}

# ---- variant 5: clang-cl (LLVM) Debug, warnings as errors -------------------
# clang-cl is clang's MSVC-compatible driver: same command-line flags/ABI as
# cl.exe (so it still links against the MSVC STL/CRT and links via link.exe on
# PATH from the imported VS dev environment) but a different frontend and
# diagnostics engine, catching a different set of warnings/UB than cl.exe or
# nvcc's cl.exe host-compiler variants above. Skipped entirely when clang-cl
# isn't installed.
if ($clangClExe) {
    Invoke-Variant -Name 'clang-cl (LLVM) Debug -WX' -BuildDir '.build.clangcl' -Config 'Debug' -ConfigureArgs @(
        '-DCMAKE_BUILD_TYPE=Debug',
        '-DCMAKE_COMPILE_WARNING_AS_ERROR=ON',
        "-DCMAKE_C_COMPILER=$clangClExe",
        "-DCMAKE_CXX_COMPILER=$clangClExe",
        # thrust/omp/test resolves Thrust from the CUDA Toolkit's bundled
        # copy (nvcc's bin dir is on PATH for the rest of this script, once
        # the CUDA variant above has run), pulling in NVIDIA's vendored
        # libcu++. Its cstdlib does `using ::aligned_alloc;`, which clang-cl
        # doesn't expose in the global namespace -- fails the build. Not a
        # real regression to catch here, so skip it via the option this
        # CMakeLists.txt already exposes rather than special-casing clang-cl
        # in the CMake files themselves.
        '-DDISABLE_THRUST_OMP=ON'
    )
}

# ---- summary ------------------------------------------------------------------
if ($script:noTests.Count -gt 0) {
    Write-Host "`nNOTE: these variants compiled but ran no tests (Boost not found by CMake): $($script:noTests -join ', ')" -ForegroundColor Yellow
    Write-Host "      Install Boost so CMake can find it, e.g.: vcpkg install boost-serialization && cmake ... -DCMAKE_TOOLCHAIN_FILE=<vcpkg>/scripts/buildsystems/vcpkg.cmake" -ForegroundColor Yellow
}
if ($script:failed.Count -gt 0) {
    Write-Host "`nFAILED: $($script:failed -join ', ')" -ForegroundColor Red
    exit 666
} else {
    Write-Host "`nAll variants passed." -ForegroundColor Green
    exit 0
}
