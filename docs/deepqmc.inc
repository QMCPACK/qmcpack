.. _deepqmc-inference:

DeepQMC wavefunction inference (experimental)
---------------------------------------------

.. warning::

   Only the helium VMC workflow below has been validated. Sign/phase propagation
   is not implemented; DMC, nonlocal pseudopotential calculations, and
   off-diagonal estimators are not supported.

QMCPACK can evaluate a pretrained DeepQMC wavefunction through an embedded
Python/JAX bridge. This is an experimental, inference-only integration: it does
not train the model, optimize its parameters, it does not require or support
a separate determinant or combination with other QMCPACK WFC. An ebedded Python intepreter
combined with Python bridge code evaluates
batches of walkers using a trained checkpoint and a hydra defined DeepQMC ansatz log amplitude.
It does not currently propagate the sign/phase from deepQMC, so it does not yet
have general nodal, complex-wavefunction, or DMC correctness. At this time it is
intended only for batched VMC use with diagonal observables.

The end-to-end recipe below was exercised for a small all-electron helium VMC
run on Linux Ubuntu 24.04 (aarch64), using Open MPI 5.0.10 with Clang 21.1.4,
Python 3.13.17, DeepQMC 1.2.0, and JAX/jaxlib 0.8.3. It is a smoke test, not a
convergence or accuracy guarantee. A complete reference package snapshot is
supplied; other Python versions and platforms require their own validation.

Python environment
~~~~~~~~~~~~~~~~~~

Use a native QMCPACK build environment with Python 3.13 development headers and
a shared embeddable Python library, in addition to the native dependencies
required by the standard :ref:`installation guide <obtaininginstalling>`. A
virtual environment does not supply missing Python development files. The
validated configuration used MPI because its HDF5 installation was parallel.

Set ``BASE_PYTHON`` to the absolute path of the Python 3.13 executable that has
those development files. The matching exact package versions are listed in
``src/QMCWaveFunctions/DeepQMC/requirements.txt``. Create an isolated
environment and install the packages as follows, starting from the QMCPACK
source directory (``uv`` must be available):

.. code-block:: console

   $ export QMCPACK_SOURCE="$PWD"
   $ export BASE_PYTHON=/absolute/path/to/python3.13

Spack's environment view is a collection of links, not necessarily a complete
copy of each installation. Headers can be present under the Python installation
prefix even if they are absent from the view. Ask the selected interpreter for
its paths rather than looking only in the view's ``include`` directory:

.. code-block:: console

   $ "$BASE_PYTHON" - <<'PY'
   from pathlib import Path
   import sys
   import sysconfig

   print("version:", sys.version.split()[0])
   print("base prefix:", sys.base_prefix)
   header = Path(sysconfig.get_path("include")) / "Python.h"
   library = Path(sysconfig.get_config_var("LIBDIR")) / sysconfig.get_config_var("LDLIBRARY")
   for path in (header, library):
       print(path, "exists:", path.is_file())
   print("Py_ENABLE_SHARED:", sysconfig.get_config_var("Py_ENABLE_SHARED"))
   PY

Both files should exist, and ``Py_ENABLE_SHARED`` should be ``1`` for this
recipe. Ways to supply the development files are:

* **Spack Python:** if the files exist in its installation prefix, use them
  directly; there is no need to reinstall Python or repair the view. If they
  are genuinely missing, complete or rebuild that Spack installation before
  creating the venv. Installing packages with ``uv`` cannot supply these files.
* **Distribution Python:** Debian/Ubuntu ``python3-dev`` or Fedora
  ``python3-devel`` supplies development files for the distribution's Python;
  use a versioned package if the distribution provides the required version.
  These are not development files for a separate Spack Python. Do not combine,
  for example, system Python 3.12 headers with a Spack Python 3.13 interpreter.
* **Site-managed Python:** a complete CPython installation built with
  ``--enable-shared`` can also supply the interpreter, headers, and library.

If the files exist but CMake does not find them, obtain the base prefix with:

.. code-block:: console

   $ export PYTHON_ROOT="$("$BASE_PYTHON" -c 'import sys; print(sys.base_prefix)')"

Then add ``-DPython3_ROOT_DIR="$PYTHON_ROOT"`` to the configure command below,
alongside its explicit venv ``Python3_EXECUTABLE``. Use the base installation,
not the venv directory, for this hint. If setting ``Python3_INCLUDE_DIR`` and
``Python3_LIBRARY`` directly, both must belong to that same Python installation.

Once the development files are available, create the isolated environment and
install the packages:

.. code-block:: console

   $ uv venv --python "$BASE_PYTHON" .venv-deepqmc
   $ export DEEPQMC_PYTHON="$PWD/.venv-deepqmc/bin/python"
   $ uv pip install --python "$DEEPQMC_PYTHON" \
       --only-binary :all: --no-binary antlr4-python3-runtime \
       --prerelease explicit -r src/QMCWaveFunctions/DeepQMC/requirements.txt
   $ uv pip check --python "$DEEPQMC_PYTHON"
   $ "$DEEPQMC_PYTHON" -c 'import deepqmc, jax; from importlib.metadata import version; print(version("deepqmc"), jax.__version__, jax.devices())'

The sole source-build exception is ``antlr4-python3-runtime`` (version 4.9.3),
a pure-Python package; the other dependencies are installed from wheels. These
commands are scoped to the virtual environment and do not change system Python.
The exact ``tfp-nightly`` development release in the snapshot is intentional;
there is no need to allow prereleases globally or upgrade the dependency stack.

Building QMCPACK
~~~~~~~~~~~~~~~~

Activate the intended native build environment and use a fresh build directory.
The following MPI+parallel-HDF5 example explicitly selects the Open MPI compiler
wrappers, backed by Clang in the validated environment. ``QMC_MPI=ON`` enables
MPI support but does not select the compiler wrappers.

.. code-block:: console

   $ cmake -S . -B build/deepqmc-llvm \
       -DCMAKE_C_COMPILER="$(command -v mpicc)" \
       -DCMAKE_CXX_COMPILER="$(command -v mpicxx)" \
       -DCMAKE_BUILD_TYPE=Release \
       -DENABLE_DEEPQMC_INFERENCE=ON \
       -DPython3_EXECUTABLE="$DEEPQMC_PYTHON" \
       -DQMC_MPI=ON -DBUILD_AFQMC=OFF \
       -DBUILD_UNIT_TESTS=ON -DBUILD_MICRO_BENCHMARKS=OFF
   $ cmake --build build/deepqmc-llvm \
       --target qmcpack test_wavefunction_trialwf -j 4

The option enables the embedded Python bridge; CMake must find both the
interpreter and development components and link ``Python3::Python``. Inspect the
CMake configure output/cache and confirm that the selected executable,
include directory, and library belong to the same Python installation as
``BASE_PYTHON``. Use that Python ABI at runtime. Keep the QMCPACK source tree
available: the bridge is imported from the source tree at runtime, and this
guide does not provide an installation/relocation mechanism for it.

Reference helium test data
~~~~~~~~~~~~~~~~~~~~~~~~~~

The reference checkpoint and its matching Hydra configuration are distributed
separately from the source repository. Download the
`reference archive <https://app.box.com/s/26h4l97drt94xar8y01ix5jebars1j8v>`_
using Box's Download button. The matching
`SHA256 checksum file <https://app.box.com/s/d80v5j4e76g5d5z7r81uer15vdylrykc>`_
is also available on Box. These are sharing pages, not direct-download URLs;
CMake does not download the files.

Place the archive in the directory where you want to keep the test data.
Verify it against the checksum pinned in this repository before extracting:

.. code-block:: console

   $ sha256sum -c "$QMCPACK_SOURCE/tests/molecules/He_deepqmc/deepqmc-he-v1.tar.gz.sha256"
   $ mkdir -p HeECP
   $ tar -xzf deepqmc-he-v1.tar.gz -C HeECP
   $ export DEEPQMC_TEST_DATA_DIR="$PWD/HeECP"

The extracted directory has this layout:

.. code-block:: text

   HeECP/
     chkpt-30000.pt
     config.yaml

The archive may also contain ``hydra.yaml`` and ``overrides.yaml`` for training
provenance. The inference bridge uses ``config.yaml``, not ``hydra.yaml``.
Once the Python environment and build above are configured, enable the
reference tests with:

.. code-block:: console

   $ cmake -S "$QMCPACK_SOURCE" -B "$QMCPACK_SOURCE/build/deepqmc-llvm" \
       -DQMC_DEEPQMC_TEST_DATA_DIR="$DEEPQMC_TEST_DATA_DIR"
   $ cmake --build "$QMCPACK_SOURCE/build/deepqmc-llvm" \
       --target qmcpack test_wavefunction_trialwf -j 4
   $ ctest --test-dir "$QMCPACK_SOURCE/build/deepqmc-llvm" \
       -R '^deepqmc-he-' --output-on-failure

CMake checks both required files against the pinned SHA256 values in
``tests/molecules/He_deepqmc/SHA256SUMS`` before enabling the tests. A missing
file or mismatched checksum is a configuration error. This option identifies
one specific reference fixture; it is not a mechanism for selecting arbitrary
trained models. Leave it unset to build without reference-data tests.

The tests run actual checkpoint inference and the short helium VMC example,
then check that its scalar output is finite. They set the Python package path
and configuration explicitly; venv activation is not required. These checks
establish that inference executed, not scientific accuracy or convergence.
With ``BUILD_UNIT_TESTS=OFF``, only the application tests are available.

Run the helium example
~~~~~~~~~~~~~~~~~~~~~~

The XML input is ``examples/molecules/He/he_deepqmc.xml``. It describes a tiny
batched particle-by-particle VMC smoke test for two-electron all-electron
helium. Supply a trusted checkpoint trained for this system and its matching
Hydra configuration; model weights are not distributed with QMCPACK.
Checkpoints are Python pickle data: load only files from a trusted source,
since unpickling can execute code.

Run from a dedicated directory so the XML's relative ``model="chkpt.pt"`` is
resolved as expected. Copy the XML and checkpoint there, and point the
configuration override at the matching configuration file:

.. code-block:: console

   $ mkdir -p run-deepqmc-he
   $ cp examples/molecules/He/he_deepqmc.xml run-deepqmc-he/
   $ cp /trusted/path/chkpt-30000.pt run-deepqmc-he/chkpt.pt
   $ cd run-deepqmc-he
   $ export DEEPQMC_HYDRA_CONFIG=/trusted/path/config.yaml
   $ export DEEPQMC_PYTHON_SITE_PACKAGES="$("$DEEPQMC_PYTHON" -c 'import sysconfig; print(sysconfig.get_path("purelib"))')"
   $ export PYTHONPATH="$DEEPQMC_PYTHON_SITE_PACKAGES${PYTHONPATH:+:$PYTHONPATH}"
   $ export JAX_PLATFORMS=cpu OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
   $ "$QMCPACK_SOURCE/build/deepqmc-llvm/bin/qmcpack" he_deepqmc.xml

Successful completion produces ``He.deepqmc.s000.scalar.dat`` with finite
energies and a ``QMCPACK execution completed successfully`` message. Two blocks
and one walker are only an installation smoke test, not an energy estimate.

Use the absolute path to the QMCPACK executable. ``DEEPQMC_HYDRA_CONFIG``
selects the matching configuration explicitly; absent that override, the
bridge can discover ``.hydra/config.yaml`` in the checkpoint directory or an
ancestor. If no configuration is found, the bridge falls back to a default
helium ansatz that need not match your checkpoint; do not rely on that fallback.
The configuration is trusted input too: Hydra instantiation and OmegaConf
resolvers can execute code. ``PYTHONPATH`` is shown as an explicit way to expose packages to
the embedded interpreter; preserve any existing entries as in the command
above. ``python_module_path`` is another optional way to add a package
directory to ``sys.path``; it is not a Python executable or environment.
For example, set ``python_module_path="/absolute/path/to/.venv-deepqmc/lib/python3.13/site-packages"``
on the ``<deepqmc>`` element to expose this environment without ``PYTHONPATH``.
Venv activation also worked on the tested system, but cannot fix a mismatched
embedded Python ABI; use explicit package discovery when launching from a
batch job or a different shell. Do not set ``PYTHONHOME`` for this recipe.

The ``<deepqmc>`` element accepts ``model`` (required), ``source`` (ion particle
set, default ``ion0``), ``name`` (default ``deepqmc``), ``mol_idx`` (default
``0``), and optional ``python_module_path`` (a package directory). Coordinates
are in bohr; the configured molecule, ion/electron species, checkpoint, and
ansatz configuration must agree. This example uses the pretrained model as-is;
it performs no training or variational parameter optimization.

CPU is the validated backend here. JAX controls its own device selection,
independently of whether QMCPACK itself was built with GPU support. For GPU
experiments, install the appropriate JAX CUDA extra in a separate environment
using the official `JAX installation instructions <https://docs.jax.dev/en/latest/installation.html>`_,
then set ``JAX_PLATFORMS=cuda`` and verify ``jax.devices()`` reports the intended
GPU before running. GPU operation is not validated by this recipe; do not infer
that a QMCPACK GPU build is required merely to use JAX on a GPU. JAX computes a
full Hessian in this bridge, and first-use JIT compilation adds warm-up cost;
no general performance claim is implied.

Focused validation and limitations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The optional real-checkpoint unit test only performs inference when both
``DEEPQMC_HE_CHECKPOINT`` and ``DEEPQMC_PYTHON_SITE_PACKAGES`` are set; otherwise
it can report success without exercising a model. With a trusted matching
checkpoint and environment, run the focused test against the built executable:

.. code-block:: console

   $ DEEPQMC_HE_CHECKPOINT=/trusted/path/chkpt-30000.pt \
     DEEPQMC_PYTHON_SITE_PACKAGES="$DEEPQMC_PYTHON_SITE_PACKAGES" \
     "$QMCPACK_SOURCE/build/deepqmc-llvm/src/QMCWaveFunctions/tests/test_wavefunction_trialwf" \
       'PythonDeepQMCBridge can load a real DeepQMC He checkpoint' -s

Keep ``DEEPQMC_HYDRA_CONFIG`` set to the matching configuration for this test.
With the environment variables set, the focused test checks the output shapes
and finite log-amplitude, gradient, and Laplacian values (12 assertions).
It does not compare against reference energies.

This bridge has been exercised only for helium VMC. Its ansatz/log-amplitude
interface and current Hessian evaluation should not be taken as validation for
other molecules, complex states, DMC, or production performance.
