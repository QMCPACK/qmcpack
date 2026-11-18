.. _deepqmc-inference:

DeepQMC wavefunction inference
------------------------------

QMCPACK can evaluate a pretrained DeepQMC wavefunction through an embedded
Python/JAX bridge. This is an experimental, inference-only integration: it does
not train the model, optimize its parameters, or require a separate determinant
or QMCPACK SPOs. The bridge evaluates batches of walkers and uses the DeepQMC
ansatz log amplitude. It does not currently provide the sign/phase, so this
guide does not establish general nodal, complex-wavefunction, or DMC correctness.
The end-to-end recipe below was exercised for a small all-electron helium VMC
run on Linux Ubuntu 24.04 (aarch64), using GCC 14.2, Python 3.13.17, DeepQMC
1.2.0, and JAX/jaxlib 0.8.3. It is a smoke test, not a convergence or accuracy
guarantee. A complete reference package snapshot is supplied; other Python
versions and platforms require their own validation.

Python environment
~~~~~~~~~~~~~~~~~~

Use a native QMCPACK build environment with Python 3.13 development headers and
a shared embeddable Python library, in addition to the native dependencies
required by the standard :ref:`installation guide <obtaininginstalling>`. A
virtual environment does not supply missing Python development files. The
validated configuration used MPI because its HDF5 installation was parallel.

Set ``BASE_PYTHON`` to the absolute path of the Python 3.13 executable that has
those development files. The matching exact package versions are listed in
``examples/molecules/He/deepqmc-requirements.txt``. Create an isolated
environment and install the packages as follows, starting from the QMCPACK
source directory (``uv`` must be available):

.. code-block:: console

   $ export QMCPACK_SOURCE="$PWD"
   $ export BASE_PYTHON=/absolute/path/to/python3.13
   $ uv venv --python "$BASE_PYTHON" .venv-deepqmc
   $ export DEEPQMC_PYTHON="$PWD/.venv-deepqmc/bin/python"
   $ uv pip install --python "$DEEPQMC_PYTHON" \
       --only-binary :all: --no-binary antlr4-python3-runtime \
       --prerelease explicit -r examples/molecules/He/deepqmc-requirements.txt
   $ uv pip check --python "$DEEPQMC_PYTHON"
   $ "$DEEPQMC_PYTHON" -c 'import deepqmc, jax; from importlib.metadata import version; print(version("deepqmc"), jax.__version__, jax.devices())'

The sole source-build exception is ``antlr4-python3-runtime`` (version 4.9.3),
a pure-Python package; the other dependencies are installed from wheels. These
commands are scoped to the virtual environment and do not change system Python.
The exact ``tfp-nightly`` development release in the snapshot is intentional;
there is no need to allow prereleases globally or upgrade the dependency stack.

Building QMCPACK
~~~~~~~~~~~~~~~~

Use a fresh build directory and your normal compiler/library flags for the
native build. For an MPI+parallel-HDF5 environment, enable MPI as usual. For
example:

.. code-block:: console

   $ cmake -S . -B build/deepqmc-cpu \
       -DCMAKE_BUILD_TYPE=Release \
       -DENABLE_DEEPQMC_INFERENCE=ON \
       -DPython3_EXECUTABLE="$DEEPQMC_PYTHON" \
       -DQMC_MPI=ON -DBUILD_AFQMC=OFF \
       -DBUILD_UNIT_TESTS=ON -DBUILD_MICRO_BENCHMARKS=OFF
   $ cmake --build build/deepqmc-cpu \
       --target qmcpack test_wavefunction_trialwf -j 4

The option enables the embedded Python bridge; CMake must find both the
interpreter and development components and link ``Python3::Python``. Inspect the
CMake configure output/cache and confirm that the selected executable,
include directory, and library belong to the same Python installation as
``BASE_PYTHON``. Use that Python ABI at runtime. Keep the QMCPACK source tree
available: the bridge is imported from the source tree at runtime, and this
guide does not provide an installation/relocation mechanism for it.

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
   $ "$QMCPACK_SOURCE/build/deepqmc-cpu/bin/qmcpack" he_deepqmc.xml

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
     "$QMCPACK_SOURCE/build/deepqmc-cpu/src/QMCWaveFunctions/tests/test_wavefunction_trialwf" \
       'PythonDeepQMCBridge can load a real DeepQMC He checkpoint' -s

Keep ``DEEPQMC_HYDRA_CONFIG`` set to the matching configuration for this test.
With the environment variables set, the focused test checks the output shapes
and finite log-amplitude, gradient, and Laplacian values (12 assertions).
It does not compare against reference energies.

This bridge has been exercised only for helium VMC. Its ansatz/log-amplitude
interface and current Hessian evaluation should not be taken as validation for
other molecules, complex states, DMC, or production performance.
