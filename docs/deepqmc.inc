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

.. _deepqmc-input-specification:

Input specification
~~~~~~~~~~~~~~~~~~~

DeepQMC is specified as the wavefunction component for the electron particle
set named by the enclosing ``wavefunction`` element's ``target`` attribute.
The component represents the complete pretrained wavefunction; do not add
Slater or Jastrow components. The ``deepqmc`` element has no child elements.

.. list-table:: Attributes of the ``deepqmc`` element
   :header-rows: 1
   :widths: 20 12 18 12 38

   * - Name
     - Datatype
     - Values
     - Default
     - Description
   * - ``name``
     - Text
     - Component name
     - ``deepqmc``
     - Name of this wavefunction component.
   * - ``source``
     - Text
     - Particle-set name
     - ``ion0``
     - Ion particle set supplying nuclear coordinates.
   * - ``model``
     - Text
     - File path
     - (required)
     - Trained checkpoint; relative paths are resolved from the run directory.
   * - ``config``
     - Text
     - File path
     - Beside checkpoint
     - Optional training configuration; defaults to ``.hydra/config.yaml``
       in the checkpoint directory.
   * - ``mol_idx``
     - Integer
     - Model-dependent
     - ``0``
     - Molecule index passed to the configured model.
   * - ``python_module_path``
     - Text
     - Directory path
     - None
     - Optional directory added to Python's ``sys.path``.

``mol_idx`` is a molecular configuration index for applicable models where a reference to a molecular
configuration is part of the model's inputs as determined by the ansatz the model uses, this is not applicable to all model ansatz and is ignored if not applicable. QMCPACK at this time does know whether this is necessary for a particular model's ansatz and so the applicability and bounds of the index are not checked.
The ``model`` file is a
full-wavefunction checkpoint, not an orbital or determinant file. The
``python_module_path`` attribute adds a directory to ``sys.path``; it does not
select a Python executable.

.. code-block:: xml
   :caption: DeepQMC wavefunction component

   <wavefunction name='psi0' target='e'>
     <deepqmc name='DNN' source='ion0' model='chkpt.pt' mol_idx='0'/>
   </wavefunction>

The optional ``config`` attribute selects a separately stored training
configuration:

.. code-block:: xml

   <deepqmc source="ion0" model="../model/chkpt-30000.pt"
            config="../model/config.yaml"/>

Explicit ``model`` and ``config`` paths may be absolute or relative to the
working directory from which QMCPACK is launched, not to the XML file.
Without ``config``, the bridge uses ``.hydra/config.yaml`` in the checkpoint
directory, following DeepQMC's conventional training output layout. A missing
configuration is an initialization error; no default molecule or ansatz is
substituted.

.. _deepqmc-system-matching:

Matching particle sets and Hamiltonian
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the standard ``MolecularHamiltonian``/``Molecule`` configuration shown
here, the resolved molecule used by the trained model needs to be determined and
ion and electron ParticleSets need to be defined to match its
particle count, spin, nuclear charges, positions, and ordering in the QMCPACK
input. Hydra configurations may instead use references, interpolation, named
molecules, or external definitions, so the relevant molecule may not appear as
an inline ``hamil.mol`` mapping in the YAML; you may need to inspect the resolved
configuration, not just the raw YAML. The hydra config.yaml initializes the ansatz; at this time we do no generate the necessary QMCPACK particle sets and Hamiltonian from it. The QMCPACK DeepQMC code does not cross check that they match and the behavior when they don't has not be investiaged.

For the supplied helium fixture, the following is an extract from the config.yaml covering the molecular configuration and hamiltonian which in this case is just the ansatz default plus an ecp_type.
usable configuration:

.. code-block:: yaml

   hamil:
     mol:
       coords: [[0, 0, 0]]
       charges: [2]
       charge: 0
       spin: 0
       unit: angstrom
     ecp_type: ccECP

``charges`` are nuclear atomic numbers, not net charges. Coordinate rows
correspond to nuclei in the same order as the XML ion particle set. DeepQMC
coordinates are converted from the configured unit to bohr; the bridge forwards
QMCPACK ion and electron coordinates without conversion, so supply them in
bohr. The XML ion particle set must have the same nuclei in the same order and
at the geometry used for this workflow. For an all-electron calculation, set
each ion species' ``charge`` to its nuclear ``Z`` and each electron species'
``charge`` to ``-1``.

For an all-electron system with total charge ``Q``, the electron count is
``N = sum(Z) - Q``. With spin ``S = Nup - Ndown``,
``Nup = (N + S)/2`` and ``Ndown = (N - S)/2``. Specify the QMCPACK electron
groups in DeepQMC's order: ``Nup`` spin-up electrons first, then ``Ndown``
spin-down electrons. The example uses group names ``u`` then ``d``; the bridge
passes coordinates in particle order, not according to the group names. For
the helium fixture, ``N=2`` and ``S=0``, so the XML has one electron in each group.

When an ECP is used, its valence electron counts (``ns_valence``) replace
nuclear ``Z`` in the electron-count relation; the actual ECP choice, mask, and
data must also be compatible. Nonlocal ECPs are not supported by this
implementation, so the fixture's ``ccECP`` setting is not an ECP usage recipe:
helium has ``Z=2`` and no active ECP mask, and the tested input is all-electron.
When ``hamil.ecp_type`` is ``None`` the model is all-electron. If it is set,
the default ECP mask selects charges greater than 2.

The corresponding QMCPACK helium particles, wavefunction, and Hamiltonian are
shown in the existing example:

.. literalinclude:: ../examples/molecules/He/he_deepqmc.xml
   :language: xml
   :start-at: <particleset name="ion0"
   :end-at: </hamiltonian>

That Hamiltonian includes electron-electron and electron-ion Coulomb terms and
uses QMCPACK's default electronic kinetic energy. A single nucleus has zero
ion-ion interaction; multi-ion systems require ion-ion repulsion and have not
been validated here.

Python environment
~~~~~~~~~~~~~~~~~~

In addtion to a native QMCPACK build environment for your platform, as described in :ref:`installation guide
<obtaininginstalling>`. using the DeepQMC wavefunction inference requires a complete
CPython 3.13 installation providing the Python executable, development headers, the Python standard
library (usually found in <prefix>/lib/python3.13/), and the shared libpython library (libpython3.13.so). A
virtual environment does not supply missing Python development files. Additionally
these instructions require the `uv` (https://github.com/astral-sh/uv) python package manager.
One way to fufill these requirements is by adding spack spec defs of:

.. code-block:: console

   - python@3.13
   - py-pip ^python@3.13
   - py-uv ^python@3.13

to your spack environment.

Once you have fulfilled these requirements an `uv` venv for deepQMC that will work with
the QMCPACK build can be created and checked with the following shell commands, they assume you are working
from the repo base directory:

.. code-block:: console

   $ export QMCPACK_SOURCE="$PWD"
   $ uv venv --python python3.13 .venv-deepqmc
   $ export DEEPQMC_PYTHON="$PWD/.venv-deepqmc/bin/python"
   $ uv pip install --python "$DEEPQMC_PYTHON" \
       --only-binary :all: --no-binary antlr4-python3-runtime \
       --prerelease explicit -r src/QMCWaveFunctions/DeepQMC/requirements.txt
   $ uv pip check --python "$DEEPQMC_PYTHON"


The ``python3.13`` command must select the intended CPython installation.
The exact package versions tested with are listed in
``src/QMCWaveFunctions/DeepQMC/requirements.txt``. This creates an isolated
python environment and installs the appropriate packages.

The sole source-build exception is ``antlr4-python3-runtime`` (version 4.9.3),
a pure-Python package; the other dependencies are installed from wheels. These
commands are scoped to the virtual environment and do not change system Python.
The exact ``tfp-nightly`` development release in the snapshot is intentional;
there is no need to allow prereleases globally or upgrade the dependency stack.

JAX backend: CPU or NVIDIA GPU
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The common install above provides the CPU JAX/jaxlib 0.8.3 packages. Choose
one backend in this same virtual environment before building or running. The
CPU route needs no extra packages:

.. code-block:: console

   $ export JAX_PLATFORMS=cpu

For the NVIDIA GPU route, install the CUDA 12 JAX extra in the same environment
(the constraint keeps the core requirement versions aligned), then explicitly
select CUDA:

.. code-block:: console

   $ uv pip install --python "$DEEPQMC_PYTHON" --only-binary :all: \
       -c "$QMCPACK_SOURCE/src/QMCWaveFunctions/DeepQMC/requirements.txt" 'jax[cuda12]==0.8.3'
   $ export JAX_PLATFORMS=cuda

The extra adds CUDA plugins and runtime packages; those additional packages
are not pinned by the CPU requirements snapshot. A compatible NVIDIA driver is
required, and the container must expose the GPUs. The CUDA runtime selected by
JAX is not necessarily the same version reported by the host driver; consult
the official `JAX installation instructions <https://docs.jax.dev/en/latest/installation.html>`_
for current compatibility guidance. After choosing either route, verify the
backend before proceeding. Import DeepQMC first so its package setup is included;
explicit ``JAX_PLATFORMS`` selection makes an unavailable requested backend an
error rather than silently falling back to CPU:

.. code-block:: console

   $ "$DEEPQMC_PYTHON" -c 'import deepqmc, jax; from importlib.metadata import version; print("deepqmc", version("deepqmc"), "jax", jax.__version__, "devices", jax.devices())'

The common CPU install may emit warnings about unavailable GPU plugins; that is
expected for CPU use. GPU execution requires the CUDA extra and a working
NVIDIA driver. QMCPACK's native ``QMC_GPU`` build options are independent and do
not install or select the JAX backend. The reference CTest workflow explicitly
requests CPU, even in an environment with GPU packages; validate GPU execution
with the backend/device check above and the manual helium run below.

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
include directory, and library belong to the base Python installation used to
create the venv. Use that Python ABI at runtime. Keep the QMCPACK source tree
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
     .hydra/
       config.yaml
       hydra.yaml
       overrides.yaml

The files ``hydra.yaml`` and ``overrides.yaml`` record training provenance.
The inference bridge uses ``.hydra/config.yaml``.
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
then check that its scalar output is finite. CMake supplies the fixture
directory and Python package path to CTest; venv activation is not required.
The configuration is read from the fixture's ``.hydra/config.yaml``. These checks
establish that inference executed, not scientific accuracy or convergence.
With ``BUILD_UNIT_TESTS=OFF``, only the application tests are available.

Run the helium example
~~~~~~~~~~~~~~~~~~~~~~

The XML input is ``examples/molecules/He/he_deepqmc.xml``. It describes a short
batched particle-by-particle VMC run for two-electron all-electron helium.
Supply a trusted checkpoint trained for this system and its matching Hydra
configuration; model weights are not distributed with QMCPACK.
Checkpoints are Python pickle data: load only files from a trusted source,
since unpickling can execute code.

Run from a dedicated directory so the XML's relative ``model="chkpt.pt"`` is
resolved as expected. Using the reference data extracted above, copy the XML,
checkpoint, and accompanying ``.hydra`` directory there:

.. code-block:: console

   $ mkdir -p run-deepqmc-he
   $ cp "$QMCPACK_SOURCE/examples/molecules/He/he_deepqmc.xml" run-deepqmc-he/
   $ cp "$DEEPQMC_TEST_DATA_DIR/chkpt-30000.pt" run-deepqmc-he/chkpt.pt
   $ cp -R "$DEEPQMC_TEST_DATA_DIR/.hydra" run-deepqmc-he/
   $ cd run-deepqmc-he
   $ export PYTHONPATH="$("$DEEPQMC_PYTHON" -c 'import sysconfig; print(sysconfig.get_path("purelib"))')${PYTHONPATH:+:$PYTHONPATH}"
   $ export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
   $ "$QMCPACK_SOURCE/build/deepqmc-llvm/bin/qmcpack" he_deepqmc.xml

Successful completion produces ``He.deepqmc.s000.scalar.dat`` with finite
energies and a ``QMCPACK execution completed successfully`` message. Two blocks
and one walker demonstrate execution, not a converged energy estimate.

Use the absolute path to the QMCPACK executable. This example uses the
checkpoint-local configuration layout described in
:ref:`deepqmc-input-specification`; no ``config`` attribute is needed.

For another system, supply its checkpoint and matching training configuration,
and adapt the QMCPACK XML particle sets and Hamiltonian as described in
:ref:`deepqmc-system-matching`.
Only the helium workflow above has been validated.

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

This example uses the pretrained model as-is; it performs no training or
variational parameter optimization.

Focused validation and limitations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Inference has been demonstrated with the supplied single-atom helium
checkpoint, whose training configuration specifies ``ccECP``. As described in
:ref:`deepqmc-system-matching`, no ECP is active for this helium system.
The integration is ready for testing with other checkpoints and larger systems
within the current VMC and observable limitations.

The CPU recipe was exercised on Linux Ubuntu 24.04 (aarch64), using Open MPI
5.0.10 with Clang 21.1.4, Python 3.13.17, DeepQMC 1.2.0, and JAX/jaxlib 0.8.3.
The GPU directions have also been exercised in a container on a host with an
NVIDIA A30 passed through and CUDA 12.6. Only CUDA 12.6 has been tested.

Only the reported Python and package versions have been tried. The package
versions are largely constrained by DeepQMC's requirements; the supplied
requirements file records the tested package versions.

The optional real-checkpoint unit test uses ``DEEPQMC_TEST_DATA_DIR`` to locate
``chkpt-30000.pt`` and its accompanying ``.hydra/config.yaml``. If the variable
is unset, the case reports success without exercising a model. With the Python
packages available as shown above, run it directly with one fixture-directory
variable:

.. code-block:: console

   $ DEEPQMC_TEST_DATA_DIR=/trusted/path/HeECP \
     "$QMCPACK_SOURCE/build/deepqmc-llvm/src/QMCWaveFunctions/tests/test_wavefunction_trialwf" \
       'PythonDeepQMCBridge can load a real DeepQMC He checkpoint' -s

When configured with ``QMC_DEEPQMC_TEST_DATA_DIR``, CTest supplies this variable
and the Python package path automatically. The focused test checks the output
shapes and finite log-amplitude, gradient, and Laplacian values.
It does not compare against reference energies.

These results demonstrate helium VMC inference, not validation of other
systems, complex states, DMC, or production performance.

Troubleshooting Python discovery
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If CMake cannot find Python's development files, ask the selected Python for
its base installation, header directory, and library directory:

.. code-block:: console

   $ python3.13 -c 'import sys; print(sys.base_prefix)'
   $ python3.13 -c 'import sysconfig; print(sysconfig.get_path("include"))'
   $ python3.13 -c 'import sysconfig; print(sysconfig.get_config_var("LIBDIR"))'

Look for ``Python.h`` and ``pyconfig.h`` in the header directory and
``libpython3.13.so`` in the library directory. Spack's view need not expose
these files if they are present in the underlying installation. In that case,
no reinstall or view repair is needed. Give CMake the base prefix as a hint:

.. code-block:: console

   $ cmake -S "$QMCPACK_SOURCE" -B "$QMCPACK_SOURCE/build/deepqmc-llvm" \
       -DPython3_ROOT_DIR="$(python3.13 -c 'import sys; print(sys.base_prefix)')"

If the files are genuinely missing, complete the Python installation before
creating the venv; installing packages with ``uv`` cannot supply them. A complete
Spack Python or site-managed CPython built with ``--enable-shared`` is suitable.
Distribution development packages, such as Debian/Ubuntu ``python3-dev`` or
Fedora ``python3-devel``, apply only to the corresponding distribution Python.
Do not combine those headers or libraries with a separate Spack interpreter.
