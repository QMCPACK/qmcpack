.. _pwscf:

Working with PWSCF
==================

PWSCF Workflows
---------------

Allowed inflowing dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TBD

Available outflowing dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TBD


PWSCF Input
-----------

Generating input
^^^^^^^^^^^^^^^^

TBD

Reading and manipulating input
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TBD

Writing input
^^^^^^^^^^^^^

TBD


PWSCF Data Analysis
-------------------

Nexus reads the results of a completed Quantum ESPRESSO ``pw.x`` calculation
with a :class:`~pwscf_analyzer.PwscfAnalyzer`.  Its query functions provide a
consistent interface to structures, energies, k-points, electronic data,
forces, and stress.  When modern Quantum ESPRESSO XML and text output are
both available, modern XML is preferred and text output fills unavailable
quantities.  Some information, such as an ionic trajectory, is unique to text
output.

The following sections show how to load an analyzer and retrieve physical
quantities.  A query normally returns ``None`` if relevant data was not written
or could not be parsed.  It raises an exception if analysis was not performed,
the quantity is inapplicable to the detected run type, or an invalid unit is
requested.


Loading data from a simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a calculation created with :func:`~pwscf.generate_pwscf`, pass the
simulation object to :func:`~nexus.analyze_output`.  The factory loads a saved
analyzer image when present; otherwise it constructs and analyzes an analyzer.

.. code-block:: python

    scf = generate_pwscf(
        ...
        )
    ao  = analyze_output(scf)

``ao`` is a ``PwscfAnalyzer``.  This is convenient in a Nexus workflow because
the simulation supplies the calculation directory and file names.


Reading data directly
^^^^^^^^^^^^^^^^^^^^^

Output can also be analyzed without a simulation object.  Give the factory the
code name and either a calculation directory, input file, or output file.  Use
explicit file names when the directory contains nonstandard names.

.. code-block:: python

    ao = analyze_output(
        'pwscf',
        path         = './scf_run',
        infile_name  = 'scf.in',
        outfile_name = 'scf.out',
        )

The analyzer discovers modern Quantum ESPRESSO XML in the calculation
directory and reads it automatically when present.


Accessing physical quantities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The query methods return final values, rather than the histories printed during
iterative calculations.  Scalar quantities are Python ``float`` or ``bool``
objects; array quantities are NumPy ``float`` arrays; and structure queries
return :class:`~structure.Structure` objects.  Use the ``units`` argument where
applicable.

Supported run types and quantities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The table lists data ordinarily produced by each major PWSCF run type.  Actual
availability also depends on PWSCF printing options and successful parsing.
Here, *common data* means ``initial_structure``, ``kpoints``, and ``kweights``;
*electronic data* means ``eigenvalues``, ``occupations``, ``Ef``, ``Evbm``,
``Ecbm``, ``band_gap``, and ``fractional_occs``.

.. list-table:: Expected query data by PWSCF run type
    :header-rows: 1
    :widths: 18 48 34

    * - Run type
      - Quantities ordinarily available
      - Not ordinarily applicable
    * - ``scf``
      - Common data, ``energy``, electronic data, ``forces``, ``stress``, and
        ``pressure``.
      - ``relaxed_structure``.
    * - ``nscf``
      - Common data, ``energy``, and electronic data.
      - ``relaxed_structure``, ``forces``, ``stress``, and ``pressure``.
    * - ``bands``
      - Common data and ``eigenvalues``; ``occupations`` and band-edge data
        when reported.
      - ``energy``, ``relaxed_structure``, ``forces``, ``stress``, and
        ``pressure``.
    * - ``relax``
      - Common data, ``energy``, electronic data, ``relaxed_structure``,
        ``forces``, ``stress``, and ``pressure``.
      - None of the listed query groups.
    * - ``vc-relax``
      - The same data as ``relax``; the final cell is included in
        ``relaxed_structure``.
      - None of the listed query groups.
    * - ``md`` and ``vc-md``
      - Common data, ``energy``, final structure, ``forces``, ``stress``, and
        ``pressure``.  Electronic data is returned when reported.
      - No query group solely because it is a dynamics run.


Checking and requiring quantities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use :meth:`~pwscf_analyzer.PwscfAnalyzer.available` to check whether all named
quantities are accessible.  It returns a single Boolean and does not initiate
additional parsing.

.. code-block:: python

    electronic_data = ao.available(
        'eigenvalues',
        'occupations',
        'Ef',
        )
    print(electronic_data)

.. code-block:: text

    True

Use :meth:`~pwscf_analyzer.PwscfAnalyzer.require` when downstream work cannot
continue without a quantity.  A later query for an unavailable required
quantity raises an exception rather than returning ``None``.

.. code-block:: python

    ao.require(
        'energy',
        'forces',
        )

    if not ao.available('forces'):
        raise RuntimeError('PWSCF did not provide final ionic forces')


Common structure and k-point data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``initial_structure()`` returns a ``Structure`` with ``axes`` shaped ``(3, 3)``
and ``pos`` shaped ``(natoms, 3)``.  Its default length unit is Angstrom; use
``units='B'`` for bohr.

.. code-block:: python

    structure = ao.initial_structure()
    print(structure.axes.shape, structure.pos.shape)

.. code-block:: text

    (3, 3) (2, 3)

``kpoints()`` returns Cartesian reciprocal coordinates as a float array of
shape ``(nkpoints, 3)``.  Its default unit is inverse bohr; use ``units='A'``
for inverse Angstrom.

.. code-block:: python

    kpoints = ao.kpoints(units = 'A')
    print(kpoints.shape)

.. code-block:: text

    (8, 3)

``kweights()`` returns the dimensionless integration weights as a float array
of shape ``(nkpoints,)``.

.. code-block:: python

    kweights = ao.kweights()
    print(kweights)

.. code-block:: text

    [0.125 0.125 0.125 0.125 0.125 0.125 0.125 0.125]


Energy and electronic quantities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``energy()`` returns the final total energy as a Python ``float``.  The
default is Hartree; ``'eV'`` and ``'Ry'`` are also accepted.

.. code-block:: python

    energy = ao.energy(units = 'eV')
    print(f'E = {energy:.3f} eV')

.. code-block:: text

    E = -129.342 eV

``eigenvalues()`` returns a float array in eV by default, shaped
``(nkpoints, nbands)`` for non-spin calculations or
``(nkpoints, 2, nbands)`` for collinear spin calculations.

.. code-block:: python

    eigenvalues = ao.eigenvalues()
    print(eigenvalues.shape)

.. code-block:: text

    (8, 24)

``occupations()`` returns a dimensionless float array with the same shape as
``eigenvalues()``.

.. code-block:: python

    occupations = ao.occupations()
    print(occupations.shape)

.. code-block:: text

    (8, 24)

``Ef()``, ``Evbm()``, ``Ecbm()``, and ``band_gap()`` return Python ``float``
values in eV by default: the Fermi energy, valence maximum, conduction
minimum, and fundamental gap.  Each accepts ``'eV'``, ``'Ha'``, or ``'Ry'``.

.. code-block:: python

    fermi_energy = ao.Ef()
    vbm          = ao.Evbm()
    cbm          = ao.Ecbm()
    band_gap     = ao.band_gap()
    print(fermi_energy, vbm, cbm, band_gap)

.. code-block:: text

    5.931 4.802 6.427 1.625

``fractional_occs()`` returns a Python ``bool`` indicating whether an
occupation is neither empty nor full within a default tolerance of ``1e-3``.

.. code-block:: python

    metallic = ao.fractional_occs(tol = 1e-4)
    print(metallic)

.. code-block:: text

    False


Final structure, forces, and stress
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``relaxed_structure()`` returns the final ``Structure`` for a relaxation or
molecular-dynamics calculation.  Its axes and positions have the same shapes
as the initial structure.

.. code-block:: python

    relaxed_structure = ao.relaxed_structure(units = 'A')
    print(relaxed_structure.pos.shape)

.. code-block:: text

    (2, 3)

``forces()`` returns final Cartesian ionic forces as a float array of shape
``(natoms, 3)``.  The default unit is eV/Angstrom; ``'Ry/B'`` and ``'Ha/B'``
are also accepted.

.. code-block:: python

    forces = ao.forces(units = 'eV/A')
    print(forces.shape)

.. code-block:: text

    (2, 3)

``stress()`` returns the final stress tensor as a float array of shape
``(3, 3)``.  The default is GPa; it also accepts ``'Pa'``, ``'bar'``,
``'kbar'``, ``'Mbar'``, ``'atm'``, ``'eV/A^3'``, ``'Ha/Bohr^3'``, and
``'Ry/Bohr^3'``.

.. code-block:: python

    stress = ao.stress(units = 'GPa')
    print(stress.shape)

.. code-block:: text

    (3, 3)

``pressure()`` returns the final hydrostatic pressure as a Python ``float``.
It accepts the same units as ``stress()`` and uses GPa by default.

.. code-block:: python

    pressure = ao.pressure(units = 'GPa')
    print(f'P = {pressure:.3f} GPa')

.. code-block:: text

    P = 0.012 GPa
