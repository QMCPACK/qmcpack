.. _user-guide-rmg:

Working with RMG
================

RMG Workflows
-------------

Allowed inflowing dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TBD

Available outflowing dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TBD


RMG Input
---------

Generating input
^^^^^^^^^^^^^^^^

TBD

Reading and manipulating input
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TBD

Writing input
^^^^^^^^^^^^^

TBD


RMG Data Analysis
-----------------

Nexus reads the results of a completed RMG calculation with an
:class:`~rmg_analyzer.RmgAnalyzer` object.  Its query functions provide a
consistent interface to structures, energies, k-points, electronic data,
forces, and stress.  RMG data analysis uses the calculation log output.

The following sections show how to load an analyzer and retrieve physical
quantities.  A query normally returns ``None`` if relevant data was not written
or could not be parsed.  It raises an exception if analysis was not performed,
the quantity is inapplicable to the detected run type, or an invalid unit is
requested.

RMG follows the same physical-quantity access patterns and return types shown
for :ref:`PWSCF data analysis <pwscf-data-analysis>`.


Loading data from a simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a calculation created with :func:`~rmg.generate_rmg`, pass the
simulation object to :func:`~nexus.analyze_output`.  The factory loads a saved
analyzer image when present; otherwise it constructs and analyzes an analyzer.

.. code-block:: python

    from nexus import analyze_output

    scf = generate_rmg(
        ...
        )
    ao  = analyze_output(scf)

``ao`` is an ``RmgAnalyzer`` object.  This is convenient in a Nexus workflow
because the simulation supplies all the information needed to perform the analysis.


Reading data directly
^^^^^^^^^^^^^^^^^^^^^

Output can also be analyzed without a simulation object.  Give
``analyze_output`` the RMG code name and its log-output path.  An input object,
input-file path, or input directory can be supplied separately.

.. code-block:: python

    from nexus import analyze_output

    ao = analyze_output(
        'rmg',
        input   = './scf_run/input',
        outfile = './scf_run/rmg.log',
        )

The analyzer reads the selected RMG log.  It can also locate and read the RMG
control input named by the log when that file is available.


Accessing physical quantities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The query methods return final values, rather than the histories printed during
iterative calculations.  Scalar quantities are Python ``float`` or ``bool``
objects; array quantities are NumPy ``float`` arrays; and structure queries
return :class:`~structure.Structure` objects.  Use the ``units`` argument where
applicable.

Supported run types and quantities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The table lists data ordinarily produced by each major RMG run type.  Actual
availability also depends on RMG printing options and successful parsing.
Here, *common data* means ``initial_structure``, ``kpoints``, and ``kweights``;
*electronic data* means ``eigenvalues``, ``occupations``, ``Ef``, ``Evbm``,
``Ecbm``, ``band_gap``, and ``fractional_occs``.

.. list-table:: Expected query data by RMG run type
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
    * - ``band``
      - Common data and ``eigenvalues``; ``occupations`` and band-edge data
        when reported.
      - ``energy``, ``relaxed_structure``, ``forces``, ``stress``, and
        ``pressure``.
    * - ``relax``
      - Common data, ``energy``, electronic data, ``relaxed_structure``,
        ``forces``, ``stress``, and ``pressure``.
      - None of the listed query groups.
    * - ``neb``
      - Common data, ``energy``, electronic data, ``relaxed_structure``,
        ``forces``, ``stress``, and ``pressure`` when reported.
      - None of the listed query groups.
    * - ``md_VE``, ``md_TE``, ``tddft``, ``exx``, and ``stm``
      - Data reported by the selected RMG mode.  Electronic quantities are
        available for ``md_VE``, ``md_TE``, and ``tddft`` when reported.
      - ``relaxed_structure`` except for ``neb`` and ``relax``.


Checking and requiring quantities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use :meth:`~rmg_analyzer.RmgAnalyzer.available` to check whether all named
quantities are accessible.  It returns a single Boolean and does not initiate
additional parsing.

.. code-block:: python

    electronic_data = ao.available('Ef','eigenvalues','occupations')
    print(electronic_data)

.. code-block:: text

    True

Use :meth:`~rmg_analyzer.RmgAnalyzer.require` when downstream work cannot
continue without a quantity.  A later query for an unavailable required
quantity raises an exception rather than returning ``None``.

.. code-block:: python

    print('Unavailable quantity:')
    print('Energy =',[ao.energy()])

    ao.require('energy','forces')

    print('\nNow required:')
    E = ao.energy()

.. code-block:: text

    Unavailable quantity:
    Energy = [None]

    Now required:
    RuntimeError: required RMG quantity "energy" is not available


Common structure and k-point data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``initial_structure()`` returns a ``Structure`` with ``axes`` shaped ``(3, 3)``
and ``pos`` shaped ``(natoms, 3)``.  Its default length unit is Angstrom; use
``units='B'`` for bohr.

.. code-block:: python

    structure = ao.initial_structure()

    print('cell axes:')
    print(structure.axes)

    print('\natomic positions:')
    print(structure.pos)

.. code-block:: text

    cell axes:
    [[2.500 0.000 0.000]
     [0.000 2.500 0.000]
     [0.000 0.000 2.500]]

    atomic positions:
    [[0.000 0.000 0.000]
     [1.250 1.250 1.250]]

``kpoints()`` returns Cartesian reciprocal coordinates as a float array of
shape ``(nkpoints, 3)``.  Its default unit is inverse bohr; use ``units='A'``
for inverse Angstrom.

.. code-block:: python

    kpoints = ao.kpoints(units = 'A')
    print(kpoints[:3])

.. code-block:: text

    [[ 0.000  0.000  0.000]
     [ 0.628  0.000  0.000]
     [ 0.000  0.628  0.000]]

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
    print(eigenvalues[0, :6])

.. code-block:: text

    [-15.382  -8.437  -1.928   2.714   6.183   9.762]

``occupations()`` returns a dimensionless float array with the same shape as
``eigenvalues()``.

.. code-block:: python

    occupations = ao.occupations()
    print(occupations[0, :6])

.. code-block:: text

    [1.000 1.000 1.000 0.000 0.000 0.000]

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

``relaxed_structure()`` returns the final ``Structure`` for a ``relax`` or
``neb`` calculation.  Its axes and positions have the same shapes as the
initial structure.

.. code-block:: python

    relaxed_structure = ao.relaxed_structure(units = 'A')
    print(relaxed_structure.pos)

.. code-block:: text

    [[0.000 0.000 0.000]
     [1.247 1.247 1.247]]

``forces()`` returns final Cartesian ionic forces as a float array of shape
``(natoms, 3)``.  The default unit is eV/Angstrom; ``'Ry/B'`` and ``'Ha/B'``
are also accepted.

.. code-block:: python

    forces = ao.forces(units = 'eV/A')
    print(forces)

.. code-block:: text

    [[-0.003  0.000  0.000]
     [ 0.003  0.000  0.000]]

``stress()`` returns the final stress tensor as a float array of shape
``(3, 3)``.  The default is GPa; it also accepts ``'Pa'``, ``'bar'``,
``'kbar'``, ``'Mbar'``, ``'atm'``, ``'eV/A^3'``, ``'Ha/Bohr^3'``, and
``'Ry/Bohr^3'``.

.. code-block:: python

    stress = ao.stress(units = 'GPa')
    print(stress)

.. code-block:: text

    [[ 0.010  0.000  0.000]
     [ 0.000  0.012  0.000]
     [ 0.000  0.000  0.014]]

``pressure()`` returns the final hydrostatic pressure as a Python ``float``.
It accepts the same units as ``stress()`` and uses GPa by default.

.. code-block:: python

    pressure = ao.pressure(units = 'GPa')
    print(f'P = {pressure:.3f} GPa')

.. code-block:: text

    P = 0.012 GPa


Strict and permissive parsing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RMG analysis reads one selected log-output file.  It may be provided directly,
or discovered as the sole ``*.log`` file under ``path``.

.. code-block:: python

    ao = analyze_output(
        'rmg',
        input   = 'input',
        outfile = 'rmg.log',
        path    = './scf_run',
        strict  = True,
        )

With the default ``strict=True``, an explicitly selected input or output must
exist, and output discovery must identify exactly one log.  If the input and
log both identify a top-level run mode, they must agree.  A disagreement raises
an exception under strict parsing.

Set ``strict=False`` to analyze an incomplete or missing log path without a
file error.  The analyzer completes with no parsed physical data, and
applicable queries return ``None``.  Any data that is successfully parsed from
an existing log remains available.

.. code-block:: python

    ao = analyze_output(
        'rmg',
        outfile = './interrupted_scf/rmg.log',
        strict = False,
        )

As with PWSCF, ``required`` quantities can be supplied at construction or
added with ``require()``.  They do not trigger additional parsing for RMG; a
query for unavailable required data raises an exception.
