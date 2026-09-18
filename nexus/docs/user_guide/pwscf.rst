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
with a :class:`~pwscf_analyzer.PwscfAnalyzer` object.  Its query functions 
provide a consistent interface to structures, energies, k-points, electronic 
data, forces, and stress.  When modern Quantum ESPRESSO XML and text output 
are both available, modern XML is preferred and text output fills unavailable
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

    from nexus import analyze_output

    scf = generate_pwscf(
        ...
        )
    ao  = analyze_output(scf)

``ao`` is a ``PwscfAnalyzer`` object.  This is convenient in a Nexus workflow 
because the simulation supplies all the information needed to perform the analysis.


Reading data directly
^^^^^^^^^^^^^^^^^^^^^

Output can also be analyzed without a simulation object.  Give `analyze_output` 
the name of the simulation code and either a calculation directory, input file, 
or output file.  Using explicit file names is recommended (see
:ref:`pwscf-source-selection` for alternatives).

.. code-block:: python

    from nexus import analyze_output

    ao = analyze_output(
        'pwscf',
        path         = './scf_run',
        infile_name  = 'scf.in',
        outfile_name = 'scf.out',
        )

    # or shortform:
    ao = analyze_output('pwscf','./scf_run','scf.in','scf.out')

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

    electronic_data = ao.available('Ef','eigenvalues','occupations')
    print(electronic_data)

.. code-block:: text

    True

Use :meth:`~pwscf_analyzer.PwscfAnalyzer.require` when downstream work cannot
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
    RuntimeError: required PWSCF quantity "energy" is not available


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

``relaxed_structure()`` returns the final ``Structure`` for a relaxation or
molecular-dynamics calculation.  Its axes and positions have the same shapes
as the initial structure.

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


.. _pwscf-source-selection:

Source selection and strict parsing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PWSCF results can be obtained from modern Quantum ESPRESSO schema XML, text
output, or both.  The ``source`` keyword selects the requested sources:
``'both'`` is the default, while ``'xml'`` and ``'out'`` select only modern XML
or text output, respectively.  Legacy ``data-file.xml`` is retained for
separate inspection but is not used by the physical-quantity queries.

The direct forms below show how :func:`~nexus.analyze_output` identifies the
calculation location and chooses candidate XML and text-output files.  The
keyword arguments shown after ``path`` are PWSCF-analyzer arguments.

.. code-block:: python

    ao = analyze_output(
        'pwscf',
        path   = './scf_run',
        source = 'both',
        strict = True,
        )

With only directory path provided (no input or output filenames), Nexus searches 
the directory for one ``*.out`` file and searches for one modern schema XML file 
named ``data-file-schema.xml`` within a ``*.save`` directory.  It checks both a
``*.save`` directory directly inside the calculation directory and one level
below it.

.. code-block:: python

    ao = analyze_output(
        'pwscf',
        path         = './scf_run',
        infile_name  = 'scf.in',
        outfile_name = 'scf.out',
        source       = 'both',
        strict       = True,
        )

An explicit ``outfile_name`` selects that text file directly.  An explicit
input file is also parsed.  When its ``CONTROL`` section specifies ``prefix``
and ``outdir``, those values identify the modern XML file at
``outdir/prefix.save/data-file-schema.xml``.  This input-derived location takes
precedence over the directory search.  If ``infile_name`` is given without
``outfile_name``, the output name defaults to the input-file stem with an
``.out`` suffix.

.. code-block:: python

    ao = analyze_output(
        'pwscf',
        path   = './scf_run/scf.in',
        source = 'xml',
        strict = True,
        )

Passing an input-file path is equivalent to providing its directory and
``infile_name``.  It enables input-guided XML discovery and infers
``scf.out`` as the text-output name if text parsing is selected.

.. code-block:: python

    ao = analyze_output(
        'pwscf',
        path   = './scf_run/scf.out',
        source = 'out',
        strict = True,
        )

Passing a text-output path selects that file directly.  No input file is
assumed, so XML discovery, if selected instead with ``source='both'`` or
``source='xml'``, uses the directory search rather than input-derived
``prefix`` and ``outdir`` values.

With ``strict=True``, all selected output files are required to exist when
analysis begins.  For ``source='xml'`` the modern XML file must exist; for
``source='out'`` the text-output file must exist.  For ``source='both'``, at
least one source must exist, so an XML-only or text-only calculation can be
analyzed without changing ``strict``.  A supplied input file is required.  
If more than one candidate modern XML or text-output file is found for a 
selected source, strict parsing raises an ambiguity error rather than choosing 
one arbitrarily.

Set ``strict=False`` to inspect an incomplete directory without file-discovery
errors.  Missing files are skipped.  An ambiguous search is likewise treated
as unavailable and is not parsed.  Applicable query functions then return
``None`` when their data was not found, while any successfully parsed data
remains accessible.

.. code-block:: python

    ao = analyze_output(
        'pwscf',
        path   = './interrupted_scf',
        source = 'both',
        strict = False,
        )

For ``source='both'``, modern XML is parsed first.  With no required
quantities, text output is parsed as well so that data unique to the log is
available.  If ``required`` quantities were supplied at construction, text
output is parsed only when XML does not provide all of them.  Regardless of
output file selection, a successfully parsed value is returned by its query 
even when it is unusual for the reconciled calculation type.
