.. _pseudo-handling:

Pseudopotential Handling in Nexus
=================================

.. dropdown:: For Developers
    :color: info
    :margin: auto auto 5 5
    :icon: code-square

    The majority of the functionality provided by :py:func:`~.generate_pseudoset` is in :py:class:`~.PseudoSet`.

.. tip::

    If you are migrating from :py:func:`~.ppset`, see :ref:`legacy-ppset`.

Nexus currently has several utilities for discovering, parsing, and copying pseudopotential files.
The latest and most complete method is through :py:func:`~.generate_pseudoset`, and is recommended for most users.

Each example has the output of the ``tree`` command in the pseudopotential directory to start.


.. _basic-pseudos:

Basic Pseudopotential Use
-------------------------

.. code-block:: none

    ccECP
    ├── C.ccECP.gamess
    ├── C.ccECP.nwchem
    ├── C.ccECP.upf
    ├── C.ccECP.xml
    ├── H.ccECP.gamess
    ├── H.ccECP.nwchem
    ├── H.ccECP.upf
    ├── H.ccECP.xml
    ├── O.ccECP.gamess
    ├── O.ccECP.nwchem
    ├── O.ccECP.upf
    └── O.ccECP.xml

    1 directory, 12 files

The most basic way to work with pseudopotentials in Nexus is to use ``settings`` to set ``pseudo_dir``, and then pass a list of file names into the relevant ``generate_`` functions.

.. code-block:: python

    from nexus import settings, run_project
    from nexus import generate_physical_system, generate_pwscf

    settings(
        pseudo_dir="/tmp/ccECP",
        ...
    )

    system = generate_physical_system(...)

    generate_pwscf(
        pseudos=["C.ccECP.upf", "H.ccECP.upf", "O.ccECP.upf"],
        system=system,
        ...
    )

    run_project()


.. _pseudoset-usage:

Using :py:func:`~.generate_pseudoset`
-------------------------------------

The variables returned by the methods described here can be passed into a generate function (e.g. :py:func:`~.generate_pwscf`), but require that you also pass the system to be simulated (created with :py:func:`~.generate_physical_system`) using the ``system`` keyword argument.

The following sections provide examples for how to use :py:func:`~.generate_pseudoset` for a variety of situations, primarily focused on the layout of a user's pseudopotential files.


Example 1 - Explicit File Name List
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: none

    ccECP
    ├── C.ccECP.gamess
    ├── C.ccECP.nwchem
    ├── C.ccECP.upf
    ├── C.ccECP.xml
    ├── H.ccECP.gamess
    ├── H.ccECP.nwchem
    ├── H.ccECP.upf
    ├── H.ccECP.xml
    ├── O.ccECP.gamess
    ├── O.ccECP.nwchem
    ├── O.ccECP.upf
    └── O.ccECP.xml

    1 directory, 12 files

Users familiar with the legacy :py:func:`~.ppset` or that have a custom pseudopotential directory can explicitly list the files they want out of their pseudopotential directory. The old style of passing ``pseudo_dir`` to ``settings`` will cause :py:func:`~.generate_pseudoset` to search in that directory for pseudopotentials.

.. important::

    The ``pseudo_dir`` given in ``settings`` is ignored if you pass ``pseudo_dir`` to :py:func:`~.generate_pseudoset`.

.. code-block:: python

    from nexus import settings, generate_pseudoset, run_project
    from nexus import generate_physical_system, generate_pwscf

    settings(pseudo_dir="/tmp/ccECP")

    ccECP = generate_pseudoset(
        qe      = ["C.ccECP.upf", "H.ccECP.upf", "O.ccECP.upf"],
        qmcpack = ["C.ccECP.xml", "H.ccECP.xml", "O.ccECP.xml"],
        pyscf   = ["C.ccECP.nwchem", "H.ccECP.nwchem", "O.ccECP.nwchem"],
        gamess  = ["C.ccECP.gamess", "H.ccECP.gamess", "O.ccECP.gamess"],
    )

    system = generate_physical_system(
        structure="structure.xsf",
        **ccECP["espresso"].get_Zeff({"C", "H", "O"})
    )

    generate_pwscf(
        pseudos=pseudos,
        system=system,
        ...,
    )

    run_project()


Example 2 - Search by Simulation Package Name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: none

    ccECP
    ├── C.ccECP.gamess
    ├── C.ccECP.nwchem
    ├── C.ccECP.upf
    ├── C.ccECP.xml
    ├── H.ccECP.gamess
    ├── H.ccECP.nwchem
    ├── H.ccECP.upf
    ├── H.ccECP.xml
    ├── O.ccECP.gamess
    ├── O.ccECP.nwchem
    ├── O.ccECP.upf
    └── O.ccECP.xml

    1 directory, 12 files

If you are using Nexus for just one code, e.g. driving high-throughput Quantum ESPRESSO calculations, the simplest option is to only get pseudos for QE.

.. tab-set::
    :sync-group: func-class-interface

    .. tab-item:: Function Interface
        :sync: func

        .. code-block:: python

            from nexus import settings, generate_pseudoset, run_project
            from nexus import generate_physical_system, generate_pwscf

            settings(pseudo_dir="/tmp/ccECP")

            pseudos = generate_pseudoset(code="quantum_espresso")

            system = generate_physical_system(
                structure="structure.xsf",
                **pseudos["espresso"].get_Zeff({"C", "H", "O"})
            )

            generate_pwscf(
                pseudos=pseudos,
                system=system,
                ...,
            )

            run_project()

    .. tab-item:: Class Interface
        :sync: class

        .. code-block:: python

            from nexus import PseudoSet, run_project
            from nexus import generate_physical_system, generate_pwscf

            pseudos = PseudoSet.from_dir(
                pseudo_dir="/tmp/ccECP",
                code="quantum_espresso",
            )

            # print(repr(pseudos)) # A single PseudoSet object
            # PseudoSet(
            #     codes = {'espresso'},
            #     pseudos = {
            #         'C': PosixPath('/tmp/ccECP/C.ccECP.upf'),
            #         'H': PosixPath('/tmp/ccECP/H.ccECP.upf'),
            #         'O': PosixPath('/tmp/ccECP/O.ccECP.upf'),
            #     },
            #     Zeff_map = {},
            # )

            system = generate_physical_system(
                structure="structure.xsf",
                **pseudos.get_Zeff({"C", "H", "O"})
            )

            generate_pwscf(
                pseudos=pseudos,
                system=system,
                ...,
            )

            run_project()


Example 3 - Filtering by File Extension
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are driving multiple codes with Nexus, e.g. running an RMG calculation for initial structure relaxation, Quantum ESPRESSO for a final relaxation calculation and orbital generation for use in a QMCPACK calculation, you can just grab all of the pseudos in the directory. This can cause problems with programs that can read several types of pseudopotentials (for example, RMG can use both ``.upf`` and ``.xml`` files); thus picking which format to use must be done manually. With :py:func:`~.generate_pseudoset`, this is done via the ``extension`` argument.

.. tab-set::
    :sync-group: func-class-interface

    .. tab-item:: Function Interface
        :sync: func

        .. code-block:: python

            from nexus import generate_pseudoset

            pseudos = generate_pseudoset(
                pseudo_dir="/tmp/ccECP",
                extension={"rmg": ".xml"},
            )

    .. tab-item:: Class Interface
        :sync: class

        .. code-block:: python

            from nexus import PseudoSet

            pseudos = PseudoSet.from_mixed_dir(
                pseudo_dir="/tmp/ccECP",
                extensions={"rmg": ".xml"}
            )


Example 4 - Searching by Inclusion Pattern
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: none

    /tmp/pseudo_dir
    ├── C.ccECP.gamess
    ├── C.ccECP.nwchem
    ├── C.ccECP.upf
    ├── C.ccECP.xml
    ├── C.NCPP.upf
    ├── C.NCPP.xml
    ├── C.USPP.upf
    ├── C.USPP.xml
    ├── H.ccECP.gamess
    ├── H.ccECP.nwchem
    ├── H.ccECP.upf
    ├── H.ccECP.xml
    ├── H.NCPP.upf
    ├── H.NCPP.xml
    ├── H.USPP.upf
    ├── H.USPP.xml
    ├── O.ccECP.gamess
    ├── O.ccECP.nwchem
    ├── O.ccECP.upf
    ├── O.ccECP.xml
    ├── O.NCPP.upf
    ├── O.NCPP.xml
    ├── O.USPP.upf
    └── O.USPP.xml

    7 directories, 30 files

.. note::

    Mixing filters, e.g. passing both ``extension`` and ``include`` will always add more filtering, never override another filter.

.. tab-set::
    :sync-group: func-class-interface

    .. tab-item:: Function Interface
        :sync: func

        .. code-block:: python

            from nexus import generate_pseudoset

            uspp = generate_pseudoset(
                pseudo_dir="/tmp/pseudo_dir",
                code="quantum_espresso",
                include="*USPP*",
            )

            ncpp = generate_pseudoset(
                pseudo_dir="/tmp/pseudo_dir",
                code="quantum_espresso",
                include="*NCPP*",
            )

            ccECP = generate_pseudoset(
                pseudo_dir="/tmp/pseudo_dir",
                code="quantum_espresso",
                include="*ccECP*",
            )

    .. tab-item:: Class Interface
        :sync: class

        .. code-block:: python

            from nexus import PseudoSet

            uspp = PseudoSet.from_dir(
                pseudo_dir="/tmp/pseudo_dir",
                code="quantum_espresso",
                include="*USPP*",
            )

            ncpp = PseudoSet.from_dir(
                pseudo_dir="/tmp/pseudo_dir",
                code="quantum_espresso",
                include="*NCPP*",
            )

            ccECP = PseudoSet.from_dir(
                pseudo_dir="/tmp/pseudo_dir",
                code="quantum_espresso",
                include="*ccECP*",
            )


Example 5 - Searching by Exclusion Pattern
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: none

    /tmp/vasp_pseudos/
    ├── C
    │   └── POTCAR
    ├── C_GW
    │   └── POTCAR
    ├── C_sv
    │   └── POTCAR
    ├── C_sv_GW
    │   └── POTCAR
    ├── H
    │   └── POTCAR
    ├── H_GW
    │   └── POTCAR
    ├── H_sv
    │   └── POTCAR
    ├── H_sv_GW
    │   └── POTCAR
    ├── O
    │   └── POTCAR
    ├── O_GW
    │   └── POTCAR
    ├── O_sv
    │   └── POTCAR
    └── O_sv_GW
        └── POTCAR

    13 directories, 12 files

A more specialized case involving VASP pseudopotentials requires the use of both the ``include`` and ``exclude`` parameters.

.. tab-set::
    :sync-group: func-class-interface

    .. tab-item:: Function Interface
        :sync: func

        .. code-block:: python

            from nexus import generate_pseudoset

            pseudos = generate_pseudoset(
                pseudo_dir="/tmp/vasp_pseudos",
                code="vasp",
                exclude="*_*",
            )
            sv_pseudos = generate_pseudoset(
                pseudo_dir="/tmp/vasp_pseudos",
                code="vasp",
                include="*_sv", # Leave out trailing asterisk to not match after 'sv'
            )
            sv_gw_pseudos = generate_pseudoset(
                pseudo_dir="/tmp/vasp_pseudos",
                code="vasp",
                include="*sv_GW",
            )
            gw_pseudos = generate_pseudoset(
                pseudo_dir="/tmp/vasp_pseudos",
                code="vasp",
                include="*_GW",
                exclude="*sv*",
            )

    .. tab-item:: Class Interface
        :sync: class

        .. code-block:: python

            from nexus import PseudoSet

            pseudos = PseudoSet.from_dir( # No `code` specified, uses auto-detect
                pseudo_dir="/tmp/vasp_pseudos",
                exclude="*_*",
            )
            sv_pseudos = PseudoSet.from_dir(
                pseudo_dir="/tmp/vasp_pseudos",
                include="*_sv", # Leave out trailing asterisk to not match after 'sv'
            )
            sv_gw_pseudos = PseudoSet.from_dir(
                pseudo_dir="/tmp/vasp_pseudos",
                include="*_sv_GW",
            )
            gw_pseudos = PseudoSet.from_dir(
                pseudo_dir="/tmp/vasp_pseudos",
                include="*_GW", # Include those ending with '_GW'
                exclude="*sv*", # But not those containing 'sv'
            )


It is important to note that the function-style interface does not support code autodetect, and if no codes are provided will create empty :py:class:`~.PseudoSet` objects for all non-VASP codes (in this case).


.. _legacy-ppset:

Migrating from :py:func:`~.ppset`
---------------------------------

The existing :py:func:`~.ppset` function is being deprecated for reasons described in its documentation.
Users are encouraged to switch to :py:func:`~.generate_pseudoset`, and try the features described above.
For those wishing for a path of least resistance migration, a short example is provided below.

.. tab-set::

    .. tab-item:: Legacy ``ppset``

        .. code-block:: python

            from nexus import settings, job, run_project, ppset
            from nexus import generate_physical_system, generate_structure
            from nexus import generate_pwscf
            from nexus import generate_pw2qmcpack, generate_qmcpack

            settings(
                pseudo_dir = "/path/to/pseudo_dir",
                ...
            )
            # Use a label
            ppset(
                label   = "ccECP",
                pwscf   = ["C.ccECP.upf", "H.ccECP.upf", "O.ccECP.upf"],
                qmcpack = ["C.ccECP.xml", "H.ccECP.xml", "O.ccECP.xml"],
            )

            system = generate_physical_system(...)

            scf = generate_pwscf(
                system  = system,
                pseudos = "ccECP",
                ...
            )

            nscf = generate_pwscf(
                system  = system,
                pseudos = "ccECP",
                ...
            )

            conv = generate_pw2qmcpack(...)

            opt = generate_qmcpack(
                system  = system,
                pseudos = "ccECP",
                ...
            )

            qmc = generate_qmcpack(
                system  = system,
                pseudos = "ccECP",
                ...
            )

            run_project()


    .. tab-item:: New ``generate_pseudoset``

        .. code-block:: python

            from nexus import settings, job, run_project, generate_pseudoset
            from nexus import generate_physical_system, generate_structure
            from nexus import generate_pwscf
            from nexus import generate_pw2qmcpack, generate_qmcpack

            settings(
                pseudo_dir = "/path/to/pseudo_dir",
                ...
            )
            # Assign to variable
            ccECP = generate_pseudoset(
                # You can also supply `pseudo_dir` here instead of in `settings`
                pwscf   = ["C.ccECP.upf", "H.ccECP.upf", "O.ccECP.upf"],
                qmcpack = ["C.ccECP.xml", "H.ccECP.xml", "O.ccECP.xml"],
            )

            system = generate_physical_system(...)

            scf = generate_pwscf(
                system  = system,
                pseudos = ccECP, # Pass variable instead of string
                ...
            )

            nscf = generate_pwscf(
                system  = system,
                pseudos = ccECP, # Pass variable instead of string
                ...
            )

            conv = generate_pw2qmcpack(...)

            opt = generate_qmcpack(
                system  = system,
                pseudos = ccECP, # Pass variable instead of string
                ...
            )

            qmc = generate_qmcpack(
                system  = system,
                pseudos = ccECP, # Pass variable instead of string
                ...
            )

            run_project()
