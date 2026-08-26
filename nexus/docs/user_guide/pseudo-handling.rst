.. _pseudo-handling:

Pseudopotential Handling in Nexus
=================================

Nexus currently has several utilities for discovering, parsing, and copying pseudopotential files.
The latest and most complete method is through the :py:class:`~.PseudoSet` class and its associated builder :py:func:`~.generate_pseudoset`, and is recommended for most users.

All objects returned by the methods described here can be passed into a generate function (e.g. :py:func:`~.generate_pwscf`), but require that you also pass the system to be simulated as a :py:class:`~.PhysicalSystem` object using the ``system`` keyword argument, which you can create with :py:func:`~.generate_physical_system`.

.. _pseudoset-usage:

Using :py:func:`~.generate_pseudoset` and :py:class:`~.PseudoSet`
-----------------------------------------------------------------

The following sections provide examples for how to use :py:class:`~.PseudoSet` for a variety of situations, primarily focused on the layout of a user's pseudopotential files.
You can switch between a class-based interface and an interface based on :py:func:`~.generate_pseudoset` by clicking on each code block below.
Each example has the output of the ``tree`` command in the pseudopotential directory to start.

Example 1 - Structured Directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

            >>> from nexus import generate_pseudoset
            >>> pseudos = generate_pseudoset(
            ...     pseudo_dir="/tmp/ccECP",
            ...     code="espresso",
            ... )
            >>> print(pseudos) # A dict with a single key
            {'espresso': PseudoSet(
                codes = {'espresso'},
                pseudos = {
                    'C': PosixPath('/tmp/ccECP/C.ccECP.upf'),
                    'H': PosixPath('/tmp/ccECP/H.ccECP.upf'),
                    'O': PosixPath('/tmp/ccECP/O.ccECP.upf'),
                },
                Zeff_map = {},
            )}

    .. tab-item:: Class Interface
        :sync: class

        .. code-block:: python

            >>> from nexus import PseudoSet
            >>> pseudos = PseudoSet.from_dir(
            ...     pseudo_dir="/tmp/ccECP",
            ...     code="espresso",
            ... )
            >>> print(repr(pseudos)) # A single PseudoSet object
            PseudoSet(
                codes = {'espresso'},
                pseudos = {
                    'C': PosixPath('/tmp/ccECP/C.ccECP.upf'),
                    'H': PosixPath('/tmp/ccECP/H.ccECP.upf'),
                    'O': PosixPath('/tmp/ccECP/O.ccECP.upf'),
                },
                Zeff_map = {},
            )


If you are driving multiple codes with Nexus, e.g. running a Quantum ESPRESSO calculation to generate orbitals and then using them in a QMCPACK calculation, you can just grab all of the pseudos in the directory.

.. tab-set::
    :sync-group: func-class-interface

    .. tab-item:: Function Interface
        :sync: func

        .. code-block:: python

            >>> from nexus import generate_pseudoset
            >>> pseudos = generate_pseudoset(
            ...     pseudo_dir="/tmp/ccECP",
            ...     extension={"rmg": ".xml"}
            ... )
            >>> for code, ps_set in pseudos.items():
            ...     print("-"*40)
            ...     print(f"Code: {code}")
            ...     print(f"{ps_set!r}")
            ----------------------------------------
            Code: espresso
            PseudoSet(
                codes = {'espresso'},
                pseudos = {
                    'C': PosixPath('/tmp/ccECP/C.ccECP.upf'),
                    'H': PosixPath('/tmp/ccECP/H.ccECP.upf'),
                    'O': PosixPath('/tmp/ccECP/O.ccECP.upf'),
                },
                Zeff_map = {},
            )
            ----------------------------------------
            Code: gamess
            PseudoSet(
                codes = {'gamess'},
                pseudos = {
                    'C': PosixPath('/tmp/ccECP/C.ccECP.gamess'),
                    'H': PosixPath('/tmp/ccECP/H.ccECP.gamess'),
                    'O': PosixPath('/tmp/ccECP/O.ccECP.gamess'),
                },
                Zeff_map = {},
            )
            ----------------------------------------
            Code: pyscf
            PseudoSet(
                codes = {'pyscf'},
                pseudos = {
                    'C': PosixPath('/tmp/ccECP/C.ccECP.nwchem'),
                    'H': PosixPath('/tmp/ccECP/H.ccECP.nwchem'),
                    'O': PosixPath('/tmp/ccECP/O.ccECP.nwchem'),
                },
                Zeff_map = {},
            )
            ----------------------------------------
            Code: qmcpack
            PseudoSet(
                codes = {'qmcpack'},
                pseudos = {
                    'C': PosixPath('/tmp/ccECP/C.ccECP.xml'),
                    'H': PosixPath('/tmp/ccECP/H.ccECP.xml'),
                    'O': PosixPath('/tmp/ccECP/O.ccECP.xml'),
                },
                Zeff_map = {},
            )
            ----------------------------------------
            Code: rmg
            PseudoSet(
                codes = {'rmg'},
                pseudos = {
                    'C': PosixPath('/tmp/ccECP/C.ccECP.xml'),
                    'H': PosixPath('/tmp/ccECP/H.ccECP.xml'),
                    'O': PosixPath('/tmp/ccECP/O.ccECP.xml'),
                },
                Zeff_map = {},
            )
            ----------------------------------------
            Code: vasp
            PseudoSet(
                codes = {'vasp'},
                pseudos = {},
                Zeff_map = {},
            )

    .. tab-item:: Class Interface
        :sync: class

        .. code-block:: python

            >>> from nexus import PseudoSet
            >>> pseudos = PseudoSet.from_mixed_dir(
            ...     pseudo_dir="/tmp/ccECP",
            ...     extensions={"rmg": ".xml"}
            ... )
            >>> for code, ps_set in pseudos.items():
            ...     print("-"*40)
            ...     print(f"Code: {code}")
            ...     print(f"{ps_set!r}")
            ----------------------------------------
            Code: espresso
            PseudoSet(
                codes = {'espresso'},
                pseudos = {
                    'C': PosixPath('/tmp/ccECP/C.ccECP.upf'),
                    'H': PosixPath('/tmp/ccECP/H.ccECP.upf'),
                    'O': PosixPath('/tmp/ccECP/O.ccECP.upf'),
                },
                Zeff_map = {},
            )
            ----------------------------------------
            Code: gamess
            PseudoSet(
                codes = {'gamess'},
                pseudos = {
                    'C': PosixPath('/tmp/ccECP/C.ccECP.gamess'),
                    'H': PosixPath('/tmp/ccECP/H.ccECP.gamess'),
                    'O': PosixPath('/tmp/ccECP/O.ccECP.gamess'),
                },
                Zeff_map = {},
            )
            ----------------------------------------
            Code: pyscf
            PseudoSet(
                codes = {'pyscf'},
                pseudos = {
                    'C': PosixPath('/tmp/ccECP/C.ccECP.nwchem'),
                    'H': PosixPath('/tmp/ccECP/H.ccECP.nwchem'),
                    'O': PosixPath('/tmp/ccECP/O.ccECP.nwchem'),
                },
                Zeff_map = {},
            )
            ----------------------------------------
            Code: qmcpack
            PseudoSet(
                codes = {'qmcpack'},
                pseudos = {
                    'C': PosixPath('/tmp/ccECP/C.ccECP.xml'),
                    'H': PosixPath('/tmp/ccECP/H.ccECP.xml'),
                    'O': PosixPath('/tmp/ccECP/O.ccECP.xml'),
                },
                Zeff_map = {},
            )
            ----------------------------------------
            Code: rmg
            PseudoSet(
                codes = {'rmg'},
                pseudos = {
                    'C': PosixPath('/tmp/ccECP/C.ccECP.xml'),
                    'H': PosixPath('/tmp/ccECP/H.ccECP.xml'),
                    'O': PosixPath('/tmp/ccECP/O.ccECP.xml'),
                },
                Zeff_map = {},
            )
            ----------------------------------------
            Code: vasp
            PseudoSet(
                codes = {'vasp'},
                pseudos = {},
                Zeff_map = {},
            )


There are some important things to note here:

#. The extension must be specified for RMG.
    - This is because RMG can use both ``.upf`` and ``.xml`` files, and since a :py:class:`~.PseudoSet` can only have a single pseudopotential for each element you must specify which you wish to use for RMG, otherwise both will be found and an error will occur.
#. The VASP entry in the outputted dictionary has an empty :py:class:`~.PseudoSet`.
    - This is because there were no VASP-compatible pseudopotentials in the provided directory.
#. Unlike the single-code example which used :py:meth:`~.PseudoSet.from_dir`, the multi-code example uses :py:meth:`~.PseudoSet.from_mixed_dir`, which can automatically filter pseudos for specific codes.
    - The function :py:meth:`~.PseudoSet.from_dir` works well with directories containing pseudopotentials only for one code, e.g. if you have only a single set of ``.upf`` pseudopotentials. It can auto-detect the codes that *could* use those pseudopotentials (it can be more than one since, for example, both QE and RMG can also read UPF-formatted pseudopotentials)

Example 2 - Unstructured Directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Like in the previous example, we start with single-code workflows.

.. tab-set::
    :sync-group: func-class-interface

    .. tab-item:: Function Interface
        :sync: func

        .. code-block:: python

            >>> from nexus import generate_pseudoset
            >>> uspp = generate_pseudoset(
            ...     pseudo_dir="/tmp/pseudo_dir",
            ...     code="espresso",
            ...     include="*USPP*",
            ... )
            >>> print(uspp) # A dict with a single key
            {'espresso': PseudoSet(
                codes = {'espresso'},
                pseudos = {
                    'C': PosixPath('/tmp/pseudo_dir/C.USPP.upf'),
                    'H': PosixPath('/tmp/pseudo_dir/H.USPP.upf'),
                    'O': PosixPath('/tmp/pseudo_dir/O.USPP.upf'),
                },
                Zeff_map = {},
            )}
            >>> ncpp = generate_pseudoset(
            ...     pseudo_dir="/tmp/pseudo_dir",
            ...     code="espresso",
            ...     include="*NCPP*",
            ... )
            >>> print(ncpp) # A dict with a single key
            {'espresso': PseudoSet(
                codes = {'espresso'},
                pseudos = {
                    'C': PosixPath('/tmp/pseudo_dir/C.NCPP.upf'),
                    'H': PosixPath('/tmp/pseudo_dir/H.NCPP.upf'),
                    'O': PosixPath('/tmp/pseudo_dir/O.NCPP.upf'),
                },
                Zeff_map = {},
            )}
            >>> ccECP = generate_pseudoset(
            ...     pseudo_dir="/tmp/pseudo_dir",
            ...     code="espresso",
            ...     include="*ccECP*",
            ... )
            >>> print(ccECP) # A dict with a single key
            {'espresso': PseudoSet(
                codes = {'espresso'},
                pseudos = {
                    'C': PosixPath('/tmp/pseudo_dir/C.ccECP.upf'),
                    'H': PosixPath('/tmp/pseudo_dir/H.ccECP.upf'),
                    'O': PosixPath('/tmp/pseudo_dir/O.ccECP.upf'),
                },
                Zeff_map = {},
            )}

    .. tab-item:: Class Interface
        :sync: class

        .. code-block:: python

            >>> from nexus import PseudoSet
            >>> uspp = PseudoSet.from_dir(
            ...     pseudo_dir="/tmp/pseudo_dir",
            ...     code="espresso",
            ...     include="*USPP*",
            ... )
            >>> print(repr(uspp)) # A single PseudoSet object
            PseudoSet(
                codes = {'espresso'},
                pseudos = {
                    'C': PosixPath('/tmp/pseudo_dir/C.USPP.upf'),
                    'H': PosixPath('/tmp/pseudo_dir/H.USPP.upf'),
                    'O': PosixPath('/tmp/pseudo_dir/O.USPP.upf'),
                },
                Zeff_map = {},
            )
            >>> ncpp = PseudoSet.from_dir(
            ...     pseudo_dir="/tmp/pseudo_dir",
            ...     code="espresso",
            ...     include="*NCPP*",
            ... )
            >>> print(repr(ncpp)) # A single PseudoSet object
            PseudoSet(
                codes = {'espresso'},
                pseudos = {
                    'C': PosixPath('/tmp/pseudo_dir/C.NCPP.upf'),
                    'H': PosixPath('/tmp/pseudo_dir/H.NCPP.upf'),
                    'O': PosixPath('/tmp/pseudo_dir/O.NCPP.upf'),
                },
                Zeff_map = {},
            )
            >>> ccECP = PseudoSet.from_dir(
            ...     pseudo_dir="/tmp/pseudo_dir",
            ...     code="espresso",
            ...     include="*ccECP*",
            ... )
            >>> print(repr(ncpp)) # A single PseudoSet object
            PseudoSet(
                codes = {'espresso'},
                pseudos = {
                    'C': PosixPath('/tmp/pseudo_dir/C.ccECP.upf'),
                    'H': PosixPath('/tmp/pseudo_dir/H.ccECP.upf'),
                    'O': PosixPath('/tmp/pseudo_dir/O.ccECP.upf'),
                },
                Zeff_map = {},
            )


Multi-code workflows become more complex since we need to specify which files to include.
Conveniently, this is made trivial with the ``include`` parameter.

.. tab-set::
    :sync-group: func-class-interface

    .. tab-item:: Function Interface
        :sync: func

        .. code-block:: python

            >>> from nexus import generate_pseudoset
            >>> ccECP = generate_pseudoset(
            ...     pseudo_dir="/tmp/ccECP",
            ...     extension={"rmg": ".xml"}
            ... )
            >>> for code, ps_set in ccECP.items():
            ...     print("-"*40)
            ...     print(f"Code: {code}")
            ...     print(f"{ps_set!r}", end="")
            ------------------------------------------------------------
            Code: espresso
            PseudoSet(
                codes = {'espresso'},
                pseudos = {
                    'C': PosixPath('/tmp/pseudo_dir/C.ccECP.upf'),
                    'H': PosixPath('/tmp/pseudo_dir/H.ccECP.upf'),
                    'O': PosixPath('/tmp/pseudo_dir/O.ccECP.upf'),
                },
                Zeff_map = {},
            )
            ------------------------------------------------------------
            Code: gamess
            PseudoSet(
                codes = {'gamess'},
                pseudos = {
                    'C': PosixPath('/tmp/pseudo_dir/C.ccECP.gamess'),
                    'H': PosixPath('/tmp/pseudo_dir/H.ccECP.gamess'),
                    'O': PosixPath('/tmp/pseudo_dir/O.ccECP.gamess'),
                },
                Zeff_map = {},
            )
            ------------------------------------------------------------
            Code: pyscf
            PseudoSet(
                codes = {'pyscf'},
                pseudos = {
                    'C': PosixPath('/tmp/pseudo_dir/C.ccECP.nwchem'),
                    'H': PosixPath('/tmp/pseudo_dir/H.ccECP.nwchem'),
                    'O': PosixPath('/tmp/pseudo_dir/O.ccECP.nwchem'),
                },
                Zeff_map = {},
            )
            ------------------------------------------------------------
            Code: qmcpack
            PseudoSet(
                codes = {'qmcpack'},
                pseudos = {
                    'C': PosixPath('/tmp/pseudo_dir/C.ccECP.xml'),
                    'H': PosixPath('/tmp/pseudo_dir/H.ccECP.xml'),
                    'O': PosixPath('/tmp/pseudo_dir/O.ccECP.xml'),
                },
                Zeff_map = {},
            )
            ------------------------------------------------------------
            Code: rmg
            PseudoSet(
                codes = {'rmg'},
                pseudos = {
                    'C': PosixPath('/tmp/pseudo_dir/C.ccECP.xml'),
                    'H': PosixPath('/tmp/pseudo_dir/H.ccECP.xml'),
                    'O': PosixPath('/tmp/pseudo_dir/O.ccECP.xml'),
                },
                Zeff_map = {},
            )
            ------------------------------------------------------------
            Code: vasp
            PseudoSet(
                codes = {'vasp'},
                pseudos = {},
                Zeff_map = {},
            )

    .. tab-item:: Class Interface
        :sync: class

        .. code-block:: python

            >>> from nexus import PseudoSet
            >>> ccECP = PseudoSet.from_mixed_dir(
            ...     pseudo_dir = "/tmp/pseudo_dir",
            ...     extensions = {"rmg": ".xml"},
            ...     include    = "*ccECP*",
            ... )
            >>> for code, ps_set in ccECP.items():
            ...     print("-"*60)
            ...     print(f"Code: {code}")
            ...     print(f"{ps_set!r}", end="")
            ------------------------------------------------------------
            Code: espresso
            PseudoSet(
                codes = {'espresso'},
                pseudos = {
                    'C': PosixPath('/tmp/pseudo_dir/C.ccECP.upf'),
                    'H': PosixPath('/tmp/pseudo_dir/H.ccECP.upf'),
                    'O': PosixPath('/tmp/pseudo_dir/O.ccECP.upf'),
                },
                Zeff_map = {},
            )
            ------------------------------------------------------------
            Code: gamess
            PseudoSet(
                codes = {'gamess'},
                pseudos = {
                    'C': PosixPath('/tmp/pseudo_dir/C.ccECP.gamess'),
                    'H': PosixPath('/tmp/pseudo_dir/H.ccECP.gamess'),
                    'O': PosixPath('/tmp/pseudo_dir/O.ccECP.gamess'),
                },
                Zeff_map = {},
            )
            ------------------------------------------------------------
            Code: pyscf
            PseudoSet(
                codes = {'pyscf'},
                pseudos = {
                    'C': PosixPath('/tmp/pseudo_dir/C.ccECP.nwchem'),
                    'H': PosixPath('/tmp/pseudo_dir/H.ccECP.nwchem'),
                    'O': PosixPath('/tmp/pseudo_dir/O.ccECP.nwchem'),
                },
                Zeff_map = {},
            )
            ------------------------------------------------------------
            Code: qmcpack
            PseudoSet(
                codes = {'qmcpack'},
                pseudos = {
                    'C': PosixPath('/tmp/pseudo_dir/C.ccECP.xml'),
                    'H': PosixPath('/tmp/pseudo_dir/H.ccECP.xml'),
                    'O': PosixPath('/tmp/pseudo_dir/O.ccECP.xml'),
                },
                Zeff_map = {},
            )
            ------------------------------------------------------------
            Code: rmg
            PseudoSet(
                codes = {'rmg'},
                pseudos = {
                    'C': PosixPath('/tmp/pseudo_dir/C.ccECP.xml'),
                    'H': PosixPath('/tmp/pseudo_dir/H.ccECP.xml'),
                    'O': PosixPath('/tmp/pseudo_dir/O.ccECP.xml'),
                },
                Zeff_map = {},
            )
            ------------------------------------------------------------
            Code: vasp
            PseudoSet(
                codes = {'vasp'},
                pseudos = {},
                Zeff_map = {},
            )


As with the single-code example, you can change out ``include`` to match different pseudopotential sets.

Example 3 - VASP Pseudopotentials
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

            >>> from nexus import generate_pseudoset
            >>> pseudos = generate_pseudoset(
            ...     pseudo_dir="/tmp/vasp_pseudos",
            ...     code="vasp",
            ...     exclude="*_*",
            ... )
            >>> print(pseudos)
            {'vasp': PseudoSet(
                codes = {'vasp'},
                pseudos = {
                    'C': PosixPath('/tmp/vasp_pseudos/C/POTCAR'),
                    'H': PosixPath('/tmp/vasp_pseudos/H/POTCAR'),
                    'O': PosixPath('/tmp/vasp_pseudos/O/POTCAR'),
                },
                Zeff_map = {},
            )}
            >>> sv_pseudos = generate_pseudoset(
            ...     pseudo_dir="/tmp/vasp_pseudos",
            ...     code="vasp",
            ...     include="*_sv", # Leave out trailing asterisk to not match after 'sv'
            ... )
            >>> print(sv_pseudos)
            {'vasp': PseudoSet(
                codes = {'vasp'},
                pseudos = {
                    'C': PosixPath('/tmp/vasp_pseudos/C_sv/POTCAR'),
                    'H': PosixPath('/tmp/vasp_pseudos/H_sv/POTCAR'),
                    'O': PosixPath('/tmp/vasp_pseudos/O_sv/POTCAR'),
                },
                Zeff_map = {},
            )}
            >>> sv_gw_pseudos = generate_pseudoset(
            ...     pseudo_dir="/tmp/vasp_pseudos",
            ...     code="vasp",
            ...     include="*sv_GW",
            ... )
            >>> print(sv_gw_pseudos)
            {'vasp': PseudoSet(
                codes = {'vasp'},
                pseudos = {
                    'C': PosixPath('/tmp/vasp_pseudos/C_sv_GW/POTCAR'),
                    'H': PosixPath('/tmp/vasp_pseudos/H_sv_GW/POTCAR'),
                    'O': PosixPath('/tmp/vasp_pseudos/O_sv_GW/POTCAR'),
                },
                Zeff_map = {},
            )}
            >>> gw_pseudos = generate_pseudoset(
            ...     pseudo_dir="/tmp/vasp_pseudos",
            ...     code="vasp",
            ...     include="*_GW",
            ...     exclude="*sv*",
            ... )
            >>> print(gw_pseudos)
            {'vasp': PseudoSet(
                codes = {'vasp'},
                pseudos = {
                    'C': PosixPath('/tmp/vasp_pseudos/C_GW/POTCAR'),
                    'H': PosixPath('/tmp/vasp_pseudos/H_GW/POTCAR'),
                    'O': PosixPath('/tmp/vasp_pseudos/O_GW/POTCAR'),
                },
                Zeff_map = {},
            )}

    .. tab-item:: Class Interface
        :sync: class

        .. code-block:: python

            >>> from nexus import PseudoSet
            >>> pseudos = PseudoSet.from_dir( # No `code` specified, uses auto-detect
            ...     pseudo_dir="/tmp/vasp_pseudos",
            ...     exclude="*_*",
            ... )
            >>> print(repr(pseudos)) # All non-sv and non-GW pseudos
            PseudoSet(
                codes = {'vasp'},
                pseudos = {
                    'C': PosixPath('/tmp/vasp_pseudos/C/POTCAR'),
                    'H': PosixPath('/tmp/vasp_pseudos/H/POTCAR'),
                    'O': PosixPath('/tmp/vasp_pseudos/O/POTCAR'),
                },
                Zeff_map = {},
            )
            >>> sv_pseudos = PseudoSet.from_dir(
            ...     pseudo_dir="/tmp/vasp_pseudos",
            ...     include="*_sv", # Leave out trailing asterisk to not match after 'sv'
            ... )
            >>> print(repr(sv_pseudos))
            PseudoSet(
                codes = {'vasp'},
                pseudos = {
                    'C': PosixPath('/tmp/vasp_pseudos/C_sv/POTCAR'),
                    'H': PosixPath('/tmp/vasp_pseudos/H_sv/POTCAR'),
                    'O': PosixPath('/tmp/vasp_pseudos/O_sv/POTCAR'),
                },
                Zeff_map = {},
            )
            >>> sv_gw_pseudos = PseudoSet.from_dir(
            ...     pseudo_dir="/tmp/vasp_pseudos",
            ...     include="*_sv_GW",
            ... )
            >>> print(repr(sv_gw_pseudos))
            PseudoSet(
                codes = {'vasp'},
                pseudos = {
                    'C': PosixPath('/tmp/vasp_pseudos/C_sv_GW/POTCAR'),
                    'H': PosixPath('/tmp/vasp_pseudos/H_sv_GW/POTCAR'),
                    'O': PosixPath('/tmp/vasp_pseudos/O_sv_GW/POTCAR'),
                },
                Zeff_map = {},
            )
            >>> gw_pseudos = PseudoSet.from_dir(
            ...     pseudo_dir="/tmp/vasp_pseudos",
            ...     include="*_GW", # Include those ending with '_GW'
            ...     exclude="*sv*", # But not those containing 'sv'
            ... )
            >>> print(repr(gw_pseudos))
            PseudoSet(
                codes = {'vasp'},
                pseudos = {
                    'C': PosixPath('/tmp/vasp_pseudos/C_GW/POTCAR'),
                    'H': PosixPath('/tmp/vasp_pseudos/H_GW/POTCAR'),
                    'O': PosixPath('/tmp/vasp_pseudos/O_GW/POTCAR'),
                },
                Zeff_map = {},
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
