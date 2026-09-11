.. raw:: html

    <style> .red {color:red; font-weight:bold} </style>
    <style> .yellow {color:yellow; font-weight:bold} </style>
    <style> .green {color:green; font-weight:bold} </style>
    <style> .lred {color:red; font-weight:bold; font-size:x-large} </style>
    <style> .lgreen {color:green; font-weight:bold; font-size:x-large} </style>

.. role:: red
.. role:: yellow
.. role:: green
.. role:: lred
.. role:: lgreen

.. _code-style:


Code Style for Developers
=========================

This section outlines requirements and recommendations regarding new code developed in Nexus. Please consult and follow this documentation prior to submitting pull requests on `GitHub <https://github.com/QMCPACK/qmcpack>`__. This will save a lot of time and effort during code reviews.


.. _text-style:

Text style
----------

Class names are **CamelCase**: ``class BaseClass:``.

All other code is **snake_case**: ``def perform_work():``.

**CamelCase** means the first letter of each word is capitalized with no separation between words. **snake_case** means each word is lowercase with single underscores separating each word.


.. _writing-strings:

Writing Strings
---------------

Split long strings with newlines on the newline characters:

:lred:`No`

.. code-block:: python

    raise ValueError("This is a really super duper long error message.\nIt will display as many lines\nand realistically there is no good reason for this to be\nover 150 characters long on one line")

:lgreen:`Yes`

.. code-block:: python

    raise ValueError(
        "This is a really super duper long error message.\n"
        "It will display as many lines\n"
        "and realistically there is no good reason for this to be\n"
        "over 150 characters long on one line"
        )

It is preferable to keep strings that do not have newlines on a single line so that users can use tools like ``grep`` to search for them in source code.

:lred:`No`

.. code-block:: python

    raise ValueError(
        "This is a string that is a single line but is actually quite "
        "long and may extend well past the suggested line length limits of Nexus"
        )

:lgreen:`Yes`

.. code-block:: python

    raise ValueError(
        "This is a string that is a single line but is actually quite long and may extend well past the suggested line length limits of Nexus"
        )

This does not always apply for f-strings or formatted strings, which are generally unable to be grepped regardless of if they are split across multiple lines.

:lred:`No`

.. code-block:: python

    raise ValueError(f"This is an f-string that has {included} substitutions but is also a single line and is {greater_than} the suggested line length limits of Nexus")

:lgreen:`Yes`

.. code-block:: python

    raise ValueError(
        f"This is single line that is {greater_than} the suggested line length "
        f"limits of Nexus {included} but "
        f"has substitutions that make it {hard} to grep."
        )

Prefer `f-strings <https://docs.python.org/3/reference/lexical_analysis.html#f-strings>`__ over :py:meth:`str.format()`:

:lred:`No`

.. code-block:: python

    msg += (
        "Values are not consistent.\n"
        "Present thing:\n"
        "{0}\n"
        "Given value:\n"
        "{1}\n"
        "Consistent value:\n"
        "{2}\n"
        "Absolute difference: {3}\n".format(
            present, given, consistent, diff
            )
        )

:lgreen:`Yes`

.. code-block:: python

    msg += (
        "Values are not consistent.\n"
        "Present thing:\n"
        f"{present}\n"
        "Given value:\n"
        f"{given}\n"
        "Consistent value:\n"
        f"{consistent}\n"
        f"Absolute difference: {diff}\n"
        )

Avoid putting non-minimal code inside an f-string, prefer to assign to a variable beforehand.

:lred:`No`

.. code-block:: python

    msg += (
        f"Present value: {sum([i for i in thing])} does not match given {sum([i for i in given])}\n"
        )

:lgreen:`Yes`

.. code-block:: python

    present = sum([i for i in thing])
    given = sum([i for i in given])
    msg += (
        f"Present value: {present} does not match given {given}\n"
        )


.. _variable-names:

Variable and function names
---------------------------

Be descriptive and generally limit names to three words or less for legibility.

:lred:`No`

.. code-block:: python

    def slice_dice_and_mince_data(data):
        ...

:lgreen:`Yes`

.. code-block:: python

    def slice_data(data):
        ...

:lred:`No`

.. code-block:: python

    total_amount_of_money_in_bag = total_number_of_pennies_in_bag + number_of_pennies_in_a_nickel*total_number_of_nickels_in_bag + number_of_pennies_in_a_dime*total_number_of_nickels_in_bag + ...

:lgreen:`Yes`

.. code-block:: python

    money_tot = pennies + 5*nickels + 10*dimes + ...


.. _encapsulation:

Encapsulation
-------------

Function and class endings are delimited with a comment:

.. code-block:: python

    def perform_work():
        ...
    #end def perform_work

    class BaseClass:
        ...
    #end class BaseClass

Enclosing comments such as ``#end if`` and ``#end for`` are not required, but are encouraged in cases of deep nesting for readability.

For closing parentheses (and similar), the closing character should have the same indentation as the function arguments or e.g. list elements.

:lred:`No`

.. code-block:: python

    def operate(key1 = 'key1',
                key2 = 'key2',
                key3 = 'key3',
    ):
        key = (key1,key2,key3)
        return key
    #end def operate

:lgreen:`Yes`

.. code-block:: python

    def operate(key1 = 'key1',
                key2 = 'key2',
                key3 = 'key3',
                ):
        key = (key1,key2,key3)
        return key
    #end def operate

:lred:`No`

.. code-block:: python

    d = dict(a = 1,
             b = 2,
    )

:lgreen:`Yes`

.. code-block:: python

    d = dict(a = 1,
             b = 2,
             )

    # or

    d = dict(a = 1,
             b = 2)

:lred:`No`

.. code-block:: python

    l = ['key1',
         'key2',
         'key3',
    ]


:lgreen:`Yes`

.. code-block:: python

    l = ['key1',
         'key2',
         'key3',
         ]

    # or

    l = ['key1',
         'key2',
         'key3']

:lred:`No`

.. code-block:: python

    def operate(argument1='arg1',other_thing='other',
                key='key',value='bar',extra_thing='bean',
                other_extra='bean2'):
        ...
    #end def operate

:lgreen:`Yes`

.. code-block:: python

    def operate(
        argument1   = 'arg1',
        other_thing = 'other',
        key         = 'key',
        value       = 'bar',
        extra_thing = 'bean',
        other_extra = 'bean2',
        ):
        ...
    #end def operate


.. _opening-files:

Opening Files
-------------

The preferred way to open files in Nexus is with a `context manager <https://docs.python.org/3/reference/compound_stmts.html#the-with-statement>`__.
This prevents accidentally leaving a file open if you don't explicitly close it.

:lred:`No`

.. code-block:: python

    fobj = open("/path/to/file", "r")
    text = fobj.read()
    fobj.close()

:lgreen:`Yes`

.. code-block:: python

    with open("/path/to/file", "r") as fobj:
        text = fobj.read()


.. _function-arguments:

Function arguments
------------------

Use keyword arguments when calling a function, unless they are few.

You should require that all boolean-type arguments are keyword-only, done by `placing an asterisk argument <https://peps.python.org/pep-3102/>`__ before the boolean argument.
This prevent's mistakes due to `Python's "truthiness" evaluation <https://docs.python.org/3/library/stdtypes.html#truth-value-testing>`__.

:lred:`No`

.. code-block:: python

    def func(arg0: int, arg1: bool, arg2: int):
        ...
    #end def func

:lgreen:`Yes`

.. code-block:: python

    def func(arg0: int, *, arg1: bool, arg2: int):
        ...
    #end def func

Functions with arbitrary arguments are allowed, but generally discouraged if there is not a good reason.

:lred:`No`

.. code-block:: python

    def func(**kwargs):
        thing1 = kwargs.get("thing1", 0)
        thing2 = kwargs.get("thing2", 0)
        return thing1 + thing2
    #end def func

:lgreen:`Yes`

.. code-block:: python

    def func(
        thing1 = 0,
        thing2 = 0,
        ):
        return thing1 + thing2
    #end def func

You should use kwargs when calling another function whose arguments may vary over time, or are excessively long.
This makes the outer function less prone to becoming out of date if the inner function changes.

:lred:`No`

.. code-block:: python

    def func_calling_func(mine1, mine2, other1, other2, ..., otherN):
        my_kwargs = [mine1, mine2]
        other_stuff = other_func(other1, other2, ..., otherN)
        return my_kwargs, other_stuff
    #end def func_calling_func

:lgreen:`Yes`

.. code-block:: python

    def func_calling_func(arg0, **kwargs):
        # Pull out your arguments, which depend on arg0
        my_kwargs = []
        for k in kwargs.items():
            if k in get_args(arg0):
                my_kwargs.append(kwargs.pop(k))
        # Remaining things in kwargs get passed to other_func
        other_stuff = other_func(**kwargs)
        return my_kwargs, other_stuff
    #end def func_calling_func


.. _classes:

Classes
-------

Avoid large/complex classes (yes, this is ironic given some of the early Nexus code).

Use tight procedural programming where possible.

Use "heavy" base classes to promote light derived classes.

Avoid set/get accessors.

New classes in Nexus should always inherit from :py:class:`DevBase`.

:lred:`No`

.. code-block:: python

    class MyClass:
        ...

:lgreen:`Yes`

.. code-block:: python

    from .developer import DevBase

    class MyClass(DevBase):
        ...


Class variables should be immutable if they are not dynamic.
If you have a dynamic class variable, annotate it with ``typing.ClassVar`` to make sure a linter doesn't flag it.

:lred:`No`

.. code-block:: python

    class MyClass(DevBase):
        unique_vals = {"a", "b", "c"} # sets are mutable
        sequence = [1, 2, 3] # lists are mutable
        mapping = {"a": 1, "b": 2} # dicts are mutable
        class_storage = obj() # Is this a class variable or an instance default?

:lgreen:`Yes`

.. code-block:: python

    from types import MappingProxyType
    from typing import ClassVar

    class MyClass(DevBase):
        unique_vals = frozenset({"a", "b", "c"}) # frozenset is built in
        sequence = (1, 2, 3) # Use a tuple
        mapping = MappingProxyType({"a": 1, "b": 2}) # Read-only view on a dict
        class_storage: ClassVar = obj() # This is a class variable, not an instance default


.. _type-hints:

Type hints
----------

Type hints are not encouraged, but are allowed in new code with the
following stipulations:

#. They appear only in function argument lists and not in code bodies.

#. If a specific type (or types) is actually required by a type-hinted argument, then explicit runtime checking is required (e.g. via ``isinstance``).

#. If the actually allowed types are significantly broader or narrower than a type hint suggests, then documentation for that argument must appear in a docstring.

In Nexus, the meaning of type hints is explicitly as loose documentation. In and of themselves, they do not represent strict requirements.

Given the limited nature of documentation that type hints provide, docstrings
are generally preferred. This follows the philosophy:

* "If you need documentation, write it."

* "If you need specific types, enforce them."


.. _type-enforcement:

Type enforcement
----------------
Python is a dynamically-typed language, so do not insist on strong typing unless it is necessary. Instead, promote `duck-typing <https://docs.python.org/3/glossary.html#term-duck-typing>`__ where possible:

:lred:`No`

.. code-block:: python

    def exponentiate(values):
        assert isinstance(values,list) # A tuple, set, frozenset, or numpy array is also valid
        return [np.exp(v) for v in values]

:lgreen:`Yes`

.. code-block:: python

    def exponentiate(values):
        return [np.exp(v) for v in values]


Use type conversion for flexibility where possible:

:yellow:`Less`: ``assert isinstance(item_counts,np.ndarray)``

:green:`More`: ``item_counts = np.asarray(item_counts)``

:yellow:`Less`: ``assert isinstance(x,float)``

:green:`More`: ``x = float(x)``


.. _accessors:

Accessors
---------

Prefer dot style over string literals:

:red:`No`: ``name = data['nested']['path']['to']['name']``

:green:`Yes`: ``name = data.nested.path.to.name``

.. note::

    Dot style places limits on using :py:class:`dict` ; use :py:class:`~.obj` or :py:class:`~.dotdict` for those cases.
    If you are creating a user-facing function that accepts a mapping, do not require :py:class:`~.obj` or :py:class:`~.dotdict`,
    either use the :py:class:`dict` interface or convert the argument internally.

In classes, avoid set/get syntax. Instead promote dot access:

:lred:`No`

.. code-block:: python

    class MyClass(DevBase):
        def __init__(self, data):
            self._data = data
        #end def __init__

        def get_data(self):
            return self._data
        #end def get_data

        def set_data(self, data):
            self._data = data
        #end def set_data
    #end class MyClass

:lgreen:`Yes`

.. code-block:: python

    class MyClass(DevBase):
        def __init__(self,data):
            self.data = data
        #end def __init__
    #end class MyClass

If function execution is needed for set/get, use properties:

:lgreen:`Yes`

.. code-block:: python

    class MyClass(DevBase):
        @property
        def data(self):
            output = repackage_data(self.internal_data)
            return output
        #end def data

        @data.setter
        def data(self, data):
            self.internal_data = process_data(data)
        #end def data
    #end class MyClass


.. _imports:

Imports
-------

Use deferred (a.k.a. "lazy") imports.

With deferred imports, Nexus requires the minimum from external libraries and only fails if a specific missing functionality is actually used during execution.
This way, Nexus will often operate successfully even if some libraries are missing or outdated.

Function/method imports:

.. code-block:: python

    def save_data():
        import h5py
        ...
    #end def save_data


.. _exceptions:

Exceptions
----------

Use ``raise`` and select an `exception <https://docs.python.org/3/library/exceptions.html>`__ that best describes the problem.

At the user interfaces (i.e. ``generate_*``) liberally check input and use one of Nexus's errors, such as :py:exc:`~.NexusError`.

Make limited use of ``try-except``.
Only use if you know exactly what the error is and how to fix it and generally only covering small/localized sections of code (e.g. not on top of deeply nested function call trees)

Allow the code to error out to reveal the actual problem (``try-except`` can mask it).


.. _auto-formatting:

Auto-formatting
---------------

Nexus code is not required to follow formats used by auto-formatting tools.

If you use such a tool for newly committed code, ensure that the formatted code is largely consistent with existing code in the Nexus repository.


.. _unit-tests:

Unit tests
----------

Developers are encouraged to include unit tests with newly committed code.

Nexus uses the `pytest <https://docs.pytest.org/en/stable/index.html>`__ framework for all testing.

It is best if new or modified functionality is tested at some level of encapsulation.

If only one test is to be written, please enclose as much functionality as possible (i.e. at a higher level in the call tree).
As an example, see ``def test_ndgrid():`` in ``nexus/tests/test_numerics.py``.

More specific unit testing is encouraged if possible. See ``nexus/tests/test_periodic_table.py`` or the tests for :py:class:`~.PseudoSet` in ``nexus/tests/test_pseudopotential.py`` for examples of good unit tests.

If you have questions about how to integrate unit tests into Nexus, please reach out to a Nexus developer, or `open an issue on GitHub <https://github.com/QMCPACK/qmcpack/issues>`__.


.. _user_documentation:

User documentation
------------------

Developers are strongly encouraged to include documentation with their code, especially if that code will be accessed at the user level.

For example, all ``generate_*`` functions should have user documentation, but sadly do not.

If functionality includes operations commonly used by materials/chemical scientists, e.g. structure manipulation, adding documentation to the user guide is also encouraged.

See :py:mod:`~nexus.periodic_table` or :py:class:`~.PseudoSet` for examples of thorough documentation.

.. _documentation-style:

Documentation style
^^^^^^^^^^^^^^^^^^^

Class and function level docstrings are encouraged.  For these,
follow NumPy's formatting for class/function docs (see `the numpydoc format <https://numpydoc.readthedocs.io/en/latest/format.html>`__):

.. code-block:: python

    def real(value):
        """Return the real part of the complex argument.

        Parameters
        ----------
        value: array_like
            Input array.

        Returns
        -------
        output: ndarray or scalar
            The real component of the complex argument.
            If value is real, the type of value is used for
            the output. If value has complex elements,
            the returned type is float.

        Examples
        --------
        >>> real(1.3 + 2.8j)
        1.03

        It works with Numpy arrays too.

        >>> real(np.array([1+2j, 3+4j]))
        array([1., 3.])
        """
        return value.real
    #end def real


.. _uv_lock:

The ``uv.lock`` file
--------------------

A recent addition to Nexus is a lockfile generated by ``uv``. This is an automatically generated file that provides a set of exact versions for anyone looking to use a reproducible Python environment, most commonly for testing and development.

To sync your ``uv`` virtual environment to this file, you can simply run

.. code-block:: console

    > uv sync

To get all of the dependencies, you can run

.. code-block:: console

    > uv sync --all-groups --all-extras

which will install all of the optional dependencies for Nexus.

Updating the lockfile
^^^^^^^^^^^^^^^^^^^^^

Updating the ``uv.lock`` file should only be done when there has been a bump in the Nexus version or if a dependency has been added/removed from Nexus.

To update the ``uv.lock`` file, run

.. code-block:: console

    > uv lock --upgrade

which will automatically determine the correct versions of each package for each version of Python that is supported.


.. _dev_utils:

Developer utilities
-------------------

Nexus allows the use of specialized Python scripts for auto-generating code that goes into the distribution form of Nexus.
These are useful when a module, e.g. :py:mod:`~nexus.pwscf_input_defs` contains information that is *only* ever gathered from an external source.
Another example is :py:mod:`~nexus.periodic_table`, which contains an enumeration with data from IUPAC and NIST.

Here we have detailed descriptions of how to use the developer utilities currently available in Nexus.


.. _pwscf_input_autogen:

``gen_pwscf_input_data.py``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This script is used to automatically generate the :py:mod:`~nexus.pwscf_input_defs` module from Quantum ESPRESSO's ``INPUT_PW.def`` files.

This script requires that you have cloned the QE source code

.. code-block:: console

    > git clone https://gitlab.com/QEF/q-e.git

Additionally, you will need to download the dependencies listed in the `README.helpdoc file <https://gitlab.com/QEF/q-e/-/blob/develop/dev-tools/README.helpdoc>`__, as well as the ``xmltodict`` and ``packaging`` Python packages.
If you are using ``uv``, the Python packages will be installed automatically when running the script, otherwise you will need to run the following command

.. code-block:: console

    > python3 -m pip install xmltodict packaging

To use this script, go to `Quantum ESPRESSO's GitLab <https://gitlab.com/QEF/q-e/-/blob/develop/PW/Doc/INPUT_PW.def>`__ and download the ``INPUT_PW.def`` file there.
Rename the file to have the version attached, for example: ``INPUT_PW_7.6.0.def``.
Do this for each version of QE that you wish to include, using the version switcher and searching in the tags for the version.
Links for versions 4.0.4 to 7.6.0 are included in :py:mod:`~dev_utils.gen_pwscf_input_data` as comments after each line in the ``QE_DOC_VERSION_DATES`` dictionary.

Once you have all of the files, run the QE helpdoc tool to convert the files to XML format.

.. code-block:: console

    > ls
    INPUT_PW_7.5.0.def  INPUT_PW_7.6.0.def
    > /path/to/q-e/dev-tools/helpdoc INPUT_PW*.def

This will create ``.txt``, ``.xml``, and potentially ``.html`` files out of the ``.def`` files.

You may get the following output if the script fails:

.. code-block::

    ***
    *** Parsing the helpdoc.schema
    ***
    
    can't read "basedir": no such variable
        while executing
    "file join $basedir helpdoc.schema"
        (in namespace eval "::helpdoc::schema" script line 1)
        invoked from within
    "namespace eval schema { ::source [file join $basedir helpdoc.schema] }"
        (procedure "readSchema" line 3)
        invoked from within
    "readSchema "
        (procedure "::helpdoc::process" line 8)
        invoked from within
    "::helpdoc::process $argv"
        (file "../../q-e/dev-tools/helpdoc" line 51)

If you get this error, simply go to ``q-e/dev-tools/helpdoc.d/helpdoc.tcl`` and replace the following code (likely around lines 180-190)

.. code-block:: tcl

    namespace eval schema { ::source [file join $basedir helpdoc.schema] }

with

.. code-block:: tcl

    namespace eval schema { ::source [file join /path/to/q-e/dev-tools helpdoc.schema] }

and rerun the command.
A successful run should produce a large amount of output text, along the lines of

.. code-block::

    ***
    *** Parsing the helpdoc.schema
    ***

       parsing ROOTELEMENT input_description ... 
          parsing ATTRIBUTE distribution ... ok
          parsing ATTRIBUTE package ... ok
          parsing ATTRIBUTE program ... ok
          parsing OPTIONAL  ... 
             parsing INTERLEAVE  ... 
                parsing ELEMENT intro ... 
                OK - parsing ELEMENT intro completed
                parsing ELEMENT toc ... 
                OK - parsing ELEMENT toc completed
             OK - parsing INTERLEAVE completed
          OK - parsing OPTIONAL completed
          parsing +  ... 
    ...

Even if you get an error at the end like

.. code-block::

    [Error]
    Execution of xsltproc failed with error message:
    
    warning: failed to load external entity "input_xx.xsl"
    error
    xsltParseStylesheetFile : cannot parse input_xx.xsl
    compilation error: file INPUT_PW_7.6.0.xml line 5 element input_description
    xsltParseStylesheetProcess : document is not a stylesheet
       Executing:  /usr/bin/xsltproc --stringparam version "" --stringparam current-date "Tue Aug 11 13:43:20 EDT 2026" INPUT_PW_7.6.0.xml > INPUT_PW_7.6.0.html ...

the ``.xml`` files should still exist, and you can proceed. The error only relates to the production of the ``.html`` files, which are not needed for this script.

At this point, assuming you had version 7.5.0 and 7.6.0 in the directory (as shown above), your directory will now look like so

.. code-block:: console

    > ls
    INPUT_PW_7.5.0.def   INPUT_PW_7.5.0.txt  INPUT_PW_7.6.0.def   INPUT_PW_7.6.0.txt
    INPUT_PW_7.5.0.html  INPUT_PW_7.5.0.xml  INPUT_PW_7.6.0.html  INPUT_PW_7.6.0.xml

You can then call the script on the directory that contains all of the files (it will automatically grab only the ``.xml`` files), like so

.. code-block:: console

    > ./gen_pwscf_input_data.py /path/to/xml_dir/

    Parsing namelist data...

    Now parsing version 7.5.0...
      Parsing namelist CONTROL...
      Parsing namelist SYSTEM...
      Parsing namelist ELECTRONS...
      Parsing namelist IONS...
      Parsing namelist CELL...
      Parsing namelist FCP...
      Parsing namelist RISM...
    Success!

    Now parsing version 7.6.0...
      Parsing namelist CONTROL...
      Parsing namelist SYSTEM...
      Parsing namelist ELECTRONS...
      Parsing namelist IONS...
      Parsing namelist CELL...
      Parsing namelist FCP...
      Parsing namelist RISM...
    Success!

    Done parsing namelist data!

    Writing output to /your/path/to/qmcpack/nexus/nexus/pwscf_input_defs.py

At this point, the namelists have been successfully written.
You are encouraged to run ``pytest`` to verify that there are no errors caused by the script.
