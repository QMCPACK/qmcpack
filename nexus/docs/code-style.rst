.. _code-style:


Code Style for Developers
=========================

This section outlines requirements and recommendations regarding new code developed in Nexus.  Please consult and follow this documentation prior to submitting pull requests on GitHub (https://github.com/QMCPACK/qmcpack).  This will save a lot of time and effort during code reviews.


.. _text-style:

Text style
----------
  
Class names are camel case: ``class BaseClass:``.


All other code is snake case: ``def perform_work():``.

Camel case means the first letter of each word is capitalized with no separation between words. Snake case means each word is lowercase with single underscores separating each word.



.. _variable-names:

Variable and function names
---------------------------

Be descriptive and generally limit names to three words or less for legibility.

No: ``def slice_dice_and_mince_data(data):``

Yes: ``def slice_data(data):``

No: ``total_amount_of_money_in_bag = total_number_of_pennies_in_bag + number_of_pennies_in_a_nickel*total_number_of_nickels_in_bag + number_of_pennies_in_a_dime*total_number_of_nickels_in_bag + ...``

Yes: ``money_tot = pennies + 5*nickels + 10*dimes + ...``



.. _documentation-style:

Documentation style
-------------------

Use inline comments, but not on every line and be brief.

Class and function level docstrings are encouraged.  For these, 
follow ``numpy`` formatting for class/function docs (see https://numpydoc.readthedocs.io/en/latest/format.html):


::

  def real(value):
      """
      Return the real part of the complex argument.
      
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
      """
      ...
  #end def real



.. _encapsulation:

Encapsulation
-------------

Function and class endings are delimited with a comment:

::

    def perform_work():
        ...
    #end def perform_work

    BaseClass:
        ...
    #end BaseClass

Enclosing comments such as ``#end if`` and ``#end for`` are not required,
but are encouraged in cases of deep nesting for readability.

For closing parentheses (and similar), the closing character should have the same indentation as the function arguments or e.g. list elements.

No:

::

    def operate(key1 = 'key1',
                key2 = 'key2',
                key3 = 'key3',
    ):            
        key = (key1,key2,key3)
        return key
    #end def operate

Yes:

::

    def operate(key1 = 'key1',
                key2 = 'key2',
                key3 = 'key3',
                ):            
        key = (key1,key2,key3)
        return key
    #end def operate

No:
 
::

    d = dict(a = 1,
             b = 2,
    )

Yes:
 
::

    d = dict(a = 1,
             b = 2,
             )

    # or

    d = dict(a = 1,
             b = 2)

No:
 
::

    l = ['key1',
         'key2',
         'key3',
    ]


Yes:
 
::

    l = ['key1',
         'key2',
         'key3',
         ]

    # or

    l = ['key1',
         'key2',
         'key3']


        


.. _function-arguments:

Function arguments
------------------

Use keyword arguments when calling a function, unless they are few.

Do not require positional or keyword arguments in function definitions.

Functions with arbitrary arguments are allowed: ``def mega_function(*args,**kwargs):``




.. _classes:

Classes
-------

Avoid large/complex classes (yes, this is ironic given some of the early Nexus code).

Use tight procedural programming where possible.

Use "heavy" base classes to promote light derived classes.

Avoid set/get accessors.



.. _type-hints:

Type hints
----------

Type hints are not encouraged, but are allowed in new code with the 
following stipulations:

  1. They appear only in function argument lists and not in code bodies.

  2. If a specific type (or types) is actually required by a type-hinted argument, then explicit runtime checking is required (e.g. via ``isinstance``).

  3. If the actually allowed types are significantly broader or narrower than a type hint suggests, then documentation for that argument must appear in a docstring. 

In Nexus, the meaning of type hints is explicitly as loose documentation. In and of themselves, they do not represent strict requirements.

Given the limited nature of documentation that type hints provide, docstrings 
are generally preferred. This follows the philosophy:

* "If you need documentation, write it."

* "If you need specific types, enforce them."



.. _type-enforcement:

Type enforcement
----------------
Python is a dynamic language, so do not insist on strong typing unless it is necessary.  Instead, promote duck-typing where possible:

No:

::

    def exponentiate(values):
        assert isinstance(values,list)
        return [np.exp(v) for v in values]

Yes:

::

    def exponentiate(values):
        return [np.exp(v) for v in values]


Use type conversion for flexibility where possible:

Less: ``assert isinstance(item_counts,np.ndarray)``

More: ``item_counts = np.asarray(item_counts)``

Less: ``assert isinstance(x,float)``

More: ``x = float(x)``


.. _accessors:

Accessors
---------

Strictly use dot style over string literals:
    
No : ``name = data['nested']['path']['to']['name']``

Yes: ``name = data.nested.path.to.name``

Note: dot style places limits on using ``dict`` ; use ``obj`` for those cases.

In classes, avoid set/get syntax. Instead promote dot access:

No:

::

    class DataClass:
        def __init__(self,data):
            self.data = data
        #end def __init__
      
        def get_data(self):
            return self.data
        #end def get_data

        def set_data(self,data):
            self.data = data
        #end def set_data
    #end class DataClass

Yes: 

::

    DataClass:
        def __init__(self,data):
            self.data = data
        #end def __init__
    #end class DataClass

If function execution is needed for set/get, use properties:

Yes:

::

    class DataClass:
        @property
        def data(self,data):
            ...
        #end def data

        @data.setter
        def data(self,data):
            ...
        #end def data
    #end class DataClass



.. _imports:

Imports
-------

Use deferred (a.k.a. "lazy") imports.  

With deferred imports, Nexus requires the minimum from external libraries and only fails if a specific missing functionality is actually used during execution.  This way, Nexus will often operate successfully even if some libraries are missing or outdated.



Top-level/header imports:

::

    try:
        import h5py
    except:
        h5py = unavailable('h5py')
  
Function/method imports:

::

    def save_data():
        import h5py
        ...
    #end def save_data



.. _exceptions:

Exceptions
----------

Use ``raise`` only for low-level utility code.

Use ``error()`` for code at higher levels, especially close to user interfaces.

At the user interfaces (i.e. ``generate_*``) liberally check input and use ``error()``.

Make limited use of ``try-except``.  Only use if you know exactly what the error is and how to fix it.

Allow the code to error out to reveal the actual problem (``try-except`` masks it).



.. _auto-formatting:

Auto-formatting
---------------

Nexus code is not required to follow formats used by auto-formatting tools.

If you use such a tool for newly committed code, ensure that the formatted code is largely consistent with existing code in the Nexus repository.



.. _unit-tests:

Unit tests
----------

Developers are encouraged to include unit tests with newly committed code.

It is best if new or modified functionality is tested at some level of encapsulation.

If only one test is to be written, please enclose as much functionality as possible (i.e. at a higher level in the call tree). 

As an example, see ``def test_ndgrid():`` in ``nexus/tests/unit/test_numerics.py``.  

If you have questions about how to integrate unit tests into Nexus, please reach out to a Nexus developer. 



.. _user_documentation:

User documentation
------------------

Developers are encouraged to include documentation in the Nexus manual for code that will be accessed at the user level.

For example, all ``generate_*`` functions should have user documentation, but sadly do not.

If functionality includes operations commonly used by materials/chemical scientists, e.g. structure manipulation, adding documentation to the manual is also encouraged. 


.. _uv_lock:

The ``uv.lock`` file
--------------------

A recent addition to Nexus is a lockfile generated by ``uv``. This is an automatically generated file that provides a set of exact versions for anyone looking to use a reproducible Python environment, most commonly for testing and development.

To sync your ``uv`` virtual environment to this file, you can simply run

.. code-block::

    uv sync

To get all of the dependencies, you can run

.. code-block::

    uv sync --all-groups --all-extras

which will install all of the optional dependencies for Nexus.

Updating the lockfile
^^^^^^^^^^^^^^^^^^^^^

Updating the ``uv.lock`` file should only be done when there has been a bump in the Nexus version or if a dependency has been added/removed from Nexus.

To update the ``uv.lock`` file, run

.. code-block::

    uv lock --upgrade

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
