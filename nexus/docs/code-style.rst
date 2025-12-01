.. _code-style:


Code Style for Developers
=========================

This section outlines requirements and recommendations regarding new code developed in Nexus.  Please consult and follow this documentation prior to submitting pull requests on GitHub (https://github.com/QMCPACK/qmcpack).  This will save a lot of time and effort during code reviews.


.. _text-style:

Text style
----------
  
Class names are camel case: ``class BaseClass:``.


All other code is snake case: ``def perform_work():``.



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
follow ``numpy`` formatting for class/function docs:


::

  def real(val):
      '''
      Return the real part of the complex argument.
      
      Parameters:
      
          val: array_like
      
             Input array.
      
      Returns:
      
          out: ndarray or scalar
      
             The real component of the complex argument. 
             If val is real, the type of val is used for 
             the output. If val has complex elements, 
             the returned type is float.
      '''
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


