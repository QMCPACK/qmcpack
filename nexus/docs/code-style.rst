.. _code-style:


Code Style for Developers
=========================

This section outlines requirements and recommendations regarding new code developed in Nexus.  Please consult and follow this documentation prior to submitting pull requests on GitHub (https://github.com/QMCPACK/qmcpack).  This will save a lot of time and effort during code reviews.


.. _text-style:

Text style
----------
  
Class names are camel case: ``class BaseClass``.


All other code is snake case: ``def perform_work():``.



.. _documentation-style

Documentation style
-------------------

Use inline comments, but not on every line and be brief.

Follow ``numpy`` formatting for class/function docs:


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



.. _encapsulation

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



.. _function-arguments

Function arguments
------------------

Use keyword arguments when calling a function, unless they are few.

Do not require positional or keyword arguments in function definitions.

Functions with arbitrary arguments are allowed: ``def mega_function(*args,**kwargs):``




.. _classes

Classes
-------

Avoid large/complex classes (yes, this is ironic given some of the early Nexus code).

Use tight procedural programming where possible.

Use "heavy" base classes to promote light derived classes.

Avoid set/get accessors.



.. _type-hints

Type hints
----------

No type hints, except when there is no other option (e.g. ``dataclass``).

If you need documentation, write it.

If you need specific types, enforce them.



.. _type-enforcement

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



.. _accessors

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
        #end def set
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



.. _imports

Imports
-------

Use lazy imports.

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



.. _exceptions

Exceptions
----------

Use ``raise`` only for low-level utility code.

Use ``error()`` for code at higher levels, especially close to user interfaces.

At the user interfaces (i.e. ``generate_*``) liberally check input and use ``error()``.

Make limited use of ``try-except``.  Only use if you know exactly what the error is and how to fix it.

Allow the code to error out to reveal the actual problem (``try-except`` masks it).



.. _auto-formatting

Auto-formatting
---------------

No auto-formatting (e.g. Black).



.. _static-analysis

Static analysis
---------------

Do not strictly enforce/require static analysis. Instead, prefer unit testing.



.. _other-generalities

Other generalities
------------------
PEP8 is not the source of all truth.
