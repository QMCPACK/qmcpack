.. _developguide:

Development Guide
=================

The section gives guidance on how to extend the functionality of QMCPACK. Future examples will likely include topics such as the
addition of a Jastrow function or a new QMC method.

QMCPACK coding standards
------------------------

This chapter presents what we collectively have agreed are best practices for the code. This includes formatting style, naming
conventions, documentation conventions, and certain prescriptions for C++ language use. At the moment only the formatting can be
enforced in an objective fashion.

New development should follow these guidelines, and contributors are expected to adhere to them as they represent an integral part
of our effort to continue QMCPACK as a world-class, sustainable QMC code. Although some of the source code has a ways to go to
live up to these ideas, new code, even in old files, should follow the new conventions not the local conventions of the file
whenever possible. Work on the code with continuous improvement in mind rather than a commitment to stasis.

The `current workflow conventions`_ for the project are described in the wiki on the GitHub repository. It will save you and all
the maintainers considerable time if you read these and ask questions up front.

A PR should follow these standards before inclusion in the mainline. You can be sure of properly following the formatting
conventions if you use clang-format.  The mechanics of clang-format setup and use can be found at
https://github.com/QMCPACK/qmcpack/wiki/Source-formatting.

The clang-format file found at ``qmcpack/src/.clang-format`` should be run over all code touched in a PR before a pull request is
prepared. We also encourage developers to run clang-tidy with the ``qmcpack/src/.clang-tidy`` configuration over all new code.

As much as possible, try to break up refactoring, reformatting, feature, and bugs into separate, small PRs. Aim for something that
would take a reviewer no more than an hour. In this way we can maintain a good collective development velocity.

.. _current workflow conventions: https://github.com/QMCPACK/qmcpack/wiki/Development-workflow

Files
-----

Each file should start with the header.

::

  //////////////////////////////////////////////////////////////////////////////////////
  // This file is distributed under the University of Illinois/NCSA Open Source License.
  // See LICENSE file in top directory for details.
  //
  // Copyright (c) 2021 QMCPACK developers
  //
  // File developed by: Name, email, affiliation
  //
  // File created by: Name, email, affiliation
  //////////////////////////////////////////////////////////////////////////////////////

If you make significant changes to an existing file, add yourself to the list of "developed by" authors.

File organization
~~~~~~~~~~~~~~~~~

Header files should be placed in the same directory as their implementations. Unit tests should be written for all new
functionality. These tests should be placed in a ``tests`` subdirectory below the implementations.

File names
~~~~~~~~~~

Each class should be defined in a separate file with the same name as the class name. Use separate ``.cpp`` implementation files
whenever possible to aid in incremental compilation.

The filenames of tests are composed by the filename of the object tested and the prefix ``test_``. The filenames of *fake* and
*mock* objects used in tests are composed by the prefixes ``fake_`` and ``mock_``, respectively, and the filename of the object
that is imitated.

Header files
~~~~~~~~~~~~

All header files should be self-contained (i.e., not dependent on following any other header when it is included). Nor should they
include files that are not necessary for their use (i.e., headers needed only by the implementation). Implementation files should
not include files only for the benefit of files they include.

There are many header files that currently violate this. Each header must use ``#define`` guards to prevent multiple inclusion.
The symbol name of the ``#define`` guards should be ``NAMESPACE(s)_CLASSNAME_H``.

Includes
~~~~~~~~

Related header files should be included without any path. Header files from external projects and standard libraries should be
includes using the ``<iostream>`` convention, while headers that are part of the QMCPACK project should be included using the
``"our_header.h"`` convention.

We are now using a new header file inclusion style following the modern CMake transition in QMCPACK, while the legacy code may
still use the legacy style. Newly written code and refactored code should be transitioned to the new style.

New style for modern CMake
^^^^^^^^^^^^^^^^^^^^^^^^^^

In QMCPACK, include paths are handled by modern CMake target dependency. Every top level folder is at least one target. For
example, ``src/Particle/CMakeLists.txt`` defines `qmcparticle` target. It propagates include path ``qmcpack/src/Particle`` to
compiling command lines in CMake via

::

  TARGET_INCLUDE_DIRECTORIES(qmcparticle PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

For this reason, the file ``qmcpack/src/Particle/Lattice/ParticleBConds3DSoa.h`` should be included as

::

  #include "Lattice/ParticleBConds3DSoa.h"

If the compiled file is not part of the same target as `qmcparticle`, the target it belongs to should have a dependency on
`qmcparticle`. For example, test source files under ``qmcpack/src/Particle/tests`` are not part of `qmcparticle` and thus requires
the following additional CMake setting

::

  TARGET_LINK_LIBRARIES(${UTEST_EXE} qmcparticle)

Legacy style
^^^^^^^^^^^^

Header files should be included with the full path based on the ``src`` directory. For example, the file
``qmcpack/src/QMCWaveFunctions/SPOSet.h`` should be included as

::

  #include "QMCWaveFunctions/SPOSet.h"

Even if the included file is located in the same directory as the including file, this rule should be obeyed.

Ordering
^^^^^^^^

For readability, we suggest using the following standard order of includes:

#. related header

#. std C library headers

#. std C++ library headers

#. Other libraries’ headers

#. QMCPACK headers

In each section the included files should be sorted in alphabetical order.

Naming
------

The balance between description and ease of implementation should be balanced such that the code remains self-documenting within a
single terminal window.  If an extremely short variable name is used, its scope must be shorter than :math:`\sim 40` lines. An
exception is made for template parameters, which must be in all CAPS. Legacy code contains a great variety of hard to read code
style, read this section and do not imitate existing code that violates it.

Namespace names
~~~~~~~~~~~~~~~

Namespace names should be one word, lowercase.

Type and class names
~~~~~~~~~~~~~~~~~~~~

Type and class names should start with a capital letter and have a capital letter for each new word. Underscores (``_``) are not
allowed. It's redundant to end these names with ``Type`` or ``_t``.

::
   \\no
   using ValueMatrix_t = Matrix<Value>;
   using RealType = double;

Variable names
~~~~~~~~~~~~~~

Variable names should not begin with a capital letter, which is reserved for type and class names. Underscores (``_``) should be
used to separate words.

Class data members
~~~~~~~~~~~~~~~~~~

Class private/protected data members names should follow the convention of variable names with a trailing underscore (``_``). The use of public member functions is discourage, rethink the need for it in the first place. Instead ``get`` and ``set`` functions are the preferred access method.

(Member) function names
~~~~~~~~~~~~~~~~~~~~~~~

Function names should start with a lowercase character and have a capital letter for each new word. The exception are the special cases for prefixed multiwalker (``mw_``) and flex (``flex_``) batched API functions. Coding convention should follow after those prefixes.

Template Parameters
~~~~~~~~~~~~~~~~~~~

Template parameters names should be in all caps with (``_``) separating words.  It's redundant to end these names with ``_TYPE``,

Lambda expressions
~~~~~~~~~~~~~~~~~~

Named lambda expressions follow the naming convention for functions:

::

  auto myWhatever = [](int i) { return i + 4; };

Macro names
~~~~~~~~~~~

Macro names should be all uppercase and can include underscores (``_``). The underscore is not allowed as first or last character.

Test case and test names
~~~~~~~~~~~~~~~~~~~~~~~~

Test code files should be named as follows:

::

  class DiracMatrix;
  //leads to
  test_dirac_matrix.cpp
  //which contains test cases named
  TEST_CASE("DiracMatrix_update_row","[wavefunction][fermion]")

where the test case covers the ``updateRow`` and  ``[wavefunction][fermion]`` indicates the test belongs to the fermion wavefunction functionality.

Comments
--------

Comment style
~~~~~~~~~~~~~

Use the ``// Comment`` syntax for actual comments.

Use

::

  /** base class for Single-particle orbital sets
   *
   * SPOSet stands for S(ingle)P(article)O(rbital)Set which contains
   * a number of single-particle orbitals with capabilities of
   * evaluating \f$ \psi_j({\bf r}_i)\f$
   */

or

::

  ///index in the builder list of sposets
  int builder_index;

Documentation
~~~~~~~~~~~~~

Doxygen will be used for source documentation. Doxygen commands should be used when appropriate guidance on this has been decided.

File docs
^^^^^^^^^

Do not put the file name after the ``\file`` Doxygen command. Doxygen will fill it in for the file the tag appears in.

::

  /** \file
   *  File level documentation
   */

Class docs
^^^^^^^^^^

Every class should have a short description (in the header of the file) of what it is and what is does. Comments for public class
member functions follow the same rules as general function comments. Comments for private members are allowed but are not
mandatory.

Function docs
^^^^^^^^^^^^^

For function parameters whose type is non-const reference or pointer to non-const memory, it should be specified if they are input
(In:), output (Out:) or input-output parameters (InOut:).

Example:

::

  /** Updates foo and computes bar using in_1 .. in_5.
   * \param[in] in_3
   * \param[in] in_5
   * \param[in,out] foo
   * \param[out] bar
   */

  //This is probably not what our clang-format would do
  void computeFooBar(Type in_1, const Type& in_2, Type& in_3,
                     const Type* in_4, Type* in_5, Type& foo,
                     Type& bar);

Variable documentation
^^^^^^^^^^^^^^^^^^^^^^

Name should be self-descriptive.  If you need documentation consider renaming first.

Golden rule of comments
~~~~~~~~~~~~~~~~~~~~~~~

If you modify a piece of code, also adapt the comments that belong to it if necessary.

Formatting and "style"
----------------------

Use the provided clang-format style in ``src/.clang-format`` to format ``.h``, ``.hpp``, ``.cu``, and ``.cpp`` files. Many of the following rules will be applied to the code by clang-format, which should allow you to ignore most of them if you always run it on your modified code.

You should use clang-format support and the ``.clangformat`` file with your editor, use a Git precommit hook to run clang-format
or run clang-format manually on every file you modify.  However, if you see numerous formatting updates outside of the code you
have modified, first commit the formatting changes in a separate PR.

Indentation
~~~~~~~~~~~

Indentation consists of two spaces. Do not use tabs in the code.

Line length
~~~~~~~~~~~

The length of each line of your code should be at most *120* characters.

Horizontal spacing
~~~~~~~~~~~~~~~~~~

No trailing white spaces should be added to any line. Use no space before a comma (``,``) and a semicolon (``;``), and add a space
after them if they are not at the end of a line.

Preprocessor directives
~~~~~~~~~~~~~~~~~~~~~~~

The preprocessor directives are not indented.
The hash is the first character of the line.

Binary operators
~~~~~~~~~~~~~~~~

The assignment operators should always have spaces around them.

Unary operators
~~~~~~~~~~~~~~~

Do not put any space between an unary operator and its argument.

Types
~~~~~

The ``using`` syntax is preferred to ``typedef`` for type aliases. If the actual type is not excessively long or complex, simply
use it; renaming simple types makes code less understandable.

Pointers and references
~~~~~~~~~~~~~~~~~~~~~~~

Pointer or reference operators should go with the type. But understand the compiler reads them from right to left.

::

  Type* var;
  Type& var;

  //Understand this is incompatible with multiple declarations
  Type* var1, var2; // var1 is a pointer to Type but var2 is a Type.

Templates
~~~~~~~~~

The angle brackets of templates should not have any external or internal padding.

::

  template<class C>
  class Class1;

  Class1<Class2<type1>> object;

Vertical spacing
~~~~~~~~~~~~~~~~

Use empty lines when it helps to improve the readability of the code, but do not use too many. Do not use empty lines after a
brace that opens a scope or before a brace that closes a scope. Each file should contain an empty line at the end of the file.
Some editors add an empty line automatically, some do not.

Variable declarations and definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Avoid declaring multiple variables in the same declaration, especially if they are not fundamental types:

  ::

    int x, y;                        // Not recommended
    Matrix a("my-matrix"), b(size);  // Not allowed

    // Preferred
    int x;
    int y;
    Matrix a("my-matrix");
    Matrix b(10);

- Use the following order for keywords and modifiers in  variable declarations:

  ::

    // General type
    [static] [const/constexpr] Type variable_name;

    // Pointer
    [static] [const] Type* [const] variable_name;

    // Integer
    // the int is not optional not all platforms support long, etc.
    [static] [const/constexpr] [signedness] [size] int variable_name;

    // Examples:
    static const Matrix a(10);


    const double* const d(3.14);
    constexpr unsigned long l(42);

Function declarations and definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The return type should be on the same line as the function name. Parameters should also be on the same line unless they do not fit
on it, in which case one parameter per line aligned with the first parameter should be used.

Also include the parameter names in the declaration of a function, that is,

::

  // calculates a*b+c
  double function(double a, double b, double c);

  // avoid
  double function(double, double, double);

  // dont do this
  double function(BigTemplatedSomething<double> a, BigTemplatedSomething<double> b,
                  BigTemplatedSomething<double> c);

  // do this
  double function(BigTemplatedSomething<double> a,
                  BigTemplatedSomething<double> b,
                  BigTemplatedSomething<double> c);

Conditionals
~~~~~~~~~~~~

Examples:

::

  if (condition)
    statement;
  else
    statement;

  if (condition)
  {
    statement;
  }
  else if (condition2)
  {
    statement;
  }
  else
  {
    statement;
  }

Switch statement
~~~~~~~~~~~~~~~~

Switch statements should always have a default case.

Example:

::

  switch (var)
  {
    case 0:
      statement1;
      statement2;
      break;

    case 1:
      statement1;
      statement2;
      break;

    default:
      statement1;
      statement2;
  }

Loops
~~~~~

Examples:

::

  for (statement; condition; statement)
    statement;

  for (statement; condition; statement)
  {
    statement1;
    statement2;
  }

  while (condition)
    statement;

  while (condition)
  {
    statement1;
    statement2;
  }

  do
  {
    statement;
  }
  while (condition);

.. _class-format:

Class format
~~~~~~~~~~~~

``public``, ``protected``, and ``private`` keywords are not indented.

Example:

::

  class Foo : public Bar
  {
  public:
    Foo();
    explicit Foo(int var);

    void function();
    void emptyFunction() {}

    void setVar(const int var)
    {
      var_ = var;
    }
    int getVar() const
    {
      return var_;
    }

  private:
    bool privateFunction();

    int var_;
    int var2_;
  };

Constructor initializer lists
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Examples:

::

  // When everything fits on one line:
  Foo::Foo(int var) : var_(var)
  {
    statement;
  }

  // If the signature and the initializer list do not
  // fit on one line, the colon is indented by 4 spaces:
  Foo::Foo(int var)
      : var_(var), var2_(var + 1)
  {
    statement;
  }

  // If the initializer list occupies more lines,
  // they are aligned in the following way:
  Foo::Foo(int var)
      : some_var_(var),
        some_other_var_(var + 1)
  {
    statement;
  }

  // No statements:
  Foo::Foo(int var)
      : some_var_(var) {}

Namespace formatting
~~~~~~~~~~~~~~~~~~~~

The content of namespaces is not indented. A comment should indicate when a namespace is closed. (clang-format will add these if
absent). If nested namespaces are used, a comment with the full namespace is required after opening a set of namespaces or an
inner namespace.

Examples:

::

  namespace ns
  {
  void foo();
  }  // ns

::

  namespace ns1
  {
  namespace ns2
  {
  // ns1::ns2::
  void foo();

  namespace ns3
  {
  // ns1::ns2::ns3::
  void bar();
  }  // ns3
  }  // ns2

  namespace ns4
  {
  namespace ns5
  {
  // ns1::ns4::ns5::
  void foo();
  }  // ns5
  }  // ns4
  }  // ns1

QMCPACK C++ guidance
--------------------

The guidance here, like any advice on how to program, should not be treated as a set of rules but rather the hard-won wisdom of
many hours of suffering development. In the past, many rules were ignored, and the absolute worst results of that will affect
whatever code you need to work with. Your PR should go much smoother if you do not ignore them.

Encapsulation
~~~~~~~~~~~~~

A class is not just a naming scheme for a set of variables and functions. It should provide a logical set of methods, could
contain the state of a logical object, and might allow access to object data through a well-defined interface related variables,
while preserving maximally ability to change internal implementation of the class.

Do not use ``struct`` as a way to avoid controlling access to the class. Only in rare cases where a class is a fully public data
structure ``struct`` is this appropriate. Ignore (or fix one) the many examples of this in QMCPACK.

Do not use inheritance primarily as a means to break encapsulation. If your class could aggregate or compose another class, do
that, and access it solely through its public interface. This will reduce dependencies.

Casting
~~~~~~~

In C++ source, avoid C style casts; they are difficult to search for and imprecise in function. An exception is made for
controlling implicit conversion of simple numerical types.

Explicit C++ style casts make it clear what the safety of the cast is and what sort of conversion is expected to be possible.

::

  int c = 2;
  int d = 3;
  double a;
  a = (double)c / d;  // Ok

  const class1 c1;
  class2* c2;
  c2 = (class2*)&c1; // NO
  SPOSetAdvanced* spo_advanced = new SPOSetAdvanced();

  SPOSet* spo = (SPOSet*)spo_advanced; // NO
  SPOSet* spo = static_cast<SPOSet*>(spo_advanced); // OK if upcast, dangerous if downcast

Pre-increment and pre-decrement
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the pre-increment (pre-decrement) operator when a variable is incremented (decremented) and the value of the expression is not
used. In particular, use the pre-increment (pre-decrement) operator for loop counters where i is not used:

::

  for (int i = 0; i < N; ++i)
  {
    doSomething();
  }

  for (int i = 0; i < N; i++)
  {
    doSomething(i);
  }

The post-increment and post-decrement operators create an unnecessary copy that the compiler cannot optimize away in the case of
iterators or other classes with overloaded increment and decrement operators.

Alternative operator representations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Alternative representations of operators and other tokens such as ``and``, ``or``, and ``not`` instead of ``&&``, ``||``, and ``!`` are not allowed.
For the reason of consistency, the far more common primary tokens should always be used.

Use of const
~~~~~~~~~~~~

- Add the ``const`` qualifier to all function parameters that are not modified in the function body.

- For parameters passed by value, add only the keyword in the function definition.

- Member functions should be specified const whenever possible.

  ::

    // Declaration
    int computeFoo(int bar, const Matrix& m)

    // Definition
    int computeFoo(const int bar, const Matrix& m)
    {
      int foo = 42;

      // Compute foo without changing bar or m.
      // ...

      return foo;
    }

    class MyClass
    {
      int count_
      ...
      int getCount() const { return count_;}
    }

Smart pointers
~~~~~~~~~~~~~~

Use of smart pointers is being adopted to help make QMCPACK memory leak free. Prior to C++11, C++ uses C-style pointers. A C-style
pointer can have several meanings and the ownership of a piece of help memory may not be clear. This leads to confusion and causes
memory leaks if pointers are not managed properly. Since C++11, smart pointers were introduced to resolve this issue. In addition,
it demands developers to think about the ownership and lifetime of declared pointer objects.

std::unique_ptr
^^^^^^^^^^^^^^^

A unique pointer is the unique owner of a piece of allocated memory. Pointers in per-walker data structure with distinct contents
should be unique pointers. For example, every walker has a trial wavefunction object which contains an SPO object pointer. Because
the SPO object has a vector to store SPO evaluation results, it cannot be shared between two trial wavefunction objects. For this
reason the SPO object pointer should be an unique pointer.

In QMCPACK, most raw pointers can be directly replaced with ``std::unique_ptr``.
Corresponding use of ``new`` operator can be replaced with ``std:make_unique``.

std::shared_ptr
^^^^^^^^^^^^^^^

A shared pointer is the shared owner of a piece of allocated memory. Moving a pointer ownership from one place to another should
not use shared pointers but C++ move semantics. Shared contents between walkers may be candidates for shared pointers. For example,
although the Jastrow factor object must be unique per walker, the pointer to the parameter data structure can be a shared pointer.
During Jastrow optimization, any update to the parameter data managed by the shared pointer will be effective immediately in all
the Jastrow objects. In another example, spline coefficients are managed by a shared pointer which achieves a single copy in
memory shared by an SPOSet and all of its clones.

.. _distance-tables:

Particles and distance tables
-----------------------------

ParticleSets
~~~~~~~~~~~~

The ``ParticleSet`` class stores particle positions and attributes
(charge, mass, etc).

The ``R`` member stores positions. For calculations, the ``R`` variable
needs to be transferred to the structure-of-arrays (SoA) storage in
``RSoA``. This is done by the ``update`` method. In the future the
interface may change to use functions to set and retrieve positions so
the SoA transformation of the particle data can happen automatically.
For now, it is crucial to call ``P.update()`` to populate ``RSoA`` anytime ``P.R`` is changed. Otherwise, the distance tables associated with ``R`` will be uninitialized or out-of-date.

::

  const SimulationCell sc;
  ParticleSet elec(sc), ions(sc);
  elec.setName("e");
  ions.setName("ion0");

  // initialize ions
  ions.create(2);
  ions.R[0] = {0.0, 0.0, 0.0};
  ions.R[1] = {0.5, 0.5, 0.5};
  ions.update(); // transfer to RSoA

  // initialize elec
  elec.create(2);
  elec.R[0] = {0.0, 0.0, 0.0};
  elec.R[1] = {0.0, 0.25, 0.0};
  const int itab = elec.addTable(ions);
  elec.update(); // update RSoA and distance tables

  // d_table is an electron-ion distance table
  const auto& d_table = elec.getDistTableAB(itab);

A particular distance table is retrieved with ``getDistTable``. Use
``addTable`` to add a ``ParticleSet`` and return the index of the
distance table. If the table already exists the index of the existing
table will be returned.

The mass and charge of each particle is stored in ``Mass`` and ``Z``.
The flag, ``SameMass``, indicates if all the particles have the same
mass (true for electrons).

Groups
^^^^^^

Particles can belong to different groups. For electrons, the groups are
up and down spins. For ions, the groups are the atomic element. The
group type for each particle can be accessed through the ``GroupID``
member. The number of groups is returned from ``groups()``. The total
number particles is accessed with ``getTotalNum()``. The number of
particles in a group is ``groupsize(int igroup)``.

The particle indices for each group are found with ``first(int igroup)``
and ``last(int igroup)``. These functions only work correctly if the
particles are packed according to group. The flag, ``IsGrouped``,
indicates if the particles are grouped or not. The particles will not be
grouped if the elements are not grouped together in the input file. This
ordering is usually the responsibility of the converters.

Code can be written to only handle the grouped case, but put an assert
or failure check if the particles are not grouped. Otherwise the code
will give wrong answers and it can be time-consuming to debug.

Distance tables
~~~~~~~~~~~~~~~

Distance tables store distances between particles. There are symmetric
(AA) tables for distance between like particles (electron-electron or
ion-ion) and asymmetric (AB) tables for distance between unlike
particles (electron-ion)

The ``Distances`` and ``Displacements`` members contain the data. The
indexing order is target index first, then source. For electron-ion
tables, the sources are the ions and the targets are the electrons.

Looping over particles
~~~~~~~~~~~~~~~~~~~~~~

Some sample code on how to loop over all the particles in an electron-ion distance table:

::


  // d_table is an electron-ion distance table

  for (int jat = 0; j < d_table.targets(); jat++) { // Loop over electrons
    for (int iat = 0; i < d_table.sources(); iat++) { // Loop over ions
       d_table.Distances[jat][iat];
    }
  }

Interactions sometimes depend on the type of group of the particles. The
code can loop over all particles and use ``GroupID[idx]`` to choose the
interaction. Alternately, the code can loop over the number of groups
and then loop from the first to last index for those groups. This method
can attain higher performance by effectively hoisting tests for group ID
out of the loop.

An example of the first approach is

::


  // P is a ParticleSet

  for (int iat = 0; iat < P.getTotalNum(); iat++) {
    int group_idx = P.GroupID[iat];
    // Code that depends on the group index
  }

An example of the second approach is

::

  // P is a ParticleSet
  assert(P.IsGrouped == true); // ensure particles are grouped

  for (int ig = 0; ig < P.groups(); ig++) { // loop over groups
    for (int iat = P.first(ig); iat < P.last(ig); iat++) { // loop over elements in each group
       // Code that depends on group
    }
  }

Wavefunction
------------

A full ``TrialWaveFunction`` is formulated as a product
of all the components. Each component derives from ``WaveFunctionComponent``.

.. math::
     \psi = \prod_c {\tilde \psi_c}

QMCPACK doesn't directly use the product form but mostly uses the log of the wavefunction.
It is a natural fit for QMC algorithms and offers a numerical advantage on computers.
The log value grows linearly instead of exponentially, beyond the range of double precision,
with respect to the electron counts in a Slater-Jastrow wave function.

The code contains an example of a
wavefunction component for a Helium atom using a simple form and is
described in :ref:`helium-wavefunction-example`

Mathematical preliminaries
~~~~~~~~~~~~~~~~~~~~~~~~~~

The wavefunction evaluation functions compute the log of the
wavefunction, the gradient and the Laplacian of the log of the
wavefunction. Expanded, the gradient and Laplacian are

.. math::
  :label: eq264

   \begin{aligned}
     {\bf G} & = \{ \nabla_i \ln(\psi) \} = \left\{ \sum_c \nabla_i \ln(\tilde \psi_c)\right\} , & {\bf \tilde G} & = \{ \nabla_i \ln(\tilde \psi) \} = \left\{ \frac{\nabla_i \tilde \psi}{\tilde \psi}\right\} \\
     {\bf L} & = \{ \nabla^2_i \ln(\psi) \} = \left\{ \sum_c \nabla^2_i \ln(\tilde \psi_c) \right\}, & {\bf \tilde L}  & = \{ \nabla^2_i \ln(\tilde \psi) \} = \left\{\frac{{\nabla^2_i} \tilde \psi}{\tilde \psi} - {\tilde G_i} \cdot {\tilde G_i} \right\}
   \end{aligned}

where :math:`i` is the electron index.
In this separable form, each wavefunction component computes its :math:`{\bf \tilde G}` ``WaveFunctionComponent::G`` and :math:`{\bf \tilde L}` ``WaveFunctionComponent::L``. The sum over components are stored in ``TrialWaveFunction::G`` and ``TrialWaveFunction::L``.
The :math:`\frac{{\nabla ^2} \psi}{\psi}` needed by kinetic part of the local energy can be computed as

.. math::
     \frac{\nabla^2_i \psi}{\psi} = {\bf L_i} + {\bf G}_i \cdot {\bf G}_i

see ``QMCHamiltonians/BareKineticEnergy.h``.

Wavefunction evaluation
~~~~~~~~~~~~~~~~~~~~~~~

The process for creating a new wavefunction component class is to derive
from WaveFunctionComponent and implement a number pure virtual
functions. To start most of them can be empty.

The following four functions evaluate the wavefunction values and
spatial derivatives:

``evaluateLog`` Computes the log of the wavefunction and the gradient
and Laplacian (of the log of the wavefunction) for all particles. The
input is the\ ``ParticleSet``\ (``P``) (of the electrons). The log of
the wavefunction should be stored in the ``LogValue`` member variable,
and used as the return value from the function.  The gradient is stored
in ``G`` and the Laplacian in ``L``.

``ratio`` Computes the wavefunction ratio (not the log) for a single
particle move (:math:`\psi_{new}/\psi_{old}`). The inputs are the
``ParticleSet``\ (``P``) and the particle index (``iat``).

``evalGrad`` Computes the gradient for a given particle. The inputs are
the ``ParticleSet``\ (``P``) and the particle index (``iat``).

``ratioGrad`` Computes the wavefunction ratio and the gradient at the
new position for a single particle move. The inputs are the
``ParticleSet``\ (``P``) and the particle index (``iat``). The output
gradient is in ``grad_iat``;

The ``updateBuffer`` function needs to be implemented, but to start it
can simply call ``evaluateLog``.

The ``put`` function should be implemented to read parameter specifics
from the input XML file.

Function use
~~~~~~~~~~~~

For debugging it can be helpful to know the under what conditions the
various routines are called.

The VMC and DMC loops initialize the walkers by calling ``evaluateLog``.
For all-electron moves, each timestep advance calls ``evaluateLog``. If
the ``use_drift`` parameter is no, then only the wavefunction value is
used for sampling. The gradient and Laplacian are used for computing the
local energy.

For particle-by-particle moves, each timestep advance

#. calls ``evalGrad``

#. computes a trial move

#. calls ``ratioGrad`` for the wavefunction ratio and the gradient at
   the trial position. (If the ``use_drift`` parameter is no, the
   ``ratio`` function is called instead.)

The following example shows part of an input block for VMC with
all-electron moves and drift.

::

   <qmc method="vmc" target="e" move="alle">
     <parameter name="use_drift">yes</parameter>
   </qmc>

Particle distances
~~~~~~~~~~~~~~~~~~

The ``ParticleSet`` parameter in these functions refers to the
electrons. The distance tables that store the inter-particle distances
are stored as an array.

To get the electron-ion distances, add the ion ``ParticleSet`` using
``addTable`` and save the returned index. Use that index to get the
ion-electron distance table.

::

   const int ei_id = elecs.addTable(ions); // in the constructor only
   const auto& ei_table = elecs.getDistTable(ei_id); // when consuming a distance table

Getting the electron-electron distances is very similar, just add the
electron ``ParticleSet`` using ``addTable``.

Only the lower triangle for the electron-electron table should be used.
It is the only part of the distance table valid throughout the run.
During particle-by-particle move, there are extra restrictions. When a
move of electron iel is proposed, only the lower triangle parts
[0,iel)[0,iel) [iel, Nelec)[iel, Nelec) and the row [iel][0:Nelec) are
valid. In fact, the current implementation of distance based two and
three body Jastrow factors in QMCPACK only needs the row [iel][0:Nelec).

In ``ratioGrad``, the new distances are stored in the ``Temp_r`` and
``Temp_dr`` members of the distance tables.

Setup
~~~~~

A builder processes XML input, creates the wavefunction, and adds it to
``targetPsi``. Builders derive from ``WaveFunctionComponentBuilder``.

The new builder hooks into the XML processing in
``WaveFunctionFactory.cpp`` in the ``build`` function.

Caching values
~~~~~~~~~~~~~~

The ``acceptMove`` and ``restore`` methods are called on accepted and
rejected moves for the component to update cached values.

Threading
~~~~~~~~~

The ``makeClone`` function needs to be implemented to work correctly
with OpenMP threading. There will be one copy of the component created
for each thread. If there is no extra storage, calling the copy
constructor will be sufficient. If there are cached values, the clone
call may need to create space.

Parameter optimization
~~~~~~~~~~~~~~~~~~~~~~

The ``checkInVariables``, ``checkOutVariables``, and ``resetParameters``
functions manage the variational parameters. Optimizable variables also
need to be registered when the XML is processed.


Variational parameter derivatives are computed in the
``evaluateDerivatives`` function. It computes the derivatives of both the log
of the wavefunction and kinetic energy with respect to optimizable parameters
and adds the results to the corresponding output arrays.

The kinetic energy derivatives are computed as

.. math::
  \sum_i -\frac{1}{2 m_i}({\partial}_\alpha {\bf L}_i + 2 {\bf G}_i \cdot {\partial}_\alpha {\bf G}_i)

with each ``WaveFunctionComponent`` contributing

.. math::
  -\frac{1}{2}{\partial}_\alpha \tilde L - G \cdot {\partial}_\alpha \tilde G

Right now :math:`1/m` factor is applied in ``TrialWaveFunction``.
This is a bug when the particle set doesn't hold equal mass particles.

.. _helium-wavefunction-example:

Helium Wavefunction Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The code contains an example of a wavefunction component for a Helium atom using STO orbitals and a Pade Jastrow.

to
The wavefunction is

.. math::
  :label: eq265

  \psi = \frac{1}{\sqrt{\pi}} \exp(-Z r_1) \exp(-Z r_2) \exp(A / (1 + B r_{12}))

where :math:`Z = 2` is the nuclear charge, :math:`A=1/2` is the
electron-electron cusp, and :math:`B` is a variational parameter. The
electron-ion distances are :math:`r_1` and :math:`r_2`, and
:math:`r_{12}` is the electron-electron distance. The wavefunction is
the same as the one expressed with built-in components in
``examples/molecules/He/he_simple_opt.xml``.

The code is in ``src/QMCWaveFunctions/ExampleHeComponent.cpp``. The
builder is in ``src/QMCWaveFunctions/ExampleHeBuilder.cpp``. The input
file is in ``examples/molecules/He/he_example_wf.xml``. A unit test
compares results from the wavefunction evaluation functions for
consistency in ``src/QMCWaveFunctions/tests/test_example_he.cpp``.

The recommended approach for creating a new wavefunction component is to
copy the example and the unit test. Implement the evaluation functions
and ensure the unit test passes.

Linear Algebra
--------------

Like in many methods which solve the Schrödinger equation, linear
algebra plays a critical role in QMC algorithms and thus is crucial to
the performance of QMCPACK. There are a few components in QMCPACK use
BLAS/LAPACK with their own characteristics.

Real space QMC
~~~~~~~~~~~~~~

Single particle orbitals
^^^^^^^^^^^^^^^^^^^^^^^^

Spline evaluation as commonly used in solid-state simulations does not use any dense linear algebra library calls.
LCAO evaluation as commonly used in molecular calculations relies on BLAS2 GEMV to compute SPOs from a basis set.

Slater determinants
^^^^^^^^^^^^^^^^^^^

Slater determinants are calculated on :math:`N \times N` Slater
matrices. :math:`N` is the number of electrons for a given spin. In the
actually implementation, operations on the inverse matrix of Slater
matrix for each walker dominate the computation. To initialize it,
DGETRF and DGETRI from LAPACK are called. The inverse matrix can be
stored out of place. During random walking, inverse matrices are updated
by either Sherman-Morrison rank-1 update or delayed update. Update
algorithms heavily relies on BLAS. All the BLAS operations require
S,C,D,Z cases.

Sherman-Morrison rank-1 update uses BLAS2 GEMV and GER on
:math:`N \times N` matrices.

Delayed rank-K update uses

-  BLAS1 SCOPY on :math:`N` array.

-  BLAS2 GEMV, GER on :math:`k \times N` and :math:`k \times k`
   matrices. :math:`k` ranges from 1 to :math:`K` when updates are
   delayed and accumulated.

-  BLAS3 GEMM at the final update.

   -  ’T’, ’N’, K, N, N

   -  ’N’, ’N’, N, K, K

   -  ’N’, ’N’, N, N, K

The optimal K depends on the hardware but it usually ranges from 32 to
256.

QMCPACK solves systems with a few to thousands of electrons. To make all
the BLAS/LAPACK operation efficient on accelerators. Batching is needed
and optimized for :math:`N < 2000`. Non-batched functions needs to be
optimized for :math:`N > 500`. Note: 2000 and 500 are only rough
estimates.

Wavefunction optimizer
^^^^^^^^^^^^^^^^^^^^^^

to be added.

Auxiliary field QMC
~~~~~~~~~~~~~~~~~~~

The AFQMC implementation in QMCPACK relies heavily on linear algebra operations from BLAS/LAPACK. The performance of the code is netirely dependent on the performance of these libraries. See below for a detailed list of the main routines used from BLAS/LAPACK. Since the AFQMC code can work with both single and double precision builds, all 4 versions of these routines (S,C,D,Z) are generally needed, for this reason we omit the data type label.

-  BLAS1: SCAL, COPY, DOT, AXPY

-  BLAS2: GEMV, GER

-  BLAS3: GEMM

-  LAPACK: GETRF, GETRI, GELQF, UNGLQ, ORGLQ, GESVD, HEEVR, HEGVX

While the dimensions of the matrix operations will depend entirely on
the details of the calculation, typical matrix dimensions range from the
100s, for small system sizes, to over 20000 for the largest calculations
attempted so far. For builds with GPU accelerators, we make use of
batched and strided implementations of these routines. Batched
implementations of GEMM, GETRF, GETRI, GELQF and UNGLQ are particularly
important for the performance of the GPU build on small to medium size
problems. Batched implementations of DOT, AXPY and GEMV would also be
quite useful, but they are not yet generally available. On GPU builds,
the code uses batched implementations of these routines when available
by default.

Slater-backflow wavefunction implementation details
---------------------------------------------------

For simplicity, consider :math:`N` identical fermions of the same spin
(e.g., up electrons) at spatial locations
:math:`\{\mathbf{r}_1,\mathbf{r}_2,\dots,\mathbf{r}_{N}\}`. Then the
Slater determinant can be written as

.. math::
  :label: eq245


   S=\det M\:,

where each entry in the determinant is an SPO evaluated at a particle
position

.. math::
  :label: eq246

   \begin{aligned}
   M_{ij} = \phi_i(\mathbf{r}_j)\:.\end{aligned}

When backflow transformation is applied to the determinant, the particle
coordinates :math:`\mathbf{r}_i` that go into the SPOs are replaced by
quasi-particle coordinates :math:`\mathbf{x}_i`:

.. math::
  :label: eq247

   \begin{aligned}
   M_{ij} = \phi_i(\mathbf{x}_j)\:, \end{aligned}

where

.. math::
  :label: eq248

   \begin{aligned}
   \mathbf{x}_i=\mathbf{r}_i+\sum\limits_{j=1,j\neq i}^N\eta(r_{ij})(\mathbf{r}_i-\mathbf{r}_j)\:. \end{aligned}

:math:`r_{ij}=\vert\mathbf{r}_i-\mathbf{r}_j\vert`. The integers i,j
label the particle/quasi-particle. There is a one-to-one correspondence
between the particles and the quasi-particles, which is simplest when
:math:`\eta=0`.

Value
~~~~~

The evaluation of the Slater-backflow wavefunction is almost identical
to that of a Slater wavefunction. The only difference is that the
quasi-particle coordinates are used to evaluate the SPOs. The actual
value of the determinant is stored during the inversion of the matrix
:math:`M` (``cgetrf``\ :math:`\rightarrow`\ ``cgetri``). Suppose
:math:`M=LU`, then :math:`S=\prod\limits_{i=1}^N L_{ii} U_{ii}`.

::

  // In DiracDeterminantWithBackflow::evaluateLog(P,G,L)
  Phi->evaluate(BFTrans->QP, FirstIndex, LastIndex, psiM,dpsiM,grad_grad_psiM);
  psiMinv = psiM;
  LogValue=InvertWithLog(psiMinv.data(),NumPtcls,NumOrbitals
    ,WorkSpace.data(),Pivot.data(),PhaseValue);

QMCPACK represents the complex value of the wavefunction in polar
coordinates :math:`S=e^Ue^{i\theta}`. Specifically, ``LogValue``
:math:`U` and ``PhaseValue`` :math:`\theta` are handled separately. In
the following, we will consider derivatives of the log value only.

Gradient
~~~~~~~~

To evaluate particle gradient of the log value of the Slater-backflow
wavefunction, we can use the :math:`\log\det` identity in
:eq:`eq249`. This identity maps the derivative of
:math:`\log\det M` with respect to a real variable :math:`p` to a trace
over :math:`M^{-1}dM`:

.. math::
  :label: eq249

   \begin{aligned}
   \frac{\partial}{\partial p}\log\det M = \text{tr}\left( M^{-1} \frac{\partial M}{\partial p} \right) .\end{aligned}

Following Kwon, Ceperley, and
Martin :cite:`Kwon1993backflow`, the particle gradient

.. math::
  :label: eq250

   \begin{aligned}
   G_i^\alpha \equiv \frac{\partial}{\partial r_i^\alpha} \log\det M = \sum\limits_{j=1}^N \sum\limits_{\beta=1}^3 F_{jj}^\beta A_{jj}^{\alpha\beta}\:, \end{aligned}

where the quasi-particle gradient matrix

.. math::
  :label: eq251

   \begin{aligned}
   A_{ij}^{\alpha\beta} \equiv \frac{\partial x_j^\beta}{\partial r_i^\alpha}\:,\end{aligned}

and the intermediate matrix

.. math::
  :label: eq252

   \begin{aligned}
   F_{ij}^\alpha\equiv\sum\limits_k M^{-1}_{ik} dM_{kj}^\alpha\:,\end{aligned}

with the SPO derivatives (w.r. to quasi-particle coordinates)

.. math::
  :label: eq253

   \begin{aligned}
   dM_{ij}^\alpha \equiv \frac{\partial M_{ij}}{\partial x_j^\alpha}\:.\end{aligned}

Notice that we have made the name change of :math:`\phi\rightarrow M`
from the notations of ref. :cite:`Kwon1993backflow`. This
name change is intended to help the reader associate M with the QMCPACK
variable ``psiM``.

::

  // In DiracDeterminantWithBackflow::evaluateLog(P,G,L)
  for(int i=0; i<num; i++) // k in above formula
  {
    for(int j=0; j<NumPtcls; j++)
    {
      for(int k=0; k<OHMMS_DIM; k++) // alpha in above formula
      {
        myG(i) += dot(BFTrans->Amat(i,FirstIndex+j),Fmat(j,j));
      }
    }
  }

:eq:`eq250` is still relatively simple to understand. The
:math:`A` matrix maps changes in particle coordinates
:math:`d\mathbf{r}` to changes in quasi-particle coordinates
:math:`d\mathbf{x}`. Dotting A into F propagates :math:`d\mathbf{x}` to
:math:`dM`. Thus :math:`F\cdot A` is the term inside the trace operator
of :eq:`eq249`. Finally, performing the trace completes the
evaluation of the derivative.

Laplacian
~~~~~~~~~

The particle Laplacian is given in
:cite:`Kwon1993backflow` as

.. math::
  :label: eq254

   \begin{aligned}
   L_i \equiv \sum\limits_{\beta} \frac{\partial^2}{\partial (r_i^\beta)^2} \log\det M = \sum\limits_{j\alpha} B_{ij}^\alpha F_{jj}^\alpha - \sum\limits_{jk}\sum\limits_{\alpha\beta\gamma} A_{ij}^{\alpha\beta}A_{ik}^{\alpha\gamma}\times\left(F_{kj}^\alpha F_{jk}^\gamma -\delta_{jk}\sum\limits_m M^{-1}_{jm} d2M_{mj}^{\beta\gamma}\right), \end{aligned}

where the quasi-particle Laplacian matrix

.. math::
  :label: eq255

   \begin{aligned}
   B_{ij}^{\alpha} \equiv \sum\limits_\beta \frac{\partial^2 x_j^\alpha}{\partial (r_i^\beta)^2}\:,\end{aligned}

with the second derivatives of the single-particles orbitals being

.. math::
  :label: eq256

   \begin{aligned}
   d2M_{ij}^{\alpha\beta} \equiv \frac{\partial^2 M_{ij}}{\partial x_j^\alpha\partial x_j^\beta}\:.\end{aligned}

Schematically, :math:`L_i` has contributions from three terms of the
form :math:`BF, AAFF, and tr(AA,Md2M)`, respectively.
:math:`A, B, M ,d2M,` and :math:`F` can be calculated and stored before
the calculations of :math:`L_i`. The first :math:`BF` term can be
directly calculated in a loop over quasi-particle coordinates
:math:`j\alpha`.

::

  // In DiracDeterminantWithBackflow::evaluateLog(P,G,L)
  for(int j=0; j<NumPtcls; j++)
    for(int a=0; a<OHMMS_DIM; k++)
      myL(i) += BFTrans->Bmat_full(i,FirstIndex+j)[a]*Fmat(j,j)[a];

Notice that :math:`B_{ij}^\alpha` is stored in ``Bmat_full``, NOT
``Bmat``.

The remaining two terms both involve :math:`AA`. Thus, it is best to
define a temporary tensor :math:`AA`:

.. math::
  :label: eq257

   \begin{aligned}
   {}_iAA_{jk}^{\beta\gamma} \equiv \sum\limits_\alpha A_{ij}^{\alpha\beta} A_{ij}^{\alpha\gamma}\:,\end{aligned}

which we will overwrite for each particle :math:`i`. Similarly, define
:math:`FF`:

.. math::
  :label: eq258

   \begin{aligned}
   FF_{jk}^{\alpha\gamma} \equiv F_{kj}^\alpha F_{jk}^\gamma\:,\end{aligned}

which is simply the outer product of :math:`F\otimes F`. Then the
:math:`AAFF` term can be calculated by fully contracting :math:`AA` with
:math:`FF`.

::

  // In DiracDeterminantWithBackflow::evaluateLog(P,G,L)
  for(int j=0; j<NumPtcls; j++)
    for(int k=0; k<NumPtcls; k++)
      for(int i=0; i<num; i++)
      {
        Tensor<RealType,OHMMS_DIM> AA = dot(transpose(BFTrans->Amat(i,FirstIndex+j)),BFTrans->Amat(i,FirstIndex+k));
        HessType FF = outerProduct(Fmat(k,j),Fmat(j,k));
        myL(i) -= traceAtB(AA,FF);
      }

Finally, define the SPO derivative term:

.. math::
  :label: eq259

   \begin{aligned}
   Md2M_j^{\beta\gamma} \equiv \sum\limits_m M^{-1}_{jm} d2M_{mj}^\beta\:,\end{aligned}

then the last term is given by the contraction of :math:`Md2M` (``q_j``)
with the diagonal of :math:`AA`.

::

  for(int j=0; j<NumPtcls; j++)
  {
    HessType q_j;
    q_j=0.0;
    for(int k=0; k<NumPtcls; k++)
      q_j += psiMinv(j,k)*grad_grad_psiM(j,k);
    for(int i=0; i<num; i++)
    {
      Tensor<RealType,OHMMS_DIM> AA = dot(
        transpose(BFTrans->Amat(i,FirstIndex+j)),
        BFTrans->Amat(i,FirstIndex+j)
      );
      myL(i) += traceAtB(AA,q_j);
    }
  }

Wavefunction parameter derivative
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To use the robust linear optimization method of
:cite:`Toulouse2007linear`, the trial wavefunction
needs to know its contributions to the overlap and hamiltonian matrices.
In particular, we need derivatives of these matrices with respect to
wavefunction parameters. As a consequence, the wavefunction :math:`\psi`
needs to be able to evaluate
:math:`\frac{\partial}{\partial p} \ln \psi` and
:math:`\frac{\partial}{\partial p} \frac{\mathcal{H}\psi}{\psi}`, where
:math:`p` is a parameter.

When 2-body backflow is considered, a wavefunction parameter :math:`p`
enters the :math:`\eta` function only (equation :eq:`eq248`).
:math:`\mathbf{r}`, :math:`\phi`, and :math:`M` do not explicitly
dependent on :math:`p`. Derivative of the log value is almost identical
to particle gradient. Namely, :eq:`eq250` applies upon the
substitution :math:`r_i^\alpha\rightarrow p`.

.. math::
  :label: eq260

   \begin{aligned}
   \frac{\partial}{\partial p} \ln\det M = \sum\limits_{j=1}^N \sum\limits_{\beta=1}^3 F_{jj}^\beta \left({}_pC_{j}^{\beta}\right)\:,\end{aligned}

where the quasi-particle derivatives are stored in ``Cmat``

.. math::
  :label: eq261

   \begin{aligned}
   {}_pC_{i}^{\alpha} \equiv \frac{\partial}{\partial p} x_{i}^{\alpha}\:.\end{aligned}

The change in local kinetic energy is a lot more difficult to calculate

.. math::
  :label: eq262

   \begin{aligned}
   \frac{\partial T_{\text{local}}}{\partial p} = \frac{\partial}{\partial p} \left\{ \left( \sum\limits_{i=1}^N \frac{1}{2m_i} \nabla^2_i \right) \ln \det M \right\} = \sum\limits_{i=1}^N \frac{1}{2m_i} \frac{\partial}{\partial p} L_i\:, \end{aligned}

where :math:`L_i` is the particle Laplacian defined in
:eq:`eq254` To evaluate :eq:`eq262`, we need to
calculate parameter derivatives of all three terms defined in the
Laplacian evaluation. Namely :math:`(B)(F)`, :math:`(AA)(FF)`, and
:math:`\text{tr}(AA,Md2M)`, where we have put parentheses around previously
identified data structures. After :math:`\frac{\partial}{\partial p}`
hits, each of the three terms will split into two terms by the product
rule. Each smaller term will contain a contraction of two data
structures. Therefore, we will need to calculate the parameter
derivatives of each data structure defined in the Laplacian evaluation:

.. math::
  :label: eq263

   \begin{aligned}
   {}_pX_{ij}^{\alpha\beta} \equiv \frac{\partial}{\partial p} A_{ij}^{\alpha\beta}\:, \\
   {}_pY_{ij}^{\alpha} \equiv \frac{\partial}{\partial p} B_{ij}^{\alpha}\:, \\
   {}_pdF_{ij}^{\alpha} \equiv \frac{\partial}{\partial p} F_{ij}^{\alpha}\:, \\
   {}_{pi}{AA'}_{jk}^{\beta\gamma} \equiv \frac{\partial}{\partial p}  {}_iAA_{jk}^{\beta\gamma}\:, \\
   {}_p {FF'}_{jk}^{\alpha\gamma} \equiv \frac{\partial}{\partial p} FF_{jk}^{\alpha\gamma}\:, \\
   {}_p {Md2M'}_{j}^{\beta\gamma} \equiv \frac{\partial}{\partial p} Md2M_j^{\beta\gamma}\:.\end{aligned}

X and Y are stored as ``Xmat`` and ``Ymat_full`` (NOT ``Ymat``) in the
code. dF is ``dFa``. :math:`AA'` is not fully stored; intermediate
values are stored in ``Aij_sum`` and ``a_j_sum``. :math:`FF'` is
calculated on the fly as :math:`dF\otimes F+F\otimes dF`. :math:`Md2M'`
is not stored; intermediate values are stored in ``q_j_prime``.

Scalar estimator implementation
-------------------------------

Introduction: Life of a specialized OperatorBase
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Almost all observables in QMCPACK are implemented as specialized derived classes of the OperatorBase base class. Each observable
is instantiated in HamiltonianFactory and added to QMCHamiltonian for tracking. QMCHamiltonian tracks two types of observables:
main and auxiliary. Main observables contribute to the local energy. These observables are elements of the simulated Hamiltonian
such as kinetic or potential energy. Auxiliary observables are expectation values of matrix elements that do not contribute to the
local energy. These Hamiltonians do not affect the dynamics of the simulation. In the code, the main observables are labeled by
“physical” flag; the auxiliary observables have“physical” set to false.

Initialization
^^^^^^^^^^^^^^

When an ``<estimator type="est_type" name="est_name" other_stuff="value"/>`` tag is present in the ``<hamiltonian/>`` section, it is first read by HamiltonianFactory. In general, the ``type`` of the estimator will determine which specialization of OperatorBase should be instantiated, and a derived class with ``myName="est_name"`` will be constructed. Then, the put() method of this specific class will be called to read any other parameters in the ``<estimator/>`` XML node. Sometimes these parameters will instead be read by HamiltonianFactory because it can access more objects than OperatorBase.

Cloning
^^^^^^^

When ``OpenMP`` threads are spawned, the estimator will be cloned by the ``CloneManager``, which is a parent class of many QMC drivers.

::

  // In CloneManager.cpp
  #pragma omp parallel for shared(w,psi,ham)
  for(int ip=1; ip<NumThreads; ++ip)
  {
    wClones[ip]=new MCWalkerConfiguration(w);
    psiClones[ip]=psi.makeClone(*wClones[ip]);
    hClones[ip]=ham.makeClone(*wClones[ip],*psiClones[ip]);
  }

In the preceding snippet, ``ham`` is the reference to the estimator on the master thread. If the implemented estimator does not allocate memory for any array, then the default constructor should suffice for the ``makeClone`` method.

::

  // In SpeciesKineticEnergy.cpp
  OperatorBase* SpeciesKineticEnergy::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return new SpeciesKineticEnergy(*this);
  }

If memory is allocated during estimator construction (usually when parsing the XML node in the ``put`` method), then the ``makeClone`` method should perform the same initialization on each thread.

::

  OperatorBase* LatticeDeviationEstimator::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    LatticeDeviationEstimator* myclone = new LatticeDeviationEstimator(qp,spset,tgroup,sgroup);
    myclone->put(input_xml);
    return myclone;
  }

Evaluate
^^^^^^^^

After the observable class (derived class of OperatorBase) is constructed and prepared (by the put() method), it is ready to be used in a QMCDriver. A QMCDriver will call ``H.auxHevaluate(W,thisWalker)`` after every accepted move, where H is the QMCHamiltonian that holds all main and auxiliary Hamiltonian elements, W is a MCWalkerConfiguration, and thisWalker is a pointer to the current walker being worked on. As shown in the following, this function goes through each auxiliary Hamiltonian element and evaluates it using the current walker configuration. Under the hood, observables are calculated and dumped to the main particle set's property list for later collection.

::

  // In QMCHamiltonian.cpp
  // This is more efficient.
  // Only calculate auxH elements if moves are accepted.
  void QMCHamiltonian::auxHevaluate(ParticleSet& P, Walker_t& ThisWalker)
  {
  #if !defined(REMOVE_TRACEMANAGER)
    collect_walker_traces(ThisWalker,P.current_step);
  #endif
    for(int i=0; i<auxH.size(); ++i)
    {
      auxH[i]->setHistories(ThisWalker);
      RealType sink = auxH[i]->evaluate(P);
      auxH[i]->setObservables(Observables);
  #if !defined(REMOVE_TRACEMANAGER)
      auxH[i]->collect_scalar_traces();
  #endif
      auxH[i]->setParticlePropertyList(P.PropertyList,myIndex);
    }
  }

For estimators that contribute to the local energy (main observables), the return value of evaluate() is used in accumulating the local energy. For auxiliary estimators, the return value is not used (``sink`` local variable above); only the value of Value is recorded property lists by the setObservables() method as shown in the preceding code snippet. By default, the setObservables() method will transfer ``auxH[i]->Value`` to ``P.PropertyList[auxH[i]->myIndex]``. The same property list is also kept by the particle set being moved by QMCDriver. This list is updated by ``auxH[i]->setParticlePropertyList(P.PropertyList,myIndex)``, where myIndex is the starting index of space allocated to this specific auxiliary Hamiltonian in the property list kept by the target particle set P.

Collection
^^^^^^^^^^

The actual statistics are collected within the QMCDriver, which owns
an EstimatorManager object. This object (or a clone in the case of
multithreading) will be registered with each mover it owns. For each mover
(such as VMCUpdatePbyP derived from QMCUpdateBase), an accumulate() call
is made, which by default, makes an accumulate(walkerset) call to the
EstimatorManager it owns. Since each walker has a property set, EstimatorManager uses that local copy to calculate statistics. The EstimatorManager performs block averaging and file I/O.

Single scalar estimator implementation guide
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Almost all of the defaults can be used for a single scalar observable. With any luck, only the put() and evaluate() methods need to be implemented. As an example, this section presents a step-by-step guide for implementing a \verb|SpeciesKineticEnergy| estimator that calculates the kinetic energy of a specific species instead of the entire particle set. For example, a possible input to this estimator can be:

``<estimator type="specieskinetic" name="ukinetic" group="u"/>``

``<estimator type="specieskinetic" name="dkinetic" group="d"/>``.

This should create two extra columns in the ``scalar.dat`` file that contains the kinetic energy of the up and down electrons in two separate columns. If the estimator is properly implemented, then the sum of these two columns should be equal to the default ``Kinetic`` column.

Barebone
^^^^^^^^

The first step is to create a barebone class structure for this simple scalar estimator. The goal is to be able to instantiate this scalar estimator with an XML node and have it print out a column of zeros in the ``scalar.dat`` file.

To achieve this, first create a header file "SpeciesKineticEnergy.h" in the QMCHamiltonians folder, with only the required functions declared as follows:

::

  // In SpeciesKineticEnergy.h
  #ifndef QMCPLUSPLUS_SPECIESKINETICENERGY_H
  #define QMCPLUSPLUS_SPECIESKINETICENERGY_H

  #include <Particle/WalkerSetRef.h>
  #include <QMCHamiltonians/OperatorBase.h>

  namespace qmcplusplus
  {

  class SpeciesKineticEnergy: public OperatorBase
  {
  public:

    SpeciesKineticEnergy(ParticleSet& P):tpset(P){ };

    bool put(xmlNodePtr cur);         // read input xml node, required
    bool get(std::ostream& os) const; // class description, required

    Return_t evaluate(ParticleSet& P);
    inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
    { // delegate responsity inline for speed
      return evaluate(P);
    }

    // pure virtual functions require overrider
    void resetTargetParticleSet(ParticleSet& P) { }                         // required
    OperatorBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi); // required

  private:
    ParticleSet& tpset;

  }; // SpeciesKineticEnergy

  } // namespace qmcplusplus
  #endif

Notice that a local reference ``tpset`` to the target particle set ``P`` is saved in the constructor. The target particle set carries much information useful for calculating observables. Next, make "SpeciesKineticEnergy.cpp," and make vacuous definitions.

::

  // In SpeciesKineticEnergy.cpp
  #include <QMCHamiltonians/SpeciesKineticEnergy.h>
  namespace qmcplusplus
  {

  bool SpeciesKineticEnergy::put(xmlNodePtr cur)
  {
    return true;
  }

  bool SpeciesKineticEnergy::get(std::ostream& os) const
  {
    return true;
  }

  SpeciesKineticEnergy::Return_t SpeciesKineticEnergy::evaluate(ParticleSet& P)
  {
    Value = 0.0;
    return Value;
  }

  OperatorBase* SpeciesKineticEnergy::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    // no local array allocated, default constructor should be enough
    return new SpeciesKineticEnergy(*this);
  }

  } // namespace qmcplusplus

Now, head over to HamiltonianFactory and instantiate this observable if an XML node is found requesting it. Look for "gofr" in HamiltonianFactory.cpp, for example, and follow the if block.

::

  // In HamiltonianFactory.cpp
  #include <QMCHamiltonians/SpeciesKineticEnergy.h>
  else if(potType =="specieskinetic")
  {
    SpeciesKineticEnergy* apot = new SpeciesKineticEnergy(*target_particle_set);
    apot->put(cur);
    targetH->addOperator(apot,potName,false);
  }

The last argument of addOperator() (i.e., the ``false`` flag) is **crucial**. This tells QMCPACK that the observable we implemented is not a physical Hamiltonian; thus, it will not contribute to the local energy. Changes to the local energy will alter the dynamics of the simulation. Finally, add "SpeciesKineticEnergy.cpp" to HAMSRCS in "CMakeLists.txt" located in the QMCHamiltonians folder. Now, recompile QMCPACK and run it on an input that requests ``<estimator type="specieskinetic" name="ukinetic"/>`` in the ``hamiltonian`` block. A column of zeros should appear in the ``scalar.dat`` file under the name "ukinetic."

Evaluate
^^^^^^^^

The evaluate() method is where we perform the calculation of the desired observable. In a first iteration, we will simply hard-code the name and mass of the particles.

::

  // In SpeciesKineticEnergy.cpp
  #include <QMCHamiltonians/BareKineticEnergy.h> // laplaician() defined here
  SpeciesKineticEnergy::Return_t SpeciesKineticEnergy::evaluate(ParticleSet& P)
  {
    std::string group="u";
    RealType minus_over_2m = -0.5;

    SpeciesSet& tspecies(P.getSpeciesSet());

    Value = 0.0;
    for (int iat=0; iat<P.getTotalNum(); iat++)
    {
      if (tspecies.speciesName[ P.GroupID(iat) ] == group)
      {
        Value += minus_over_2m*laplacian(P.G[iat],P.L[iat]);
      }
    }
    return Value;

    // Kinetic column has:
    // Value = -0.5*( Dot(P.G,P.G) + Sum(P.L) );
  }

*Voila*---you should now be able to compile QMCPACK, rerun, and see that the values in the "ukinetic" column are no longer zero. Now, the only task left to make this basic observable complete is to read in the extra parameters instead of hard-coding them.

Parse extra input
^^^^^^^^^^^^^^^^^

The preferred method to parse extra input parameters in the XML node is
to implement the put() function of our specific observable. Suppose we
wish to read in a single string that tells us whether to record the
kinetic energy of the up electron (group=“u") or the down electron
(group=“d"). This is easily achievable using the OhmmsAttributeSet
class,

::

  // In SpeciesKineticEnergy.cpp
  #include <OhmmsData/AttributeSet.h>
  bool SpeciesKineticEnergy::put(xmlNodePtr cur)
  {
    // read in extra parameter "group"
    OhmmsAttributeSet attrib;
    attrib.add(group,"group");
    attrib.put(cur);

    // save mass of specified group of particles
    SpeciesSet& tspecies(tpset.getSpeciesSet());
    int group_id  = tspecies.findSpecies(group);
    int massind   = tspecies.getAttribute("mass");
    minus_over_2m = -1./(2.*tspecies(massind,group_id));

    return true;
  }

where we may keep "group" and "minus\_over\_2m" as local variables to our specific class.

::

  // In SpeciesKineticEnergy.h
  private:
    ParticleSet& tpset;
    std::string  group;
    RealType minus_over_2m;

Notice that the previous operations are made possible by the saved reference to the target particle set. Last but not least, compile and run a full example (i.e., a short DMC calculation) with the following XML nodes in your input:

``<estimator type="specieskinetic" name="ukinetic" group="u"/>``

``<estimator type="specieskinetic" name="dkinetic" group="d"/>``

Make sure the sum of the “ukinetic" and “dkinetic" columns is
**exactly** the same as the Kinetic columns at **every block**.

For easy reference, a summary of the complete list of changes follows:

::

  // In HamiltonianFactory.cpp
  #include "QMCHamiltonians/SpeciesKineticEnergy.h"
  else if(potType =="specieskinetic")
  {
  	SpeciesKineticEnergy* apot = new SpeciesKineticEnergy(*targetPtcl);
  	apot->put(cur);
  	targetH->addOperator(apot,potName,false);
  }

::

  // In SpeciesKineticEnergy.h
  #include <Particle/WalkerSetRef.h>
  #include <QMCHamiltonians/OperatorBase.h>

  namespace qmcplusplus
  {

  class SpeciesKineticEnergy: public OperatorBase
  {
  public:

    SpeciesKineticEnergy(ParticleSet& P):tpset(P){ };

    // xml node is read by HamiltonianFactory, eg. the sum of following should be equivalent to Kinetic
    // <estimator name="ukinetic" type="specieskinetic" target="e" group="u"/>
    // <estimator name="dkinetic" type="specieskinetic" target="e" group="d"/>
    bool put(xmlNodePtr cur);         // read input xml node, required
    bool get(std::ostream& os) const; // class description, required

    Return_t evaluate(ParticleSet& P);
    inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
    { // delegate responsity inline for speed
      return evaluate(P);
    }

    // pure virtual functions require overrider
    void resetTargetParticleSet(ParticleSet& P) { }                         // required
    OperatorBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi); // required

  private:
    ParticleSet&       tpset; // reference to target particle set
    std::string        group; // name of species to track
    RealType   minus_over_2m; // mass of the species !! assume same mass
    // for multiple species, simply initialize multiple estimators

  }; // SpeciesKineticEnergy

  } // namespace qmcplusplus
  #endif

::

  // In SpeciesKineticEnergy.cpp
  #include <QMCHamiltonians/SpeciesKineticEnergy.h>
  #include <QMCHamiltonians/BareKineticEnergy.h> // laplaician() defined here
  #include <OhmmsData/AttributeSet.h>

  namespace qmcplusplus
  {

  bool SpeciesKineticEnergy::put(xmlNodePtr cur)
  {
    // read in extra parameter "group"
    OhmmsAttributeSet attrib;
    attrib.add(group,"group");
    attrib.put(cur);

    // save mass of specified group of particles
    int group_id  = tspecies.findSpecies(group);
    int massind   = tspecies.getAttribute("mass");
    minus_over_2m = -1./(2.*tspecies(massind,group_id));

    return true;
  }

  bool SpeciesKineticEnergy::get(std::ostream& os) const
  { // class description
    os << "SpeciesKineticEnergy: " << myName << " for species " << group;
    return true;
  }

  SpeciesKineticEnergy::Return_t SpeciesKineticEnergy::evaluate(ParticleSet& P)
  {
    Value = 0.0;
    for (int iat=0; iat<P.getTotalNum(); iat++)
    {
      if (tspecies.speciesName[ P.GroupID(iat) ] == group)
      {
        Value += minus_over_2m*laplacian(P.G[iat],P.L[iat]);
      }
    }
    return Value;
  }

  OperatorBase* SpeciesKineticEnergy::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  { //default constructor
    return new SpeciesKineticEnergy(*this);
  }

  } // namespace qmcplusplus

Multiple scalars
~~~~~~~~~~~~~~~~

It is fairly straightforward to create more than one column in the ``scalar.dat`` file with a single observable class. For example, if we want a single SpeciesKineticEnergy estimator to simultaneously record the kinetic energies of all species in the target particle set, we only have to write two new methods: addObservables() and setObservables(), then tweak the behavior of evaluate(). First, we will have to override the default behavior of addObservables() to make room for more than one column in the ``scalar.dat`` file as follows,

::

  // In SpeciesKineticEnergy.cpp
  void SpeciesKineticEnergy::addObservables(PropertySetType& plist, BufferType& collectables)
  {
    myIndex = plist.size();
    for (int ispec=0; ispec<num_species; ispec++)
    { // make columns named "$myName_u", "$myName_d" etc.
      plist.add(myName + "_" + species_names[ispec]);
    }
  }

where “num_species” and “species_name” can be local variables
initialized in the constructor. We should also initialize some local
arrays to hold temporary data.

::

  // In SpeciesKineticEnergy.h
  private:
    int num_species;
    std::vector<std::string> species_names;
    std::vector<RealType> species_kinetic,vec_minus_over_2m;

  // In SpeciesKineticEnergy.cpp
  SpeciesKineticEnergy::SpeciesKineticEnergy(ParticleSet& P):tpset(P)
  {
    SpeciesSet& tspecies(P.getSpeciesSet());
    int massind = tspecies.getAttribute("mass");

    num_species = tspecies.size();
    species_kinetic.resize(num_species);
    vec_minus_over_2m.resize(num_species);
    species_names.resize(num_species);
    for (int ispec=0; ispec<num_species; ispec++)
    {
      species_names[ispec] = tspecies.speciesName[ispec];
      vec_minus_over_2m[ispec] = -1./(2.*tspecies(massind,ispec));
    }
  }

Next, we need to override the default behavior of ``setObservables()`` to transfer multiple values to the property list kept by the main particle set, which eventually goes into the ``scalar.dat`` file.

::

  // In SpeciesKineticEnergy.cpp
  void SpeciesKineticEnergy::setObservables(PropertySetType& plist)
  { // slots in plist must be allocated by addObservables() first
    copy(species_kinetic.begin(),species_kinetic.end(),plist.begin()+myIndex);
  }

Finally, we need to change the behavior of evaluate() to fill the local
vector “species_kinetic” with appropriate observable values.

::

  SpeciesKineticEnergy::Return_t SpeciesKineticEnergy::evaluate(ParticleSet& P)
  {
    std::fill(species_kinetic.begin(),species_kinetic.end(),0.0);

    for (int iat=0; iat<P.getTotalNum(); iat++)
    {
      int ispec = P.GroupID(iat);
      species_kinetic[ispec] += vec_minus_over_2m[ispec]*laplacian(P.G[iat],P.L[iat]);
    }

    Value = 0.0; // Value is no longer used
    return Value;
  }

That's it! The SpeciesKineticEnergy estimator no longer needs the "group" input and will automatically output the kinetic energy of every species in the target particle set in multiple columns. You should now be able to run with
``<estimator type="specieskinetic" name="skinetic"/>`` and check that the sum of all columns that start with "skinetic" is equal to the default "Kinetic" column.

HDF5 output
~~~~~~~~~~~

If we desire an observable that will output hundreds of scalars per simulation step (e.g., SkEstimator), then it is preferred to output to the ``stat.h5`` file instead of the ``scalar.dat`` file for better organization. A large chunk of data to be registered in the ``stat.h5`` file is called a "Collectable" in QMCPACK. In particular, if a OperatorBase object is initialized with ``UpdateMode.set(COLLECTABLE,1)``, then the "Collectables" object carried by the main particle set will be processed and written to the ``stat.h5`` file, where "UpdateMode" is a bit set (i.e., a collection of flags) with the following enumeration:

::

  // In OperatorBase.h
  ///enum for UpdateMode
  enum {PRIMARY=0,
    OPTIMIZABLE=1,
    RATIOUPDATE=2,
    PHYSICAL=3,
    COLLECTABLE=4,
    NONLOCAL=5,
    VIRTUALMOVES=6
  };

As a simple example, to put the two columns we produced in the previous section into the ``stat.h5`` file, we will first need to declare that our observable uses "Collectables."

::

  // In constructor add:
  hdf5_out = true;
  UpdateMode.set(COLLECTABLE,1);

Then make some room in the "Collectables" object carried by the target particle set.

::

  // In addObservables(PropertySetType& plist, BufferType& collectables) add:
  if (hdf5_out)
  {
    h5_index = collectables.size();
    std::vector<RealType> tmp(num_species);
    collectables.add(tmp.begin(),tmp.end());
  }

Next, make some room in the ``stat.h5`` file by overriding the registerCollectables() method.

::

  // In SpeciesKineticEnergy.cpp
  void SpeciesKineticEnergy::registerCollectables(std::vector<observable_helper>& h5desc, hid_t gid) const
  {
    if (hdf5_out)
    {
      std::vector<int> ndim(1,num_species);
      observable_helper h5o(myName);
      h5o.set_dimensions(ndim,h5_index);
      h5o.open(gid);
      h5desc.push_back(h5o);
    }
  }

Finally, edit evaluate() to use the space in the "Collectables" object.

::

  // In SpeciesKineticEnergy.cpp
  SpeciesKineticEnergy::Return_t SpeciesKineticEnergy::evaluate(ParticleSet& P)
  {
    RealType wgt = tWalker->Weight; // MUST explicitly use DMC weights in Collectables!
    std::fill(species_kinetic.begin(),species_kinetic.end(),0.0);

    for (int iat=0; iat<P.getTotalNum(); iat++)
    {
      int ispec = P.GroupID(iat);
      species_kinetic[ispec] += vec_minus_over_2m[ispec]*laplacian(P.G[iat],P.L[iat]);
      P.Collectables[h5_index + ispec] += vec_minus_over_2m[ispec]*laplacian(P.G[iat],P.L[iat])*wgt;
    }

    Value = 0.0; // Value is no longer used
    return Value;
  }

There should now be a new entry in the ``stat.h5`` file containing the same columns of data as the ``stat.h5`` file. After this check, we should clean up the code by

- making "hdf5_out" and input flag by editing the put() method and

- disabling output to ``scalar.dat`` when the "hdf5_out" flag is on.

Estimator output
----------------

Estimator definition
~~~~~~~~~~~~~~~~~~~~

For simplicity, consider a local property :math:`O(\bf{R})`, where
:math:`\bf{R}` is the collection of all particle coordinates. An
*estimator* for :math:`O(\bf{R})` is a weighted average over walkers:

.. math::
  :label: eq242

   \begin{aligned}
   E[O] = \left(\sum\limits_{i=1}^{N^{tot}_{walker}} w_i O(\bf{R}_i) \right) / \left( \sum \limits_{i=1}^{N^{tot}_{walker}} w_i \right). \end{aligned}

:math:`N^{tot}_{walker}` is the total number of walkers collected in the
entire simulation. Notice that :math:`N^{tot}_{walker}` is typically far
larger than the number of walkers held in memory at any given simulation
step. :math:`w_i` is the weight of walker :math:`i`.

In a VMC simulation, the weight of every walker is 1.0. Further, the
number of walkers is constant at each step. Therefore,
:eq:`eq242` simplifies to

.. math::
  :label: eq243

   \begin{aligned}
   E_{VMC}[O] = \frac{1}{N_{step}N_{walker}^{ensemble}} \sum_{s,e} O(\bf{R}_{s,e})\:.\end{aligned}

Each walker :math:`\bf{R}_{s,e}` is labeled by *step index* s and
*ensemble index* e.

In a DMC simulation, the weight of each walker is different and may
change from step to step. Further, the ensemble size varies from step to
step. Therefore, :eq:`eq242` simplifies to

.. math::
  :label: eq244

   \begin{aligned}
   E_{DMC}[O] = \frac{1}{N_{step}} \sum_{s} \left\{ \left(\sum_e w_{s,e} O(\bf{R}_{s,e})  \right) / \left( \sum \limits_{e} w_{s,e} \right)  \right\}\:.\end{aligned}

We will refer to the average in the :math:`\{\}` as *ensemble average*
and to the remaining averages as *block average*. The process of
calculating :math:`O(\mathbf{R})` is *evaluate*.

Class relations
~~~~~~~~~~~~~~~

A large number of classes are involved in the estimator collection process. They often have misleading class or method names. Check out the document gotchas in the following list:

#. ``EstimatorManager`` is an unused copy of ``EstimatorManagerBase``.
   ``EstimatorManagerBase`` is the class used in the QMC drivers. (PR
   #371 explains this.)

#. ``EstimatorManagerBase::Estimators`` is completely different from
   ``QMCDriver::Estimators``, which is subtly different from
   ``OperatorBase::Estimators``. The first is a list of pointers to
   ``ScalarEstimatorBase``. The second is the master estimator (one per
   MPI group). The third is the slave estimator that exists one per
   OpenMP thread.

#. ``QMCHamiltonian`` is NOT a parent class of ``OperatorBase``.
   Instead, ``QMCHamiltonian`` owns two lists of ``OperatorBase`` named
   ``H`` and ``auxH``.

#. ``QMCDriver::H`` is NOT the same as ``QMCHamiltonian::H``. The first
   is a pointer to a ``QMCHamiltonian``. ``QMCHamiltonian::H`` is a
   list.

#. ``EstimatorManager::stopBlock(std::vector)`` is completely different
   from ``EstimatorManager::`` ``stopBlock(RealType)``, which is the
   same as ``stopBlock(RealType, true)`` but that is subtly different
   from ``stopBlock(RealType, false)``. The first three methods are
   intended to be called by the master estimator, which exists one per
   MPI group. The last method is intended to be called by the slave
   estimator, which exists one per OpenMP thread.

Estimator output stages
~~~~~~~~~~~~~~~~~~~~~~~

Estimators take four conceptual stages to propagate to the output files: evaluate, load ensemble, unload ensemble, and collect. They are easier to understand in reverse order.

Collect stage
^^^^^^^^^^^^^

File output is performed by the master ``EstimatorManager`` owned by ``QMCDriver``. The first 8+ entries in ``EstimatorManagerBase::AverageCache`` will be written to ``scalar.dat``. The remaining entries in ``AverageCache`` will be written to ``stat.h5``. File writing is triggered by ``EstimatorManagerBase`` ``::collectBlockAverages`` inside ``EstimatorManagerBase::stopBlock``.

::

  // In EstimatorManagerBase.cpp::collectBlockAverages
    if(Archive)
    {
      *Archive << std::setw(10) << RecordCount;
      int maxobjs=std::min(BlockAverages.size(),max4ascii);
      for(int j=0; j<maxobjs; j++)
        *Archive << std::setw(FieldWidth) << AverageCache[j];
      for(int j=0; j<PropertyCache.size(); j++)
        *Archive << std::setw(FieldWidth) << PropertyCache[j];
      *Archive << std::endl;
      for(int o=0; o<h5desc.size(); ++o)
        h5desc[o]->write(AverageCache.data(),SquaredAverageCache.data());
      H5Fflush(h_file,H5F_SCOPE_LOCAL);
    }

``EstimatorManagerBase::collectBlockAverages`` is triggered from the master-thread estimator via either ``stopBlock(std::vector)`` or ``stopBlock(RealType, true)``. Notice that file writing is NOT triggered by the slave-thread estimator method ``stopBlock(RealType, false)``.

::

  // In EstimatorManagerBase.cpp
  void EstimatorManagerBase::stopBlock(RealType accept, bool collectall)
  {
    //take block averages and update properties per block
    PropertyCache[weightInd]=BlockWeight;
    PropertyCache[cpuInd] = MyTimer.elapsed();
    PropertyCache[acceptInd] = accept;
    for(int i=0; i<Estimators.size(); i++)
      Estimators[i]->takeBlockAverage(AverageCache.begin(),SquaredAverageCache.begin());
    if(Collectables)
    {
      Collectables->takeBlockAverage(AverageCache.begin(),SquaredAverageCache.begin());
    }
    if(collectall)
      collectBlockAverages(1);
  }

::

  // In ScalarEstimatorBase.h
  template<typename IT>
  inline void takeBlockAverage(IT first, IT first_sq)
  {
    first += FirstIndex;
    first_sq += FirstIndex;
    for(int i=0; i<scalars.size(); i++)
    {
      *first++ = scalars[i].mean();
      *first_sq++ = scalars[i].mean2();
      scalars_saved[i]=scalars[i]; //save current block
      scalars[i].clear();
    }
  }

At the collect stage, ``calarEstimatorBase::scalars`` must be populated with ensemble-averaged data. Two derived classes of ``ScalarEstimatorBase`` are crucial: ``LocalEnergyEstimator`` will carry ``Properties``, where as ``CollectablesEstimator`` will carry ``Collectables``.

Unload ensemble stage
^^^^^^^^^^^^^^^^^^^^^

``LocalEnergyEstimator::scalars`` are populated by
``ScalarEstimatorBase::accumulate``, whereas
``CollectablesEstimator::scalars`` are populated by
``CollectablesEstimator::`` ``accumulate_all``. Both accumulate methods
are triggered by ``EstimatorManagerBase::accumulate``. One confusing
aspect about the unload stage is that
``EstimatorManagerBase::accumulate`` has a master and a slave call
signature. A slave estimator such as ``QMCUpdateBase::Estimators``
should unload a subset of walkers. Thus, the slave estimator should call
``accumulate(W,it,it_end)``. However, the master estimator, such as
``SimpleFixedNodeBranch::myEstimator``, should unload data from the
entire walker ensemble. This is achieved by calling ``accumulate(W)``.

::

  void EstimatorManagerBase::accumulate(MCWalkerConfiguration& W)
  { // intended to be called by master estimator only
    BlockWeight += W.getActiveWalkers();
    RealType norm=1.0/W.getGlobalNumWalkers();
    for(int i=0; i< Estimators.size(); i++)
      Estimators[i]->accumulate(W,W.begin(),W.end(),norm);
    if(Collectables)//collectables are normalized by QMC drivers
      Collectables->accumulate_all(W.Collectables,1.0);
  }

::

  void EstimatorManagerBase::accumulate(MCWalkerConfiguration& W
   , MCWalkerConfiguration::iterator it
   , MCWalkerConfiguration::iterator it_end)
  { // intended to be called slaveEstimator only
    BlockWeight += it_end-it;
    RealType norm=1.0/W.getGlobalNumWalkers();
    for(int i=0; i< Estimators.size(); i++)
      Estimators[i]->accumulate(W,it,it_end,norm);
    if(Collectables)
      Collectables->accumulate_all(W.Collectables,1.0);
  }

::

  // In LocalEnergyEstimator.h
  inline void accumulate(const Walker_t& awalker, RealType wgt)
  { // ensemble average W.Properties
    // expect ePtr to be W.Properties; expect wgt = 1/GlobalNumberOfWalkers
    const RealType* restrict ePtr = awalker.getPropertyBase();
    RealType wwght= wgt* awalker.Weight;
    scalars[0](ePtr[WP::LOCALENERGY],wwght);
    scalars[1](ePtr[WP::LOCALENERGY]*ePtr[WP::LOCALENERGY],wwght);
    scalars[2](ePtr[LOCALPOTENTIAL],wwght);
    for(int target=3, source=FirstHamiltonian; target<scalars.size(); ++target, ++source)
      scalars[target](ePtr[source],wwght);
  }

::

  // In CollectablesEstimator.h
  inline void accumulate_all(const MCWalkerConfiguration::Buffer_t& data, RealType wgt)
  { // ensemble average W.Collectables
    // expect data to be W.Collectables; expect wgt = 1.0
    for(int i=0; i<data.size(); ++i)
      scalars[i](data[i], wgt);
  }

At the unload ensemble stage, the data structures ``Properties`` and ``Collectables`` must be populated by appropriately normalized values so that the ensemble average can be correctly taken. ``QMCDriver`` is responsible for the correct loading of data onto the walker ensemble.

Load ensemble stage
^^^^^^^^^^^^^^^^^^^

| ``Properties`` in the MC ensemble of walkers ``QMCDriver::W`` is
  populated by ``QMCHamiltonian``
| ``::saveProperties``. The master ``QMCHamiltonian::LocalEnergy``,
  ``::KineticEnergy``, and ``::Observables`` must be properly populated
  at the end of the evaluate stage.

::

  // In QMCHamiltonian.h
    template<class IT>
    inline
    void saveProperty(IT first)
    { // expect first to be W.Properties
      first[LOCALPOTENTIAL]= LocalEnergy-KineticEnergy;
      copy(Observables.begin(),Observables.end(),first+myIndex);
    }

``Collectables``'s load stage is combined with its evaluate stage.


Evaluate stage
^^^^^^^^^^^^^^

| The master ``QMCHamiltonian::Observables`` is populated by slave
  ``OperatorBase`` ``::setObservables``. However, the call signature
  must be ``OperatorBase::setObservables`` ``(QMCHamiltonian::``
| ``Observables)``. This call signature is enforced by
  ``QMCHamiltonian::evaluate`` and ``QMCHamiltonian::``
| ``auxHevaluate``.

::

  // In QMCHamiltonian.cpp
  QMCHamiltonian::Return_t
  QMCHamiltonian::evaluate(ParticleSet& P)
  {
    LocalEnergy = 0.0;
    for(int i=0; i<H.size(); ++i)
    {
      myTimers[i]->start();
      LocalEnergy += H[i]->evaluate(P);
      H[i]->setObservables(Observables);
  #if !defined(REMOVE_TRACEMANAGER)
      H[i]->collect_scalar_traces();
  #endif
      myTimers[i]->stop();
      H[i]->setParticlePropertyList(P.PropertyList,myIndex);
    }
    KineticEnergy=H[0]->Value;
    P.PropertyList[WP::LOCALENERGY]=LocalEnergy;
    P.PropertyList[LOCALPOTENTIAL]=LocalEnergy-KineticEnergy;
    // auxHevaluate(P);
    return LocalEnergy;
  }

::

  // In QMCHamiltonian.cpp
  void QMCHamiltonian::auxHevaluate(ParticleSet& P, Walker_t& ThisWalker)
  {
  #if !defined(REMOVE_TRACEMANAGER)
    collect_walker_traces(ThisWalker,P.current_step);
  #endif
    for(int i=0; i<auxH.size(); ++i)
    {
      auxH[i]->setHistories(ThisWalker);
      RealType sink = auxH[i]->evaluate(P);
      auxH[i]->setObservables(Observables);
  #if !defined(REMOVE_TRACEMANAGER)
      auxH[i]->collect_scalar_traces();
  #endif
      auxH[i]->setParticlePropertyList(P.PropertyList,myIndex);
    }
  }

Estimator use cases
~~~~~~~~~~~~~~~~~~~

VMCSingleOMP pseudo code
^^^^^^^^^^^^^^^^^^^^^^^^

::

  bool VMCSingleOMP::run()
  {
    masterEstimator->start(nBlocks);
    for (int ip=0; ip<NumThreads; ++ip)
      Movers[ip]->startRun(nBlocks,false);  // slaveEstimator->start(blocks, record)

    do // block
    {
      #pragma omp parallel
      {
        Movers[ip]->startBlock(nSteps);  // slaveEstimator->startBlock(steps)
        RealType cnorm = 1.0/static_cast<RealType>(wPerNode[ip+1]-wPerNode[ip]);
        do // step
        {
          wClones[ip]->resetCollectables();
          Movers[ip]->advanceWalkers(wit, wit_end, recompute);
          wClones[ip]->Collectables *= cnorm;
          Movers[ip]->accumulate(wit, wit_end);
        } // end step
        Movers[ip]->stopBlock(false);  // slaveEstimator->stopBlock(acc, false)
      } // end omp
      masterEstimator->stopBlock(estimatorClones);  // write files
    } // end block
    masterEstimator->stop(estimatorClones);
  }

DMCOMP  pseudo code
^^^^^^^^^^^^^^^^^^^

::

  bool DMCOMP::run()
  {
    masterEstimator->setCollectionMode(true);

    masterEstimator->start(nBlocks);
    for(int ip=0; ip<NumThreads; ip++)
      Movers[ip]->startRun(nBlocks,false);  // slaveEstimator->start(blocks, record)

    do // block
    {
      masterEstimator->startBlock(nSteps);
      for(int ip=0; ip<NumThreads; ip++)
        Movers[ip]->startBlock(nSteps);  // slaveEstimator->startBlock(steps)

      do // step
      {
        #pragma omp parallel
        {
        wClones[ip]->resetCollectables();
        // advanceWalkers
        } // end omp

        //branchEngine->branch
        { // In WalkerControlMPI.cpp::branch
        wgt_inv=WalkerController->NumContexts/WalkerController->EnsembleProperty.Weight;
        walkers.Collectables *= wgt_inv;
        slaveEstimator->accumulate(walkers);
        }
        masterEstimator->stopBlock(acc)  // write files
      }  // end for step
    }  // end for block

    masterEstimator->stop();
  }

Summary
~~~~~~~

Two ensemble-level data structures, ``ParticleSet::Properties`` and
``::Collectables``, serve as intermediaries between evaluate classes and
output classes to ``scalar.dat`` and ``stat.h5``. ``Properties`` appears
in both ``scalar.dat`` and ``stat.h5``, whereas ``Collectables`` appears
only in ``stat.h5``. ``Properties`` is overwritten by
``QMCHamiltonian::Observables`` at the end of each step.
``QMCHamiltonian::Observables`` is filled upon call to
``QMCHamiltonian::evaluate`` and ``::auxHevaluate``. ``Collectables`` is
zeroed at the beginning of each step and accumulated upon call to
``::auxHevaluate``.

| Data are output to ``scalar.dat`` in four stages: evaluate, load,
  unload, and collect. In the evaluate stage,
  ``QMCHamiltonian::Observables`` is populated by a list of
  ``OperatorBase``. In the load stage, ``QMCHamiltonian::Observables``
  is transferred to ``Properties`` by ``QMCDriver``. In the unload stage,
  ``Properties`` is copied to ``LocalEnergyEstimator::scalars``. In the
  collect stage, ``LocalEnergyEstimator::scalars`` is block-averaged to
  ``EstimatorManagerBase``
| ``::AverageCache`` and dumped to file. For ``Collectables``, the
  evaluate and load stages are combined in a call to
  ``QMCHamiltonian::auxHevaluate``. In the unload stage,
  ``Collectables`` is copied to ``CollectablesEstimator::scalars``. In
  the collect stage, ``CollectablesEstimator``
| ``::scalars`` is block-averaged to
  ``EstimatorManagerBase::AverageCache`` and dumped to file.

Appendix: dmc.dat
~~~~~~~~~~~~~~~~~

| There is an additional data structure,
  ``ParticleSet::EnsembleProperty``, that is managed by
  ``WalkerControlBase::EnsembleProperty`` and directly dumped to
  ``dmc.dat`` via its own averaging procedure. ``dmc.dat`` is written by
  ``WalkerControlBase::measureProperties``, which is called by
  ``WalkerControlBase::branch``, which is called by
  ``SimpleFixedNodeBranch``
| ``::branch``, for example.

.. bibliography:: /bibs/developing.bib
