.. _reading:

Recommended Reading
===================

The sections below contain information, or at least links to information,
that should be helpful for anyone who wants to use Nexus, but who
is not an expert in one of the following areas: installing python and related
modules, installing PWSCF and QMCPACK, the Python programming language, and
the theory and practice of Quantum Monte Carlo.

.. _install-python:

Helpful Links for Installing Python Modules
-------------------------------------------

Python itself
   |
   | Download: http://www.python.org/download/
   | Be sure to get Python 2.x, not 3.x.

Numpy and Scipy
   |
   | Download and installation: http://www.scipy.org/Installing_SciPy.

Matplotlib
   |
   | Download: http://matplotlib.org/downloads.html
   | Installation: http://matplotlib.org/users/installing.html

H5py
   |
   | Download and installation: http://www.h5py.org/

.. _install-code:

Helpful Links for Installing Electronic Structure Codes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PWSCF: pw.x, pw2qmcpack.x, pw2casino.x
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Download: http://www.quantum-espresso.org/ ,
  http://qmcpack.org/downloads/
| Installation instructions: See the ``README`` file in
  ``qmcpack/external_codes/quantum_espresso`` for instructions on
  patching Quantum Espresso version 5.1.

QMCPACK: qmcpack, qmcpack_complex, convert4qmc, ppconvert, sqd
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Download: http://qmcpack.org/downloads/
| Install: http://docs.qmcpack.org/ug/a00002.html

Wfconvert: wfconvert
~~~~~~~~~~~~~~~~~~~~

| Download: svn co http://qmctools.googlecode.com/svn/trunk/wfconvert
| See also: https://code.google.com/p/qmctools/

VASP
~~~~

| Download: Proprietary, see https://www.vasp.at/
| Install: http://cms.mpi.univie.ac.at/vasp/vasp/vasp.html

GAMESS
~~~~~~

Download: http://www.msg.ameslab.gov/gamess/License_Agreement.html
Install: See the “readme.unix” file in the GAMESS source distribution
(gamess/machines/readme.unix).

.. _learn-python:

Brushing up on Python
---------------------

Python
~~~~~~

Python is a flexible, multi-paradigm, interpreted programming language with
powerful intrinsic datatypes and a large library of modules that greatly expand
its functionality.  A good way to learn the language is through the extensive
Documentation provided on the python.org website. If you have never worked with
Python before, be sure to go through the Tutorial. To learn more about the
intrinsic data types and standard libraries look at Library Reference.
A very short introduction to Python is in :ref:`python-basics`.

================= ============================================
Documentation     http://docs.python.org/2/
================= ============================================
Tutorial          http://docs.python.org/2/tutorial/index.html
Library Reference http://docs.python.org/2/library/index.html
================= ============================================

NumPy
~~~~~

Other than the Python Standard Library, the main library/module Nexus
makes heavy use of is NumPy.  NumPy provides a convenient and fairly
fast implementation of multi-dimensional arrays and related functions, much like
MATLAB.  If you want to learn about NumPy arrays, the NumPy
Tutorial is recommended.  For more detailed information, see the NumPy User Guide
and the NumPy Reference Manual. If MATLAB is one of your native languages, check out
NumPy for MATLAB Users.

========== ====================================================
Tutorial   http://www.scipy.org/Tentative_NumPy_Tutorial
========== ====================================================
User Guide http://docs.scipy.org/doc/numpy/user/index.html#user
Reference  http://docs.scipy.org/doc/numpy/reference/
MATLAB     http://www.scipy.org/NumPy_for_Matlab_Users
========== ====================================================

Matplotlib
~~~~~~~~~~

Plotting in Nexus is currently handled by Matplotlib. If you want to
learn more about plotting with Matplotlib, the Pyplot Tutorial is a good
place to start. More detailed information is in the User’s Guide.
Sometimes Examples provide the fastest way to learn.

============ ================================================
Tutorial     http://matplotlib.org/users/pyplot_tutorial.html
============ ================================================
User’s Guide http://matplotlib.org/users/index.html
Examples     http://matplotlib.org/examples/
============ ================================================

Scipy and H5Py
~~~~~~~~~~~~~~

Nexus also occasionally uses functionality from SciPy and H5Py. Learning
more about them is unlikely to help you interact with Nexus. However,
they are quite valuable on their own. SciPy provides access to special
functions, numerical integration, optimization, interpolation, fourier
transforms, eigenvalue solvers, and statistical analysis. To get an
overview, try the SciPy Tutorial. More detailed material is found in the
Scipy Reference. H5Py provides a NumPy-like interface to HDF5 data
files, which QMCPACK creates. To learn more about interacting with HDF5
files through H5Py, try the Quick Start Guide. For more information, see
the General Documentation.

+-------------------+--------------------------------------------------------------+
| SciPy Tutorial    | http://docs.scpy.org/doc/scipy/reference/tutorial/index.html |
|                   |                                                              |
+-------------------+--------------------------------------------------------------+
| SciPy Reference   | http://docs.scipy.org/doc/scipy/reference/                   |
+-------------------+--------------------------------------------------------------+
| H5Py Quick Guide  | http://docs.h5py.org/en/stable/quick.html                    |
+-------------------+--------------------------------------------------------------+
| H5Py General Docs | http://docs.h5py.org/en/stable/                              |
+-------------------+--------------------------------------------------------------+

.. _learn-qmc:

Quantum Monte Carlo: Theory and Practice
----------------------------------------

Currently, review articles may be the best way to get an overview of
Quantum Monte Carlo methods and practices. The review article by
Foulkes, *et al.* from 2001 remains quite relevant and is lucidly
written. Other review articles also provide a broader perspective on
QMC, including more recent developments. Another resource that can be
useful for newcomers (and needs to be updated) is the QMC Wiki. If you
are aware of resources that fill a gap in the information presented here
(almost a certainty), please contact the developer at krogeljt@ornl.gov
to add your contribution.

===================   ======================================================
QMC Review Articles
===================   ======================================================
Foulkes, 2001         http://rmp.aps.org/abstract/RMP/v73/i1/p33_1
Bajdich, 2009         http://www.physics.sk/aps/pub.php?y=2009&pub=aps-09-02
Needs, 2010           http://iopscience.iop.org/0953-8984/22/2/023201/
Kolorenc, 2011        http://iopscience.iop.org/0034-4885/74/2/026502/
===================   ======================================================

======================  =======================================================
Online Resources
======================  =======================================================
QMCWiki &               http://www.qmcwiki.org
QMC Summer School 2012  http://www.mcc.uiuc.edu/summerschool/2012/program.html
======================  =======================================================
