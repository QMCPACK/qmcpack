.. _python-basics:

Basic Python Constructs
=======================

Basic Python data types (``int``, ``float``, ``str``, ``tuple``,
``list``, ``array``, ``dict``, ``obj``) and programming constructs
(``if`` statements, ``for`` loops, functions w/ keyword arguments) are
briefly overviewed below. All examples can be executed interactively in
Python. To do this, type “``python``” at the command line and paste any
of the shaded text below at the “``>>>``” prompt. For more information
about effective use of Python, consult the detailed online
documentation: https://docs.python.org/2/.

Intrinsic types: ``int, float, str``
------------------------------------

.. code:: rest

  #this is a comment
  i=5                     # integer
  f=3.6                   # float
  s='quantum/monte/carlo' # string
  n=None                  # represents "nothing"

  f+=1.4                  # add-assign (-,*,/ also): 5.0
  2**3                    # raise to a power: 8
  str(i)                  # int to string: '5'
  s+'/simulations'        # joining strings: 'quantum/monte/carlo/simulations'
  'i={0}'.format(i)       # format string: 'i=5'

Container types: ``tuple, list, array, dict, obj``
--------------------------------------------------

.. code:: rest

  from numpy import array  # get array from numpy module
  from generic import obj  # get obj from generic module

  t=('A',42,56,123.0)     # tuple

  l=['B',3.14,196]        # list

  a=array([1,2,3])        # array

  d={'a':5,'b':6}         # dict

  o=obj(a=5,b=6)          # obj

                          # printing
  print t                 #  ('A', 42, 56, 123.0)
  print l                 #  ['B', 3.1400000000000001, 196]
  print a                 #  [1 2 3]
  print d                 #  {'a': 5, 'b': 6}
  print o                 #    a               = 5
                          #    b               = 6

  len(t),len(l),len(a),len(d),len(o) #number of elements: (4, 3, 3, 2, 2)

  t[0],l[0],a[0],d['a'],o.a  #element access: ('A', 'B', 1, 5, 5)

  s = array([0,1,2,3,4])  # slices: works for tuple, list, array
  s[:]                    #   array([0, 1, 2, 3, 4])
  s[2:]                   #   array([2, 3, 4])
  s[:2]                   #   array([0, 1])
  s[1:4]                  #   array([1, 2, 3])
  s[0:5:2]                #   array([0, 2, 4])

                          # list operations
  l2 = list(l)            #   make independent copy
  l.append(4)             #   add new element: ['B', 3.14, 196, 4]
  l+[5,6,7]               #   addition: ['B', 3.14, 196, 4, 5, 6, 7]
  3*[0,1]                 #   multiplication:  [0, 1, 0, 1, 0, 1]

  b=array([5,6,7])        # array operations
  a2 = a.copy()           #   make independent copy
  a+b                     #   addition: array([ 6, 8, 10])
  a+3                     #   addition: array([ 4, 5, 6])
  a*b                     #   multiplication: array([ 5, 12, 21])
  3*a                     #   multiplication: array([3, 6, 9])

                          # dict/obj operations
  d2 = d.copy()           #   make independent copy
  d['c'] = 7              #   add/assign element
  d.keys()                #   get element names: ['a', 'c', 'b']
  d.values()              #   get element values: [5, 7, 6]

                          # obj-specific operations
  o.c = 7                 #   add/assign element
  o.set(c=7,d=8)          #   add/assign multiple elements

An important feature of Python to be aware of is that assignment is most
often by reference, *i.e.* new values are not always created. This point
is illustrated below with an ``obj`` instance, but it also holds for
``list``, ``array``, ``dict``, and others.

.. code:: rest

  >>> o = obj(a=5,b=6)
  >>>
  >>> p=o
  >>>
  >>> p.a=7
  >>>
  >>> print o
    a               = 7
    b               = 6

  >>> q=o.copy()
  >>>
  >>> q.a=9
  >>>
  >>> print o
    a               = 7
    b               = 6

Here ``p`` is just another name for ``o``, while ``q`` is a fully
independent copy of it.

Conditional Statements: ``if/elif/else``
----------------------------------------

.. code:: rest

  a = 5
  if a is None:
      print 'a is None'
  elif a==4:
      print 'a is 4'
  elif a<=6 and a>2:
      print 'a is in the range (2,6]'
  elif a<-1 or a>26:
      print 'a is not in the range [-1,26]'
  elif a!=10:
      print 'a is not 10'
  else:
      print 'a is 10'
  #end if

The “``#end if``” is not part of Python syntax, but you will see text
like this throughout the Project Suite for clear encapsulation.

Iteration: ``for``
------------------

.. code:: rest

  from generic import obj

  l = [1,2,3]
  m = [4,5,6]
  s = 0
  for i in range(len(l)):  # loop over list indices
      s += l[i] + m[i]
  #end for

  print s                  # s is 21

  s = 0
  for v in l:              # loop over list elements
      s += v
  #end for

  print s                  # s is 6

  o = obj(a=5,b=6)
  s = 0
  for v in o:              # loop over obj elements
      s += v
  #end for

  print s                  # s is 11

  d = {'a':5,'b':4}
  for n,v in o.items():# loop over name/value pairs in obj
      d[n] += v
  #end for

  print d                  # d is {'a': 10, 'b': 10}

Functions: ``def``, argument syntax
-----------------------------------

.. code:: rest

  def f(a,b,c=5):          # basic function, c has a default value
      print a,b,c
  #end def f

  f(1,b=2)                 # prints: 1 2 5


  def f(*args,**kwargs):   # general function, returns nothing
      print args           #     args: tuple of positional arguments
      print kwargs         #   kwargs: dict of keyword arguments
  #end def f

  f('s',(1,2),a=3,b='t')   # 2 pos., 2 kw. args, prints:
                           #   ('s', (1, 2))
                           #   {'a': 3, 'b': 't'}

  l = [0,1,2]
  f(*l,a=6)                # pos. args from list, 1 kw. arg, prints:
                           #   (0, 1, 2)
                           #   {'a': 6}
  o = obj(a=5,b=6)
  f(*l,**o)                # pos./kw. args from list/obj, prints:
                           #   (0, 1, 2)
                           #   {'a': 5, 'b': 6}

  f(                       # indented kw. args, prints
      blocks   = 200,      #   ()
      steps    = 10,       #   {'steps': 10, 'blocks': 200, 'timestep': 0.01}
      timestep = 0.01
      )

  o = obj(                 # obj w/ indented kw. args
      blocks   = 100,
      steps    =  5,
      timestep = 0.02
      )

  f(**o)                   # kw. args from obj, prints:
                           #   ()
                           #   {'timestep': 0.02, 'blocks': 100, 'steps': 5}
