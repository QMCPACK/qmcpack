.. _contrib:

Contributing to the Manual
==========================

This section briefly describes how to contribute to the manual and is primarily "by developers, for developers."   This section should iterate until a consistent view on style/contents is reached.

**Desirable:**

-  Use the following table templates when describing XML input.

-  Unicode rules

   -  Do not use characters for which well-established idioms
      exists, especially dashes, quotes, and apostrophes.

   -  Use math mode markup instead of unicode characters for equations.

   -  Be cautious of WYSIWYG word processors; cutting and pasting can
      pickup characters promoted to unicode by the program.

   -  Take a look at your text multibyte expanded; that is open it in
      (emacs and ‘esc-x toggle-enable-multibyte-characters’)—see any
      unicode you did not intend?

-  Newly added entries to a Bib file should be as complete as possible.
   Use a tool such as JabRef or Zotero that can automate creation of
   these entries from just a DOI.

**Forbidden:**

-  Including images instead of text tables.

-  Saving files in encodings other than UTF8. Some may
   report being ASCII encoded since they contain no unicode characters.

**Missing sections (these are opinions, not decided priorities):**

-  Description of XML input in general. Discuss XML format, use of
   attributes and ``<parameter/>`` s in general, case sensitivity
   (input is generally case sensitive), and behavior of when
   unrecognized XML elements are encountered (they are generally ignored
   without notification).

-  Overview of the input file in general, broad structure, and at least
   one full example that works in isolation.

**Information currently missing for a complete reference
specification:**

-  Noting how many instances of each child element are allowed.
   Examples: ``simulation``–1 only, ``method``–1 or more, ``jastrow``–0
   or more.

Table templates follow for describing XML elements in reference fashion.
A number of examples can be found in, for example,
:ref:`hamiltobs`. Preliminary style is (please weigh in with
opinions): typewriter text (``\texttt\{}``) for XML elements, attributes, and
parameter names; normal text for literal information in the datatype,
values, and default columns; bold (``\textbf{}``) text if an attribute or parameter
must take on a particular value (values column); italics (``\textit{}``) for
descriptive (nonliteral) information in the values column (e.g.,
*anything*, *non-zero*); and required/optional attributes or parameters
noted by ``some_attr`` :math:`^r`/``some_attr`` :math:`^r` superscripts. Valid datatypes
are text, integer, real, Boolean, and arrays of each. Fixed length
arrays can be noted, for example, by “real array(3).”

Template for a generic XML element:

``generic`` element:

+------------------+----------------------------------+
| parent elements: | ``parent1`` ``parent2``          |
+------------------+----------------------------------+
| child elements:  | ``child1`` ``child2`` ``child3`` |
+------------------+----------------------------------+

  attributes:

  +---------------------------------+---------------+------------+-------------+-----------------+
  | **Name**                        | **Datatype**  | **Values** | **Default** | **Description** |
  +=================================+===============+============+=============+=================+
  | ``attr1``\ :math:`^r`           | text          |            |             |                 |
  +---------------------------------+---------------+------------+-------------+-----------------+
  | ``attr2``\ :math:`^r`           | integer       |            |             |                 |
  +---------------------------------+---------------+------------+-------------+-----------------+
  | ``attr3``\ :math:`^r`           | real          |            |             |                 |
  +---------------------------------+---------------+------------+-------------+-----------------+
  | ``attr4``\ :math:`^r`           | boolean       |            |             |                 |
  +---------------------------------+---------------+------------+-------------+-----------------+
  | ``attr5``\ :math:`^r`           | text array    |            |             |                 |
  +---------------------------------+---------------+------------+-------------+-----------------+
  | ``attr6``\ :math:`^r`           | integer array |            |             |                 |
  +---------------------------------+---------------+------------+-------------+-----------------+
  | ``attr7``\ :math:`^r`           | real array    |            |             |                 |
  +---------------------------------+---------------+------------+-------------+-----------------+
  | ``attr8``\ :math:`^r`           | boolean array |            |             |                 |
  +---------------------------------+---------------+------------+-------------+-----------------+

  parameters:

  +----------------------------------+---------------+------------+-------------+-----------------+
  | **Name**                         | **Datatype**  | **Values** | **Default** | **Description** |
  +==================================+===============+============+=============+=================+
  | ``param1``\ :math:`^r`           | text          |            |             |                 |
  +----------------------------------+---------------+------------+-------------+-----------------+
  | ``param2``\ :math:`^r`           | integer       |            |             |                 |
  +----------------------------------+---------------+------------+-------------+-----------------+
  | ``param3``\ :math:`^r`           | real          |            |             |                 |
  +----------------------------------+---------------+------------+-------------+-----------------+
  | ``param4``\ :math:`^r`           | boolean       |            |             |                 |
  +----------------------------------+---------------+------------+-------------+-----------------+
  | ``param5``\ :math:`^r`           | text array    |            |             |                 |
  +----------------------------------+---------------+------------+-------------+-----------------+
  | ``param6``\ :math:`^r`           | integer array |            |             |                 |
  +----------------------------------+---------------+------------+-------------+-----------------+
  | ``param7``\ :math:`^r`           | real array    |            |             |                 |
  +----------------------------------+---------------+------------+-------------+-----------------+
  | ``param8``\ :math:`^r`           | boolean array |            |             |                 |
  +----------------------------------+---------------+------------+-------------+-----------------+

  body text: Long form description of body text format

“Factory” elements are XML elements that share a tag but whose contents
change based on the value an attribute, or sometimes multiple
attributes, take. The attribute(s) that determines the allowed content
is subsequently referred to as the “type selector” (e.g., for
``<estimator/>`` elements, the type selector is usually the ``type``
attribute). These types of elements are frequently encountered as they
correspond (sometimes loosely, sometimes literally) to polymorphic
classes in QMCPACK that are built in “factories.” This name is true to
the underlying code but may be obscure to the general user (is there a
better name to retain the general meaning?).

The following template should be provided each time a new “factory” type
is encountered (such as ``<estimator/>``). The table lists all types of
possible elements (see “type options” in the template) and any
attributes that are common to all possible related elements. Specific
“derived” elements are then described one at a time with the previous
template, noting the type selector in addition to the XML tag (e.g.,
“``estimator type=density`` element”).

Template for shared information about “factory” elements.

``generic`` factory element:

+------------------+----------------------------------+
| parent elements: | ``parent1`` ``parent2``          |
+------------------+----------------------------------+
| child elements:  | ``child1`` ``child2`` ``child3`` |
+------------------+----------------------------------+
| type selector    | ``some`` attribute               |
+------------------+----------------------------------+
| type options     | Selection 1                      |
+------------------+----------------------------------+
|                  | Selection 2                      |
+------------------+----------------------------------+
|                  | Selection 3                      |
+------------------+----------------------------------+
|                  | ...                              |
+------------------+----------------------------------+

  shared attributes:

  +---------------------------------+---------------+------------+-------------+-----------------+
  | **Name**                        | **Datatype**  | **Values** | **Default** | **Description** |
  +=================================+===============+============+=============+=================+
  | ``attr1``\ :math:`^r`           | text          |            |             |                 |
  +---------------------------------+---------------+------------+-------------+-----------------+
  | ``attr2``\ :math:`^r`           | integer       |            |             |                 |
  +---------------------------------+---------------+------------+-------------+-----------------+
  | ...                             |               |            |             |                 |
  +---------------------------------+---------------+------------+-------------+-----------------+
