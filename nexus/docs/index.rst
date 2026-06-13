.. NEXUS documentation documentation master file, created by
   sphinx-quickstart on Mon Jul 13 18:25:36 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

:html_theme.sidebar_secondary.remove:

####################
The Nexus User Guide
####################

:Release: |release|
:Date: |today|

**Helpful Links**:
`GitHub <https://github.com/QMCPACK/qmcpack>`_ |
`Issue Tracking <https://github.com/QMCPACK/qmcpack/issues>`_

The Nexus Python package provides both basic and advanced tooling for running calculations on supercomputing clusters for several density functional theory (DFT) and wavefunction theory (WFT) programs, as well as for driving QMCPACK calculations.

.. toctree::
   :hidden:

   user_guide/index
   api_docs/index
   extra_reading/index
   code-style


.. grid:: 2
    :gutter: 2 3 4 4

    .. grid-item-card::
        :text-align: center

        **Getting Started**
        ^^^

        New user guide for Nexus

        +++

        .. button-ref:: user_guide/index
            :ref-type: doc
            :expand:
            :color: secondary
            :click-parent:

            To the User Guide

    .. grid-item-card::
        :text-align: center

        **API Documentation**
        ^^^

        API Documentation for Nexus

        +++

        .. button-ref:: api_docs/index
            :ref-type: doc
            :expand:
            :color: secondary
            :click-parent:

            To the API Docs


.. grid:: 2
    :gutter: 2 3 4 4

    .. grid-item-card::
        :text-align: center

        **Extra Reading**
        ^^^

        Recommended reading for those who are new to Python or QMC practice.

        +++

        .. button-ref:: extra_reading/index
            :ref-type: doc
            :expand:
            :color: secondary
            :click-parent:

            To the Extra Reading

    .. grid-item-card::
        :text-align: center

        **Contributing to Nexus**
        ^^^

        Guidelines for those who wish to contribute to Nexus

        +++

        .. button-ref:: code-style
            :ref-type: doc
            :expand:
            :color: secondary
            :click-parent:

            To the Contributor Guide


Supported Code APIs
===================

.. card::
    :margin: 3

    .. button-ref:: qmcpack-api
        :ref-type: ref
        :color: primary
        :expand:
        :click-parent:

        QMCPACK

.. grid:: 2
    :gutter: 1

    .. grid-item::

        .. grid:: 1 1 1 1
            :gutter: 3

            .. grid-item-card:: 

                .. button-ref:: espresso-api
                    :ref-type: ref
                    :color: primary
                    :expand:
                    :click-parent:

                    Quantum ESPRESSO

            .. grid-item-card:: 

                .. button-ref:: vasp-api
                    :ref-type: ref
                    :color: primary
                    :expand:
                    :click-parent:

                    VASP

            .. grid-item-card:: 

                .. button-ref:: quantum_package-api
                    :ref-type: ref
                    :color: primary
                    :expand:
                    :click-parent:

                    Quantum Package

    .. grid-item::

        .. grid:: 1 1 1 1
            :gutter: 3

            .. grid-item-card:: 

                .. button-ref:: pyscf-api
                    :ref-type: ref
                    :color: primary
                    :expand:
                    :click-parent:

                    PySCF

            .. grid-item-card:: 

                .. button-ref:: gamess-api
                    :ref-type: ref
                    :color: primary
                    :expand:
                    :click-parent:

                    GAMESS

            .. grid-item-card:: 
                
                .. button-ref:: rmg-api
                    :ref-type: ref
                    :color: primary
                    :expand:
                    :click-parent:

                    RMG

.. `homepage <https://qmcpack.org/>`__
.. `homepage <https://www.quantum-espresso.org/>`__
.. `homepage <https://pyscf.org/>`__
.. `homepage <https://www.vasp.at/>`__
.. `homepage <https://www.msg.chem.iastate.edu/gamess/documentation.html>`__
.. `homepage <https://quantumpackage.github.io/qp2/>`__
.. `homepage <http://www.rmgdft.org/>`__