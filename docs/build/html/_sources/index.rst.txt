.. DICToolbox documentation master file, created by
   sphinx-quickstart on Tue Nov 20 10:23:06 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to DIC Toolbox's documentation!
========================================

DIC Toolbox is a python tools created to work on 2D strain field map such as the one obtained using 7D software from Pierre Vacher (Annecy). It is an open source code under CC-BY-CC licence (https://creativecommons.org/licenses/by/2.0/fr/deed.en). There is no guarantee if you are using it. It has been tested with python 3.7


.. toctree::
   :maxdepth: 1
   :caption: Documentation:


Installation
============

For simple user
***************

If you want just to use the toolbox without doing any devellopement

.. code:: bash

	pip install git+https://github.com/ThomasChauve/DICToolbox


For develloper
**************

If you want to access the code and be able to do some modification

.. code:: bash

    git clone https://github.com/ThomasChauve/DICToolbox
    cd DICToolbox/
    pip install .

Then you will find all the package in python using

.. code:: python

    import DICToolbox

Uninstall
*********

.. code:: bash
    
    pip uninstall DICToolbox


.. toctree::
    :maxdepth: 1
    :numbered:
    :caption: Documentation

    Documentation/Documentation
    CLASS/dic

.. toctree::
    :maxdepth: 1
    :numbered:
    :caption: CLASS

    CLASS/image2d
    CLASS/sTM
    CLASS/utils


Contact
=======
:Author: Thomas Chauve
:Contact: thomas.chauve@univ-grenoble-alpes.fr

:organization: UiO
:status: This is a "work in progress"
:version: 2.0



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
