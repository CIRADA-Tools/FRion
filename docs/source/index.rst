.. FRion documentation master file, created by
   sphinx-quickstart on Fri Apr 16 15:03:12 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _homepage:

Welcome to FRion's documentation!
=================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   :hidden:

   predict
   correct

FRion is a Python3 package for the prediction and correction of ionospheric
Faraday rotation, which can be useful in certain radio astronomy applications.
FRion focuses on time-averaged effects of the ionosphere, for cases where data
gets time-averaged before an ionospheric correction can be applied.

Users interested in time-dependent corrections should instead look at other 
packages, such as `RMextract <https://github.com/lofar-astron/RMextract/>`_ 
(which is used internally by FRion), `ionFR <https://github.com/csobey/ionFR>`_,
or `ALBUS <https://github.com/twillis449/ALBUS_ionosphere>`_.

A mathematical derivation of how the time-independent ionospheric Faraday 
rotation correction is defined, along with some remarks on its use, 
can be found :download:`here <./Ionospheric_Correction.pdf>`.

The package consists of two parts:

* :ref:`predict`: Functions for predicting the ionospheric Faraday rotation and computing the time-average.
* :ref:`correct`: Functions for correcting Stokes Q and U data cubes for the effects of ionospheric Faraday rotation.

Each part can be imported into Python scripts, or the basic functionality can
be accessed through the following terminal commands:

- ``frion_predict`` 
    Runs the prediction script, given user-supplied observation time, location, 
    direction, and frequencies. Can read this information from a 
    FITS header.
- ``frion_correct`` 
    Runs the correction script, applying a correction
    to a pair of Stokes Q and U cubes to removed the predicted ionospheric 
    modulation and depolarization.

Use the ``-h`` flag to get detailed usage instructions for each.


Usage
------------





Installation
------------

FRion will be released on PyPi once it is complete and tested, which will
enable easy installation through pip.

In the mean time, it can be installed by downloading the code from 
`this link <https://github.com/Cameron-Van-Eck/FRion/archive/refs/heads/main.zip>`_, 
unzipping, moving the code directory somewhere convenient, 
going into the code directory, then running ``pip install -e .``.
This will install the package to the Python packages directory.

RMextract must be installed seperately it. It is now available through pip, 
using ``pip install RMextract``, but this will try to install casacore as a dependency. 
casacore can be difficult to install on some systems, so if this causes a problem
you can install it without casacore by using ``pip install --no-deps RMextract``.
casacore is not required: if casacore is not installed, RMextract will use the ephem package instead, which installs
automatically with FRion.


It should then
be importable using the statements ``import FRion.predict as predict`` and
``import FRion.correct as correct``, or runable on the terminal with the commands
``frion_predict`` and ``frion_correct``.









Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
