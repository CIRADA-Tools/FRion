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
gets time-averaged before an ionospheric correction can be applied, but also
has functions to produce time-series of the ionospheric Faraday rotation.

This package uses `RMExtract <https://github.com/lofar-astron/RMextract/>`_ 
for the underlying ionospheric calculations. Users interested in alternative
ionospheric Faraday rotation packages can look at
`ionFR <https://github.com/csobey/ionFR>`_,
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
    Runs the time-averaged prediction script, given user-supplied observation 
    time, location, direction, and frequencies. Can read this information from 
    a FITS header.
- ``frion_correct`` 
    Runs the correction script, applying a correction (generated by the predict 
    functions) to a pair of Stokes Q and U cubes to removed the predicted 
    ionospheric modulation and depolarization.

Use the ``-h`` flag to get detailed usage instructions for each.

**Notice:** The use of FRion has changed somewhat significantly in version 1.1.
See the :ref:`notes below<ver1.1>` for more information about the changes. The most
important change is that FRion now requires an account with CDDIS and a corresponding
.netrc file.


Installation
------------

FRion has been released on PyPi, so it can be installed with pip as ``pip install FRion``.

Alternatively, it can be installed by downloading the code from 
`this link <https://github.com/Cameron-Van-Eck/FRion/archive/refs/heads/main.zip>`_, 
unzipping, moving the code directory somewhere convenient, 
going into the code directory, then running ``pip install -e .``.
This will install the package to the Python packages directory.

RMExtract must be installed seperately. It is now available through pip, 
using ``pip install RMextract``, but this will try to install casacore as a dependency. 
casacore can be difficult to install on some systems, so my reccomenadation is
to install it without casacore by using ``pip install --no-deps RMextract``.
casacore is not required for RMExtract: if casacore is not installed, RMExtract will use the pyephem package instead, which installs
automatically with FRion.


It should then
be importable using the statements ``import FRion.predict as predict`` and
``import FRion.correct as correct``, or runable on the terminal with the commands
``frion_predict``, ``frion_timeseries``, and ``frion_correct``.

Note that to use the default mode (including the command line tools) requires
an account with CDDIS and a corresponding .netrc file in order to download TEC 
data; see below for details.

In some cases users may encounter an error 
``RuntimeError: Cannot convert due to missing frame information``.
This occurs when RMextract finds casacore and tries to use it, but is missing
the casadata module. This can be solved by installing casadata
(``pip install --index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casadata``),
and updating the ``.casarc`` file to point to the updated install.
Removing casacore 
(``pip uninstall python-casacore``) solves the issue by forcing RMextract to
rely on the ephem module instead.


Usage
------------
To generate ionospheric Faraday rotation predictions as a function of time, 
the user can use the predict timeseries functions. These can be used from
the command line using ``frion_timeseries``. This requires the user to supply the
start and end times of the observation, the location of the telescope, and the 
sky coordinates. These values can be supplied from a FITS or PSRFITS file if
the correct keywords are in the file's header. The user can choose to save the
values to a file and/or generate a plot.

These features can be accessed in a Python script using the functions 
available in the :ref:`predict` module, starting with :py:func:`FRion.predict.get_RM()`.

Generating time-averaged predictions can be done similarly from the command line
using the ``frion_predict`` command or the functions in the :ref:`predict` module.
These predictions require the same information, plus the frequencies of each 
channel.

Warnings about reduced accuracy when using PyEphem can be safely ignored 
(accuracy is approximately 1 arcsecond, and ionospheric data is so coarse
that this has no effect on results).

Stokes Q and U FITS cubes can be corrected for the time-averaged Faraday rotation,
using the :ref:`correct` module. Tools for time-dependent corrections to 
different data types (pulsar observation, visibilities, etc) are outside the 
scope of this package.

The correction functions can be used from the command line using ``frion_correct``,
or within a script using :py:func:`FRion.correct.apply_correction_to_files()`.
The correction relies on the predictions generated by the predict module, so
the user must run the predict tools first and save the ionospheric modulation 
to a text file; this text file is used as an input by the correct tools.

Note that the correct tools will create a new pair of Stokes Q and U FITS cubes,
with the same size as the input cubes. The user must ensure that sufficient disk
space is available, otherwise this step will fail.

The default correct functions require holding the full Q and U cubes in RAM
while processing, which may overflow some systems (requiring the use of much 
slower virtual memory). A large-file version has also been developed to
reduce this memory footprint, and can be enabled in the command-line tool by
setting the ``-L`` flag or in scripts by using the 
:py:func:`FRion.correct.apply_correction_large_cube()` function.
 
 


.. _ver1.1:

Version 1.1
--------------------
Versions of FRion prior to 1.1 used the default settings in RMextract for
downloading ionospheric TEC data, which was to download the CODE global TEC
maps from CODE. Unfortunately, CODE stopped producing new TEC maps in January
2023, so it can no longer be used for observations made after that time.

The best alternative source of TEC data was the `NASA CDDIS archive`_,
but this requires a (free) account to access the data. At time of writing 
(Aug 2023), RMextract does not support this. A workaround was developed by
Anna Ordog and Art Davydov at DRAO that allows RMextract to download from CDDIS,
but it is not clear when this will be deployed into the official release.

As a workaround, I have modified FRion to download the TEC data itself before
invoking RMextract. This bypasses the need for RMextract to do any downloading,
which is currently not supported for CDDIS data. This can be bypassed by advanced
users, by setting the keyword ``pre_download = False`` in the predict functions.
Also, I have changed the default TEC source to be the JPL global TEC maps. This
was selected based on the results of `Porayko et al. 2019`_
who found that the JPL maps produced the lowest residual errors in LOFAR observations.

Downloading data from CDDIS requires an account. Information about applying
for an account, as well as creating a .netrc file that carries your login 
credentials, can be `found here`_.

.. _NASA CDDIS archive: https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/atmospheric_products.html
.. _Porayko et al. 2019: https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.4100P/abstract
.. _found here: https://cddis.nasa.gov/Data_and_Derived_Products/CDDIS_Archive_Access.html

ver 1.1.1: Fixes to propagation of Stokes headers in Stokes U (was previously copying Stokes Q header)  
ver 1.1.2: Updated to new filename convention used in CDDIS for new data.  
ver 1.1.3: JPLG data is no longer being produced (as of Aug 2023); switching default data to UQRG.  



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
