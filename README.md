# Ionospheric Faraday rotation corrections.

Scripts for calculating and applying ionospheric Faraday rotation corrections for radio polarization data, primarily time-independent corrections.

**Full documentation can be found [here](https://frion.readthedocs.io/en/latest/).**

The core underlying tool is [RMextract](https://github.com/lofar-astron/RMextract/) which calculates the time-dependent ionospheric Faraday rotation. This script uses that to generate the time-integrated correction for a given observation.

This is being written for the POSSUM pipeline, but with hopes that it can be general enough to apply to other data sets.


Cameron Van Eck, April 2021
