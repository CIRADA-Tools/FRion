# Ionospheric Faraday rotation corrections.

Scripts for calculating and applying ionospheric Faraday rotation corrections for radio polarization data, primarily time-independent corrections.

**Full documentation can be found [here](https://frion.readthedocs.io/en/latest/).**
Note that there are significant changes in **Version 1.1**, which are documented there.

The core underlying tool is [RMextract](https://github.com/lofar-astron/RMextract/) which calculates the time-dependent ionospheric Faraday rotation. This script uses that to generate the time-integrated correction for a given observation.

This is being written for the POSSUM pipeline, but with hopes that it can be general enough to apply to other data sets.

This packages original author, and maintainer as of Mar 2022, is Cameron Van Eck (cameron.van.eck (at) utoronto.ca).
Please submit bug reports and feature requests fo the GitHub issues page, and feel free to email me with questions or comments.

More information on the Canadian Initiative for Radio Astronomy Data Analysis (CIRADA) can be found at cirada.ca.

