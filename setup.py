#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
import sys
from shutil import rmtree

from setuptools import find_packages, setup, Command

NAME = 'FRion'
DESCRIPTION = 'Ionospheric Faraday rotation prediction and correction for radio astronomy polarization cubes.'
URL = 'https://github.com/Cameron-Van-Eck/FRion'
REQUIRES_PYTHON = '>=3.5.0'
VERSION = '1.1.4'
DOWNLOAD_URL = 'https://github.com/Cameron-Van-Eck/FRion/archive/refs/heads/main.zip'

REQUIRED = [
    'numpy', 'astropy', 'pyephem', 'requests', 'scipy'
]

extras_require={}

here = os.path.abspath(os.path.dirname(__file__))

try:
    with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    python_requires=REQUIRES_PYTHON,
    url=URL,
    download_url=DOWNLOAD_URL,
    packages=['FRion'],
    entry_points={
        'console_scripts': ['frion_predict=FRion.predict:predict',
                            'frion_timeseries=FRion.predict:timeseries',
                            'frion_correct=FRion.correct:command_line'],
    },
    install_requires=REQUIRED,
    include_package_data=True,
    license='MIT',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    maintainer='Cameron Van Eck',
    maintainer_email='cameron.vaneck@anu.edu.au',
)
