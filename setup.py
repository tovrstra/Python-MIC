#!/usr/bin/env python
# Python-MIC is a Python MD demo with energy (and forces) computed in C++.
# Copyright (C) 2015 Toon Verstraelen
#
# This file is part of Python-MIC.
#
# Python-MIC is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# Python-MIC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http:#www.gnu.org/licenses/>
#
#--

import numpy as np
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


setup(
    name='Python-MIC',
    version='0.0.0',
    description='Python-MIC is a Python MD demo with energy (and forces) computed in C++.',
    author='Toon Verstraelen',
    author_email='Toon.Verstraelen@UGent.be',
    py_modules=['verlet.py'],
    cmdclass = {'build_ext': build_ext},
    ext_modules=[
        Extension("wrapper",
            sources=['mic.cpp', 'wrapper.pyx'],
            depends=['mic.h', 'mic.pxd'],
            include_dirs=[np.get_include()],
            language="c++")],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2',
        'Programming Language :: Cython',
        'Programming Language :: C++',
        'Topic :: Science/Engineering :: Molecular Science'
    ],
)
