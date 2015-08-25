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
cimport numpy as np
np.import_array()


cimport mic


__all__ = ['compute_lj_mic']


def _check_c_safe(array, shape, name):
    if not array.flags['C_CONTIGUOUS']:
        raise TypeError('The argument %s must be a C-contiguous array.' % name)
    for ishape in xrange(len(shape)):
        if shape[ishape] is not None and shape[ishape] != array.shape[ishape]:
            raise TypeError('Array %s has wrong dimension in axis %i. '
                            'Expected %i but got %i' % (
                            name, ishape, shape[ishape], array.shape[ishape]))

def compute_lj_mic(np.ndarray[double, ndim=2] pos not None,
                   np.ndarray[double, ndim=1] box not None,
                   np.ndarray[double, ndim=2] force,
                   double sigma, double epsilon):
    # Good old type declaration
    cdef np.ndarray[double, ndim=2] tmp

    _check_c_safe(pos, (None, 3), 'pos')
    natom = pos.shape[0]
    _check_c_safe(box, (3, ), 'box')
    if force is None:
        return mic.compute_lj_mic(&pos[0,0], &box[0], NULL, sigma, epsilon, natom)
    else:
        _check_c_safe(force, (natom, 3), 'force')
        # Make sure force is assigned to a variable with guaranteed type.
        tmp = force
        return mic.compute_lj_mic(&pos[0,0], &box[0], &tmp[0,0], sigma, epsilon, natom)
