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


cdef extern from "mic.h":
    double compute_lj_mic(double* pos, double* box, double* force, double sigma,
        double epsilon, long natom) except +
