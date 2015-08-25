// Python-MIC is a Python MD demo with energy (and forces) computed in C++.
// Copyright (C) 2015 Toon Verstraelen
//
// This file is part of Python-MIC.
//
// Python-MIC is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// Python-MIC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>
//
//--


#ifndef PYTHON_MIC_MIC_H
#define PYTHON_MIC_MIC_H

/**
    @brief
        Computes Lennard-Jones energy (and forces) between atoms in a periodic
        box, using the minimum image convention.

    @param pos
        A pointer to an array with atomic positions. The shape is assumed to be
        (natom, 3).

    @param box
        A pointer to the lengths of the orthorombic box. The shape is assumed to
        be (3,).

    @param force
        A pointer to an output array with atomic forces. Same shape as pos. If
        NULL, the forces are not computed.

    @param sigma, epsilon
        Duhh.

    @paran natom
        The number of atoms. This must be strictly positive.

    This routine does not assume that the atoms are all located within the
    orthorombic box. Relative vectors between a pair of atoms are always reduced
    using the minimum image convention to the shortes possible relative vector.

 */
double compute_lj_mic(double* pos, double* box, double* force, double sigma,
    double epsilon, long natom);

#endif
