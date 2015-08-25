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


#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <stdexcept>
#include "mic.h"


double compute_lj_mic(double* pos, double* box, double* force, double sigma,
    double epsilon, long natom)
{
    if (natom <= 0) {
        throw std::domain_error("The number of atoms must be strictly positive.");
    }
    double energy = 0.0;

    // Double loop over all pairs of atoms.
    for (long iatom0=0; iatom0<natom; iatom0++) {
        for (long iatom1=0; iatom1<iatom0; iatom1++) {
            // Compute relative vector.
            double delta[3];
            delta[0] = pos[3*iatom0  ] - pos[3*iatom1  ];
            delta[1] = pos[3*iatom0+1] - pos[3*iatom1+1];
            delta[2] = pos[3*iatom0+2] - pos[3*iatom1+2];
            // Apply minimum image convention.
            delta[0] -= box[0]*round(delta[0]/box[0]);
            delta[1] -= box[1]*round(delta[1]/box[1]);
            delta[2] -= box[2]*round(delta[2]/box[2]);
            // Compute contribution to the total energy.
            double dsq = (delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);
            double xsq = sigma*sigma/dsq;
            double x6 = xsq*xsq*xsq;
            energy += 4*epsilon*(x6*x6 - x6);
            // Compute forces if needed.
            if (force!=NULL) {
                double gd = -24*epsilon/dsq*(2*x6*x6 - x6);
                delta[0] *= gd;
                delta[1] *= gd;
                delta[2] *= gd;
                force[3*iatom0  ] -= delta[0];
                force[3*iatom0+1] -= delta[1];
                force[3*iatom0+1] -= delta[2];
                force[3*iatom1  ] += delta[0];
                force[3*iatom1+1] += delta[1];
                force[3*iatom1+1] += delta[2];
            }
        }
    }

    return energy;
}
