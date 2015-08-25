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
from wrapper import compute_lj_mic


__all__ = ['initialize_configuration', 'verlet_driver']


def dist_mic(pos0, pos1, box):
    '''Compute distance between two points in minimum image convention

    Parameters
    ----------
    pos0, pos1 : numpy array with shape (3,), float
                 Two positions (of atoms).
    box : numpy array with shape (3,), float
          The box edge lengths.

    Returns
    -------
    d : float
        The distance between the two points in the minumum image convention.
    '''
    delta = pos0 - pos1
    delta -= box*np.round(delta/box)
    return np.linalg.norm(delta)


class InitializationError(Exception):
    pass


def initialize_configuration(box, temp, sigma, natom, mass, maxtry=1000):
    '''Create random initial configuration that is not completely stupid

    The initial positions are guaranteed to be separated by distance sigma.
    The initial velocities are sampled from the Maxwell-Botlzmann distribution.

    Parameters
    ----------
    box : numpy array with shape (3,), float
          The box edge lengths.
    temp : float
           The temperature in reduced units. (The Boltzmann constant is 1 in
           reduced units.)
    sigma : float
            The minimum distance between two atoms.
    natom : int
            The number of atoms.
    maxtry : int, optional
             The maximum number of attempts to add a new atom that is not too
             close to the previous ones, before giving up. (default=100)

    Returns
    -------
    pos : numpy array with shape (natom, 3), type float
          The random initial positions.
    vel : numpy array with shape (natom, 3), type float
          The random initial velocities.

    Raises
    ------
    InitializationError
        If an atom can not be added after maxtry attempts without being too
        close to any previously added atom.
    '''
    vel = np.random.normal(0, np.sqrt(temp/mass), (natom, 3))
    pos = np.zeros((natom, 3))

    iatom = 0
    ntry = 0
    while iatom < natom:
        pos[iatom] = np.random.uniform(0, 1)*box
        too_close = False
        for iother in xrange(iatom):
            d = dist_mic(pos[iatom], pos[iother], box)
            if d < sigma:
                too_close = True
                break
        if not too_close:
            iatom += 1
            ntry = 0
        else:
            ntry += 1
        if ntry > maxtry:
            raise InitializationError('Could not add new atom in less than %i attempts.' % maxtry)

    return pos, vel


def verlet_driver(energy_force_fn, pos, vel, mass, h, nstep):
    '''Run a simple NVE MD simulation

    This function prints a line for every MD step with the time step, kinetic
    energy, potential energy and the total energy (conserved quantity).

    Parameters
    ----------
    energy_force_fn : function
                      A function that takes positions as argument and returns an
                      energy. It must also have an optional output argument
                      that, when given, is used to compute the forces.
    pos : numpy array with shape (natom, 3), type float
          The initial positions. This array will be modified in-place.
    vel : numpy array with shape (natom, 3), type float
          The initial velocities. This array will be modified in-place.
    mass : float
           The mass of each particle.
    h : float
        The time step.
    nstep : int
            The number of time steps.

    Local variables
    ---------------
    force : numpy array with shape (natom, 3), type float
            The array with forces.
    acc : numpy array with shape (natom, 3), type float
          The array with particle accelarations.
    delta : numpy array with shape (natom, 3), type float
            The displacement in the verlet algorithm
    epot : float
           The potential energy.
    ekin : float
           The kinetic energy.
    time : float
           The current time.
    istep : int
            The step number.
    '''
    # Allocation and initialization of some arrays
    force = np.zeros(pos.shape)
    epot = energy_force_fn(pos, force)
    ekin = 0.5*mass*(vel**2).sum()
    acc = force/mass
    time = 0

    print ' Step     Time       E_kin       E_pot      E_cons'
    print '--------------------------------------------------'
    print ' Init  %7.1f  %10.5f  %10.5f  %10.5f' % (time, ekin, epot, ekin+epot)

    for istep in xrange(nstep):
        delta = h*vel + (0.5*h**2)*acc
        pos += delta
        epot = energy_force_fn(pos, force)
        acc_new = force/mass
        vel += 0.5*(acc+acc_new)*h
        ekin = 0.5*mass*(vel**2).sum()
        acc = acc_new
        time += h

        print '%5i  %7.1f  %10.5f  %10.5f  %10.5f' % (istep, time, ekin, epot, ekin+epot)
    print '--------------------------------------------------'


class LJEnergyForceFn(object):
    def __init__(self, box, epsilon, sigma):
        self.box = box
        self.epsilon = epsilon
        self.sigma = sigma

    def __call__(self, pos, force):
        force[:] = 0.0
        return compute_lj_mic(pos, self.box, force, self.sigma, self.epsilon)


if __name__ == '__main__':
    box = np.array([100.0, 100.0, 100.0])
    natom = 100
    epsilon = 1.0
    sigma = 1.0
    temp = 1.0
    mass = 1.0
    h = 0.01
    nstep = 1000
    energy_force_fn = LJEnergyForceFn(box, epsilon, sigma)
    pos, vel = initialize_configuration(box, temp, sigma, natom, mass)
    verlet_driver(energy_force_fn, pos, vel, mass, h, nstep)
