#!/usr/bin/env python3

# This file is part of Planets.
#
# Planets is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Planets is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Planets.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from gravity import gravity

def test_preserve_direction():
    m = np.array([1., 2.]) # M_Earth
    r = np.array([0.001, 0.001]) # AU

    for x1,x2 in [(-1.,1.), (-2.,1.), (-1.,2.)]:
        for v1,v2 in [(0.,0.), (-2.,1.), (-1.,2.)]:
            x = np.asfortranarray([[x1, 0.], [x2, 0.]]) # AU
            v = np.asfortranarray([[v1, 0.], [v2, 0.]]) # AU/day

            gravity.integrate(1e5, 1.e-5, x, v, m, r, len(m))

            assert(all(v[:,1] == 0))
            assert(all(x[:,1] == 0))

            x = np.asfortranarray([[0., x1], [0., x2]]) # AU
            v = np.asfortranarray([[0., v1], [0., v2]]) # AU/day

            gravity.integrate(1e5, 1.e-5, x, v, m, r, len(m))

            assert(all(v[:,0] == 0))
            assert(all(x[:,0] == 0))

def test_preserve_symmetry():
    m = np.array([1., 1.]) # M_Earth
    r = np.array([0.001, 0.001]) # AU
    x = np.asfortranarray([[-1., 1.], [1., -1.]]) # AU
    v = np.asfortranarray([[0., 0.], [0., 0.]]) # AU/day

    gravity.integrate(1000, 1.e-5, x, v, m, r, len(m))

    assert(abs(sum(x[:,0]) / x[0,0]) < 1.e-15)
    assert(abs(sum(x[:,1]) / x[0,1]) < 1.e-15)
    assert(abs(sum(v[:,0]) / v[0,0]) < 1.e-15)
    assert(abs(sum(v[:,1]) / v[0,1]) < 1.e-15)

def test_orbit():
    m = np.array([332953., 1.]) # Sun and Earth
    r = np.array([0.001, 0.001]) # AU
    x = np.asfortranarray([[0., 0.], [1., 0.]]) # AU
    v = np.asfortranarray([[0., 0.], [0., 0.0172]]) # AU/day

    L0 = m[1] * v[1,1]
    E0 = 0.5 * v[1,1]**2 - gravity.g*(m[0]+m[1])/x[1,0]

    gravity.integrate(9129, 1.e-2, x, v, m, r, len(m))

    L1 = m[1] * (x[1,0]*v[1,1] - x[1,1]*v[1,0]) + m[0] * (x[0,0]*v[0,1] - x[0,1]*v[0,0])
    E1 = 0.5 * (v[1,0]**2 + v[1,1]**2) - gravity.g*(m[0]+m[1])/np.sqrt(x[1,0]**2 + x[1,1]**2)
    assert(abs(L0 - L1) / L0 < 1.e-8)
    assert(abs(E0 - E1) / E0 < 1.e-15)

    assert(abs(x[1,0]) < 1.e-5)
    assert(abs(x[1,1] - 1.) < 1.e-3)
    assert(abs(v[1,0] + 0.0172) < 1.e-5)
    assert(abs(v[1,1]) < 1.e-5)

def test_collision():
    m = np.array([1., 1.]) # M_Earth
    r = np.array([1., 1.]) # AU
    x = np.asfortranarray([[-1., 0.], [1., 0.]]) # AU
    v = np.asfortranarray([[0., 0.], [0., 0.]]) # AU/day

    gravity.integrate(1, 1., x, v, m, r, len(m))

    assert(m[0] == 2.)
    assert(m[1] == 0.)
    assert(r[0] > 1.)
    assert(all(x[0,:] == 0.))
    assert(all(v[0,:] == 0.))

test_preserve_direction()
test_preserve_symmetry()
test_orbit()
test_collision()
