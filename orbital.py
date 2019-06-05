# Copyright (C) 2019 Greenweaves Software Limited

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this software.  If not, see <http://www.gnu.org/licenses/>

from numpy import matmul,sign
from math import cos,sin,radians,isclose,degrees,pi

def compose(omega=0,I=0,Omega=0):
    cos_omega = cos(radians(omega))
    sin_omega = sin(radians(omega))
   
    P1        = [[cos_omega, -sin_omega, 0],\
                 [sin_omega,  cos_omega, 0],\
                 [0,          0,         1]]
    
    cos_I     = cos(radians(I))
    sin_I     = sin(radians(I))
   
    P2        = [[1,          0,      0],     \
                 [0,          cos_I, -sin_I], \
                 [0,          sin_I,  cos_I]]
    
    cos_Omega = cos(radians(Omega))
    sin_Omega = sin(radians(Omega))
    
    P3        =  [[cos_Omega, -sin_Omega, 0],\
                  [sin_Omega,  cos_Omega, 0],\
                  [0,          0,         1]]
   
    return matmul(P3, matmul(P2, P1))


def kepler(eccentricy=0,mean_anomaly=0,tolerance=0.1e-7,N=10000,k=0.85):
    M = mean_anomaly % (2 * pi)
    
    E = M + sign(sin(M)*k*eccentricy)
    for i in range(N):
        correction = (E - eccentricy*sin(E) - M)/(1-eccentricy*cos(E))
        if abs(correction)<tolerance: return E
        E -= correction
    
if __name__=='__main__':
    import unittest
    
    class TestJupiter(unittest.TestCase):
        def test_mult(self):
            print(compose(omega=14.7392,I=1.30537,Omega=100.535))
            
    class TestKepler(unittest.TestCase):
        def test_kepler(self):
            eccentric_anomaly = kepler(eccentricy=0.205635,mean_anomaly=1.2)
            M = eccentric_anomaly - 0.205635 * sin(eccentric_anomaly)
            self.assertAlmostEqual(1.2,M)
        def test_kepler_jupiter(self):
            self.assertAlmostEqual(radians(189.059),
                                   kepler(mean_anomaly=radians(189.495),
                                          eccentricy=0.0484007),
                                   places=4)
            
    unittest.main()