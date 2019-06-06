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

from numpy import matmul,sign,arctan2
from math import cos,sin,radians,isclose,degrees,pi,sqrt

def compose(omega=0,I=0,Omega=0):
    cos_omega = cos(omega)
    sin_omega = sin(omega)
   
    P1        = [[cos_omega, -sin_omega, 0],\
                 [sin_omega,  cos_omega, 0],\
                 [0,          0,         1]]
    
    cos_I     = cos(I)
    sin_I     = sin(I)
   
    P2        = [[1,          0,      0],     \
                 [0,          cos_I, -sin_I], \
                 [0,          sin_I,  cos_I]]
    
    cos_Omega = cos(Omega)
    sin_Omega = sin(Omega)
    
    P3        =  [[cos_Omega, -sin_Omega, 0],\
                  [sin_Omega,  cos_Omega, 0],\
                  [0,          0,         1]]
   
    return matmul(P3, matmul(P2, P1))

def get_XYZ(omega=0,I=0,Omega=0,f=0,r=1):
    cos_Omega   = cos(Omega)
    sin_Omega   = sin(Omega)
    cos_I       = cos(I)
    sin_I       = sin(I)
    cos_omega_f = cos(omega + f)
    sin_omega_f = sin(omega + f)
    return (r * (cos_Omega*cos_omega_f - sin_Omega*sin_omega_f*cos_I),
            r * (sin_Omega*cos_omega_f + cos_Omega*sin_omega_f*cos_I),
            r * sin_omega_f*sin_I)
    
def kepler(eccentricy=0,mean_anomaly=0,tolerance=0.1e-12,N=10000,k=0.85):
    M = mean_anomaly % (2 * pi)
    E = M + sign(sin(M)*k*eccentricy)
    
    for i in range(N):
        correction = (E - eccentricy*sin(E) - M)/(1 - eccentricy * cos(E))
        if abs(correction)<tolerance: return E
        E -= correction
    
if __name__=='__main__':
    import unittest
    
    class TestJupiter(unittest.TestCase):
        
        # test_mult
        #
        # This test uses Murray & Dermott (2.124)
        def test_compose(self):
            rotation = compose(omega=radians(14.7392 -100.535 ),I=radians(1.30537),Omega=radians(100.535))
            self.assertAlmostEqual(0.966839,   rotation[0][0], places=5)
            self.assertAlmostEqual(-0.254401,  rotation[0][1], places=6)
            self.assertAlmostEqual(0.0223971,  rotation[0][2], places=6)
            self.assertAlmostEqual(0.254373,   rotation[1][0], places=5)
            self.assertAlmostEqual(0.967097,   rotation[1][1], places=6)
            self.assertAlmostEqual(0.00416519, rotation[1][2], places=6)
            self.assertAlmostEqual(-0.0227198, rotation[2][0], places=6)
            self.assertAlmostEqual(0.00167014, rotation[2][1], places=7)
            self.assertAlmostEqual(0.99974,    rotation[2][2], places=6)
        
        # test_mult
        #
        # This test uses Murray & Dermott (2.124)            
        def test_get_XYZ(self):
            x     = -5.39027
            y     = -0.818277
            r     = sqrt(x*x + y * y)
            f     = arctan2(y,x)
            X,Y,Z = get_XYZ(omega=radians(14.7392 -100.535 ),I=radians(1.30537),Omega=radians(100.535),r=r,f=f)
            self.assertAlmostEqual(-5.00336, X, places=5)
            self.assertAlmostEqual(-2.16249, Y, places=5)
            self.assertAlmostEqual(0.121099, Z, places=6)
            
    class TestKepler(unittest.TestCase):
        def test_kepler_inverse(self):
            eccentric_anomaly = kepler(eccentricy=0.205635,mean_anomaly=1.2)
            M = eccentric_anomaly - 0.205635 * sin(eccentric_anomaly)
            self.assertAlmostEqual(1.2,M)
        def test_kepler_jupiter(self):
            self.assertAlmostEqual(radians(189.059),
                                   kepler(mean_anomaly=radians(189.495),
                                          eccentricy=0.0484007),
                                   places=4)
            
    unittest.main()