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

import numpy as np
from math import cos,sin,radians,isclose

def compose(omega,I,Omega):
    cos_omega = cos(radians(omega))
    sin_omega = sin(radians(omega))
    assert isclose(cos_omega*cos_omega + sin_omega*sin_omega,1,abs_tol=1e-10)
    P1        = [[cos_omega, -sin_omega, 0],\
                 [sin_omega,  cos_omega, 0],\
                 [0,          0,         1]]
    assert isclose(np.linalg.det(P1),1,abs_tol=1e-10)
    cos_I     = cos(radians(I))
    sin_I     = sin(radians(I))
    assert isclose(cos_I*cos_I + sin_I*sin_I,1,abs_tol=1e-10)
    P2        = [[1,          0,      0],     \
                 [0,          cos_I, -sin_I], \
                 [0,          sin_I,  cos_I]]
    assert isclose(np.linalg.det(P2),1,abs_tol=1e-10)
    cos_Omega = cos(radians(Omega))
    sin_Omega = sin(radians(Omega))
    assert isclose(cos_Omega*cos_Omega + sin_Omega*sin_Omega,1,abs_tol=1e-10)
    P3        =  [[cos_Omega, -sin_Omega, 0],\
                  [sin_Omega,  cos_Omega, 0],\
                  [0,          0,         1]]
    assert isclose(np.linalg.det(P3),1,abs_tol=1e-10)
    Pj = np.matmul(P3, np.matmul(P2, P1))
    assert isclose(np.linalg.det(Pj),1,abs_tol=1e-10)
    return Pj
    
if __name__=='__main__':
    import unittest
    
    class Test1(unittest.TestCase):
        def test1(self):
            a = [[1, 1],\
                 [0, 1]]
            b = [[4, 1], \
                 [2, 2]]
            print (np.matmul(a, b))            
            self.assertEqual(1,1)
            print(compose(14.7392,1.30537,100.535))

    unittest.main()