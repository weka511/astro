#!/usr/bin/env python

#   Copyright (C) 2026 Simon Crase

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

'''
Euler's method of integrating an ODE
'''

import numpy as np

class Euler:
    '''
    Euler's method of integrating an ODE
    
    A simple workhorse integrator for dx/dt = F(x)
    '''
    def __init__(self,F,n=1, atol=1.0e-8):
        '''
        Parameters:
            F       The function on the RHS
            n       Controls whether to use Exlicit Euler (n==1) or
                    iterate implicit Euler
            atol    Maximum error we can tolerate when solving implicit equations
        '''
        self.F = F
        self.n = n
        self.atol = atol
        
    def __str__(self):
        '''
        Used to display details of integrator
        '''
        if self.n == 1:
            return 'Euler'
        else:
            return f'Implicit Euler, n={self.n}'
        
    def step(self, h, x):
        '''
        Perform one step of integration
        
        Parameters:
            h       Step size
            x       Current value
        '''
        x1 = x.copy()
        for i in range(self.n):
            delta_x = h*self.F(x1)
            x1 = x + delta_x
            if np.linalg.norm(delta_x) < self.atol: return x1
        return x1