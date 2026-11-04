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
    Adams Bashforth Predictor with Adams Moulton Corrector
    
    Based on Numerical Initial Value Problems in Ordinary Differential Equations,
    C. William Gear,
    ISBN 0-13-626606-1
'''

from typing import Callable
from unittest import TestCase, main
import numpy as np
import numpy.typing as npt
from numpy.testing import assert_array_almost_equal

__version__ = '1.0'
__author__ = 'Simon Crase'


class Predictor:
    '''
    This class implements the Adams-Bashforth predictor
    
    Attributes:
        dx          Function for right hand side of ODE
        beta        Vector of numerators of coefficients from Table 7.3
        divisor     Common denominator of coefficients from Table 7.3
        order       The order of the integrator
        state       Indicates whether DX has been fully populated
        DX          Matrix of x values from previous steps, used to calculate next X
        starter     Integrator used to populate DX initially
    '''
    @staticmethod
    def get_coefficients(order:int):
        '''
        Calculate the numerators for the coefficients
        
        Copied from Gear, Table 7.3
        
        Parameters:
            order      The order of the integrator
        '''
        beta = [
            [1],
            [3, -1],
            [23, -16, 5],
            [55, -59, 37, -9],
            [1901, -2774, 2616, -1274, 251],
            [4277, -7923, 9982, -7298, 2877, -475]
        ]
        return np.array(beta[order-1])
        
    def __init__(self, dx:Callable[[npt.NDArray[float]], npt.NDArray[float]], order:int, starter):
        '''
        Parameters:
            dx       Function for right hand side of ODE
            order    The order of the integrator
        '''
        self.dx = dx
        self.beta = Predictor.get_coefficients(order)
        self.divisor = self.beta.sum()
        self.order = order 
        self.state = -1
        self.DX = np.zeros((order+1,0))
        self.starter = starter
          
    def __str__(self):
        '''
        Used to display details of integrator
        '''
        return f'{type(self).__name__}: order={self.order}'

    def step(self, h:float, x:npt.NDArray[float]):
        '''
        Compute x after next step
        
        Parameters:
           h    Step size
           x    Current value of x
        '''
        if self.state == -1:
            n = len(x)
            self.DX = np.full((self.order+1, n),np.nan)
            self.state += 1
            
        if self.state < self.order:        
            self.DX[self.state, :] = self.dx(x)
            self.state += 1    
            return self.starter.step(h,x)              
        else:
            x_new = x + h* np.dot(self.beta,self.DX[:-1])/self.divisor
            self.DX[-1,:] = self.dx(x_new)
            return x_new
    
    def has_filled(self):
        '''
        Indicates whether DX has been fully populated (except, possibly, the last row)
        '''
        return self.state >= self.order
    
    def roll(self):
        '''
        Used after a step to make room for the latest fx(x)
        '''
        if self.has_filled():
            self.DX = np.roll(self.DX,-1,axis=0)

class Corrector:
    '''
    This class implements the Adams-Moulton corrector
    
    Attributes:
        beta        Vector of numerators of coefficients from Table 7.4
        divisor     Common denominator of coefficients from Table 7.4
        predictor   Used to make initial prediction of x
        n           Number of correction steps
        atol        Maximum tolerable error when we iterate
    '''
    @staticmethod
    def get_coefficients(order:int):
        '''
        Calculate the numerators for the coefficients
        
               Copied from Gear, Table 7.4
        
        Parameters:
            order      The order of the integrator
        '''        
        beta = [
            [1],
            [1, 1],
            [5, 8, -1],
            [9, 19, -5, 1],
            [251, 646, -264, 106, -19],
            [475, 1427, -798, 482, -173, 27]
        ]
        return np.array(beta[order-1])
    
    def __init__(self, predictor,n:int=1,atol=1e-6):
        '''
        Parameters:
            dx         Function for right hand side of ODE
            predictor  Used to make initial prediction of x
            n          Maximum number of correction steps
        '''
        self.beta = Corrector.get_coefficients(predictor.order)
        self.divisor = self.beta.sum()
        self.predictor = predictor
        self.n = n
        self.atol = atol
        
    def __str__(self):
        '''
        Used to display details of integrator
        '''
        return f'Adams-Bashforth-Moulton: order={self.predictor.order}'
    
    def step(self, h:float, x:npt.NDArray[float]):
        '''
        Compute x after next step
        
        Parameters:
           h    Step size
           x    Current value of x
        '''
        x_previous = np.full_like(x,np.inf)
        prediction = self.predictor.step(h,x)
        if self.predictor.has_filled():
            for i in range(self.n):
                x_new = x + h* np.dot(self.beta,self.predictor.DX[:-1])/self.divisor
                self.predictor.DX[-1,:] = self.predictor.dx(x_new)
                error = np.linalg.norm(x_new-x_previous)
                if error < self.atol: break
                x_previous = x_new
            self.predictor.roll()
            return x_new
        else:
            return prediction

class MockStarter:
    '''
    Used when testing to initialize Predictor.DX
    '''
    def __init__(self,values):
        self.values = values
        self.index = 0
        
    def step(self,h,x):
        self.index += 1
        return self.values[self.index]
    
class TestAdamsBashforth(TestCase):
    '''
    Tests to verify predictor
    '''
    def test_coefficients(self):
        '''
        Verify coefficients by comparing sums with expected values in Gear, Table 7.3
        '''
        expected_sums = [1, 2, 12, 24, 720, 1440]
        for i in range(len(expected_sums)):
            self.assertEqual(expected_sums[i],Predictor.get_coefficients(i+1).sum())
            
    def test_gear_7_1(self):
        self.helper(0.0001,100,decimal=4)
        
    def test_gear_7_1(self):
        self.helper(0.0001,100,decimal=4)        
        
    def helper(self,h,N,order=6,decimal=3):
        X = np.full((N,1),np.nan)
        Target = np.zeros_like(X)
        for i in range(N):
            Target[i,0] = np.exp(-h*i)        
        starter = MockStarter(Target)
        predictor = Predictor(lambda x: -x,order,starter)
        X[0,:] = np.array([1])
        for i in range(1,N):
            X[i] = predictor.step(h,X[i-1])
 
        assert_array_almost_equal(Target,X,decimal=decimal)
 
class TestAdamsMoulton(TestCase):  
    '''
    Tests to verify corrector
    '''    
    def test_coefficients(self):
        '''
        Verify coefficients by comparing sums with expected values in Gear, Table 7.4
        '''        
        expected_sums = [1, 2, 12, 24, 720, 1440]
        for i in range(len(expected_sums)):
            self.assertEqual(expected_sums[i],Corrector.get_coefficients(i+1).sum())    

    def test_1(self):
        self.helper(h=0.0001,N=100000,order=6,decimal=6)
    
    def helper(self,h=0.001,N=100,order=5,decimal=3):
        X = np.full((N,1),np.nan)
        Target = np.zeros_like(X)
        for i in range(N):
            Target[i,0] = np.exp(-h*i)        
        starter = MockStarter(Target)
        predictor = Predictor(lambda x: -x,order,starter)        
        corrector = Corrector(predictor, n=5)
        X[0,:] = np.array([1])
        for i in range(1,N):
            X[i] = corrector.step(h,X[i-1])
 
        assert_array_almost_equal(Target,X,decimal=3)
        
if __name__ == '__main__':
    main()
