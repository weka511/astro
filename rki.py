#!/usr/bin/env python

# Copyright (C) 2015-2026 Greenweaves Software Pty Ltd

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

'''
Implicit Runge Ketta (symplectic) integrators
'''

from abc import ABC
from argparse import ArgumentParser
from pathlib import Path
from time import time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np

__version__ = '1.0'
__author__ = 'Simon Crase'

class Driver:
    '''
    Use an integrator to solve an ODE to within a specified error 
    
    Attributes:
        integrator    Used for integrating ODE
        h_minimum     Step size cannot be increased beyond this value
        h             Initial step size (variable)
        h_maximum     Step size cannot be increased beyond this value
        atol          Maximum tolerable error
        mult          Used to set a lower bound for error (fraction of epsilon)
    '''
    def __init__(self, integrator, 
                 h_minimum=1.0e-9, h=0.5, h_maximum=1.0, atol=1.0e-9, mult=0.01, 
                 get_distance=lambda k, k_new: max([abs(a - b) for (a, b) in zip(k, k_new)])):
        '''
        Initialize Driver
        
        Parameters:
            integrator    Used for integrating ODE
            h_minimum     Step size cannot be increased beyond this value
            h             Initial step size (variable)
            h_maximum     Step size cannot be increased beyond this value
            atol          Maximum tolerable error
            mult          Used to set a lower bound for error (fraction of epsilon)     
        '''
        self.integrator = integrator
        self.h_minimum = h_minimum
        self.h_maximum = h_maximum
        self.epsilon = atol
        self.h = h
        self.min_epsilon = mult * atol
 
    def step(self, y):
        '''
        Solve for one step, varying step size to keep truncation error within bounds
        
        Parameters:
           y        Current value of dependent variable
           
        Returns:
            Value of dependent variable at the end of the step
        '''
        try:
            y1 = self.integrator.step(self.h, y)  # one step estimate of next value
            
            y11 = self.integrator.step(0.5 * self.h,
                                       self.integrator.step(0.5 * self.h, y)) # Two half steps
            
            error = np.linalg.norm(y1- y11) # Estimate error by comparing the two estimates for next y

            if error > self.epsilon:  # too large - reduce step size
                self.h *= (self.epsilon/error)**(1.0/self.integrator.order)
                return self.step(y)
            
            self.adjust_stepsize(error)
            if self.h < self.h_minimum:
                self.h = self.h_minimum                
                    
            return y11
        except ImplicitRungeKutta.Failed:
            self.h *= 0.5
            return self.step(y)

    def adjust_stepsize(self,error):
        '''
        Verify that our stepsize isn't too small or too large
        
        Parameters:
            error
        '''
        if error > self.min_epsilon: return
        
        self.h *= ((self.min_epsilon/error)**(1.0/self.integrator.order) if error > 0 else 2.0)

        if self.h > self.h_maximum:
            self.h = self.h_maximum    

class ImplicitRungeKutta(ABC):
    '''
    Parent class for Implicit Runge-Kutta integrators
    
    see https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Gauss.E2.80.93Legendre_methods
    
    The Butcher tableau is stored in three members: a, b, and c.
    '''
    class Failed(Exception):
        ''' 
        Exception thrown when we can't solve implicit equations
        '''
        def __init__(self, value):
                self.value = value

        def __str__(self):
            return repr(self.value)

    def __init__(self, dy, max_iterations, atol, order, a,b,c):
        '''
        Initialize
        
        Parameters:
           dy              Function f for ODE = dy/dx=f(y)
           max_iterations  Maximum number of iterations for solving implicit equations
           atol            Maximum error we can tolerate when solving implicit equations
           order           Order of solver
           a               The matrix from the Butcher tableau https://mathworld.wolfram.com/ButcherTableau.html
           b               The b vector from the Butcher tableau
           c               The c vector from the Butcher tableau
        '''
        self.dy = dy
        self.max_iterations = max_iterations
        self.atol = atol
        self.order = order
        self.a = a
        self.b = b
        self.c = c
        self.s = len(b)
        
    def __str__(self):
        '''
        Used to display details of integrator
        '''
        return f'{type(self).__name__}: order={self.order},max_iterations={self.max_iterations},atol={self.atol}'

    def step(self, h, y):
        '''
        Compute y after next step
        
        Parameters:
           h    Step size
           y    Current value of y
        '''
        diff = 0.0
        k = np.zeros((len(self.b),len(y)))
        for _ in range(self.max_iterations):
            k_new = self.iterate(h, y, k)
            diff = np.linalg.norm(k - k_new)
            if diff < self.atol:
                return y + h*np.inner(self.b,k_new.T)
            else:
                k = k_new
                
        raise ImplicitRungeKutta.Failed(
            f'Failed to converge within {self.atol} after {self.max_iterations} iterations: error={diff}')

    def iterate(self, h, y, k):
        '''Iterate One step in solution of iterative equations for k
        
        Parameters:
            h       Step size
            y       Current value of solution
            k       The k term in Eunge Kutta
        '''
        result = np.zeros((self.s,len(y)))
        for i in range(self.s):
            result[i,:] = self.dy(y + h*np.inner(self.a[i,:],k.T))
        return result

class ImplicitRungeKutta2(ImplicitRungeKutta):
    '''
    4th order Gauss-Legendre
    '''
    def __init__(self, dy, max_iterations, atol):
        '''
        Parameters:
           dy              Function f for ODE = dy/dx=f(y)
           max_iterations  Maximum number of iterations for solving implicit equations
           atol            Maximum error we can tolerate when solving implicit equations
        '''
        super().__init__(dy, max_iterations, atol, 4,
                np.array([
                    [0.25, 0.25 - np.sqrt(3.0)/6.0],
                    [0.25 + np.sqrt(3.0)/6.0, 0.25],
                ]),
                np.array([
                    0.5, 0.5
                ]),
                np.array([
                    0.5 - np.sqrt(3.0)/6,
                    0.5 + np.sqrt(3.0)/6
                ]))
   
class ImplicitRungeKutta4(ImplicitRungeKutta):
    '''
    6th order Gauss-Legendre
    '''

    def __init__(self, dy, max_iterations=200, atol=1.0e-18):
        '''
        Parameters:
           dy              Function f for ODE = dy/dx=f(y)
           max_iterations  Maximum number of iterations for solving implicit equations
           atol            Maximum error we can tolerate when solving implicit equations
        '''
        super().__init__(dy, max_iterations, atol, 6,
                         np.array([
                                [5.0/36.0, 2.0/9.0 -np.sqrt(15.0)/15.0, 5.0/36.0 - np.sqrt(15.0)/30.0],
                                [5.0/36.0 + np.sqrt(15.0)/24.0, 2.0/9.0, 5.0/36.0 - np.sqrt(15.0)/24.0],
                                [5.0/36.0 + np.sqrt(15.0)/30.0, 2.0/9.0 + np.sqrt(15.0)/15.0, 5.0/36.0]
                            ]),
                         np.array([
                                5.0/18.0,
                                4.0/9.0,
                                5.0/18.0
                            ]),
                         np.array([
                                0.5 - np.sqrt(15.0)/10.0,
                                0.5,
                                0.5 + np.sqrt(15.0)/10.0
                            ]))

def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
    return parser.parse_args()

def run_once(ax,rk,nn = 100):
    driver = Driver(rk, 0.000000001, 0.5, 1.0, 0.000000001)
    try:
        y = np.array([1, 0])
        xs = []
        ys = []
        for i in range(nn):
            y = driver.step(y)
            xs.append(y[0])
            ys.append(y[1])           
        ax.plot(xs, ys)
        ax.set_title(f'{type(rk).__name__}: nn={nn}')
    except ImplicitRungeKutta.Failed as e:
        print(f'Exception {e}')
    
def main():
    args = parse_args()
    rcParams['text.usetex'] = True
    start  = time()    
    fig = figure(figsize=(12, 12))
    fig.suptitle(Path(__file__).stem)
    run_once(ax=fig.add_subplot(1,2,1,adjustable='box',aspect=1.0),
             rk=ImplicitRungeKutta2(lambda y: [y[1], -y[0]], 10, 0.000000001),
             nn=150)
    run_once(ax=fig.add_subplot(1,2,2,adjustable='box',aspect=1.0),
             rk=ImplicitRungeKutta4(lambda y: [y[1], -y[0]], 10, 0.000000001),
             nn=50)    
    
    fig.tight_layout(h_pad=2)
    fig.savefig(Path(args.figs)/Path(__file__).stem)    
    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
    
    if args.show:
        show()

if __name__ == '__main__':
    main()
