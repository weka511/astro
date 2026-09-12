#!/usr/bin/env python

# Copyright (C) 2015-2017,2026 Greenweaves Software Pty Ltd

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
        integrator  Used for integrating ODE
        h_minimum   Step size cannot be increased beyond this value
        h           Initial step size (variable)
        h_maximum   Step size cannot be increased beyond this value
        atol        Maximum tolerable error
        mult        Used to set a lower bound for error (fraction of epsilon)
        distance    Distance function used to compute error when we solve equations  
    '''
    def __init__(self, integrator, 
                 h_minimum=1.0e-9, h=0.5, h_maximum=1.0, atol=1.0e-9, mult=0.01, 
                 distance=lambda k, k_new: max([abs(a - b) for (a, b) in zip(k, k_new)])):
        '''
        Initialize Driver
        
        Parameters:
            integrator  Used for integrating ODE
            h_minimum   Step size cannot be increased beyond this value
            h           Initial step size (variable)
            h_maximum   Step size cannot be increased beyond this value
            atol        Maximum tolerable error
            mult        Used to set a lower bound for error (fraction of epsilon)
            distance    Distance function used to compute error when we solve equations       
        '''
        self.integrator = integrator
        self.h_minimum = h_minimum
        self.h_maximum = h_maximum
        self.epsilon = atol
        self.h = h
        self.min_epsilon = mult * atol
        self.distance = distance

    def step(self, y):
        '''
        Solve for one step, varying step size to keep weeor within bounds
        
        Parameters:
           y        Current value of dependent variable
        '''
        try:
            y1 = self.integrator.step(self.h, y)  # one step estimate of next value
            
            y11 = self.integrator.step(0.5 * self.h,
                                       self.integrator.step(0.5 * self.h, y)) # Two half steps
            
            error = self.integrator.distance(y1, y11) # Estimate error by comparing the two estimates for next y

            if error > self.epsilon:  # too large - reduce step size
                self.h *= (self.epsilon/error)**(1.0/self.integrator.order)
                return self.step(y)

            if error < self.min_epsilon:  # Is our stepsize too small?
                if error > 0:
                    self.h *= (self.min_epsilon/error)**(1.0/self.integrator.order)
                else:
                    self.h *= 2.0

                if self.h > self.h_maximum:
                    self.h = self.h_maximum
            return y11
        except ImplicitRungeKutta.Failed:
            self.h *= 0.5
            return self.step(y)


class ImplicitRungeKutta(ABC):
    '''
    Parent class for Implicit Ringe-Kutta integrators
    
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

    def __init__(self, dy, max_iterations, max_iteration_error, order, a,b,c,
                 distance=lambda k0, k1: max([abs(a - b) for (a, b) in zip(k0, k1)])
                 ):
        '''
        Initialize
        
        Parameters:
           dy                  Function f for ODE = dy/dx=f(y)
           max_iterations      Maximum number of iterations for solving implicit equations
           max_iteration_error Maximum error we can tolerate when solving implicit equations
           order               Order of solver
           distance            Distance function used to compute error when we solve equations
           a
           b
           c
        '''
        self.dy = dy
        self.max_iterations = max_iterations
        self.max_iteration_error = max_iteration_error
        self.order = order
        self.distance = distance
        self.a = a
        self.b = b
        self.c = c
        self.s = len(b)

    def step(self, h, y):
        '''
        Compute y after next step
        
        Parameters:
           h    Step size
           y    Current value of y
        '''
        k = np.zeros((len(self.b),len(y)))
        for i in range(self.max_iterations):
            k_new = self.iterate(h, y, k)
            if min([self.distance(k0, k1) for (k0, k1) in zip(k, k_new)]) < self.max_iteration_error:
                yy = [y0 for y0 in y]
                for l in range(len(yy)):
                    yy[l] += h * np.dot(self.b,k_new[:,l])
                return yy
            else:
                k = k_new
        self.fail()

    def iterate(self, h, y, k):
        '''Iterate One step in solution of iterative equations for k
        
        Parameters:
            h
            y
            k
            
        assume that each row of the matrix is a single vector a = [[row1,...]] 
        '''
        result = np.zeros((self.s,len(y)))
        for i in range(self.s):
            yy = y.copy()
            for l in range(len(k[i])):
                yy[l] += h * np.dot(self.a[i,:], k[:,l])
  
            result[i] = self.dy(yy)
        return result

    def fail(self):
        '''
        Used  to throw Exception if we cannot solve implicit equations
        '''
        raise ImplicitRungeKutta.Failed(
            f'Failed to converge within {self.max_iteration_error} after {self.max_iterations} iterations')


class ImplicitRungeKutta2(ImplicitRungeKutta):
    '''
    4th order Gauss-Legendre
    '''
    def __init__(self, dy, max_iterations, max_iteration_error):
        super().__init__(dy, max_iterations, max_iteration_error, 4,
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
        super().__init__(dy, max_iterations, atol, 6,
                         np.array([
                                [5.0/36.0, 2.0/9.0 -np.sqrt(15.0)/15.0, 5.0/36.0 -np.sqrt(15.0)/30.0],
                                [5.0/36.0 +np.sqrt(15.0)/24.0, 2.0/9.0, 5.0/36.0 -np.sqrt(15.0)/24.0],
                                [5.0/36.0 +np.sqrt(15.0)/30.0, 2.0/9.0 +np.sqrt(15.0)/15.0, 5.0/36.0]
                            ]),
                         np.array([
                                5.0/18.0,
                                4.0/9.0,
                                5.0/18.0
                            ]),
                         np.array([
                                0.5 -np.sqrt(15.0)/10.0,
                                0.5,
                                0.5 +np.sqrt(15.0)/10.0
                            ]))


def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
    return parser.parse_args()

def run_once(ax,rk,):
    driver = Driver(rk, 0.000000001, 0.5, 1.0, 0.000000001)
    try:
        nn = 1000
        y = np.array([1, 0])
        xs = []
        ys = []
        for i in range(nn):
            y = driver.step(y)
            xs.append(y[0])
            ys.append(y[1])           
        ax.plot(xs, ys)
        ax.set_title(type(rk).__name__)
    except ImplicitRungeKutta.Failed as e:
        print(f'Exception {e}')
    
def main():
    args = parse_args()
    rcParams['text.usetex'] = True
    start  = time()    
    fig = figure(figsize=(12, 12))
    fig.suptitle(Path(__file__).stem)
    run_once(ax=fig.add_subplot(1,2,1,adjustable='box',aspect=1.0),
             rk=ImplicitRungeKutta2(lambda y: [y[1], -y[0]], 10, 0.000000001))
    run_once(ax=fig.add_subplot(1,2,2,adjustable='box',aspect=1.0),
             rk=ImplicitRungeKutta4(lambda y: [y[1], -y[0]], 10, 0.000000001))    
    
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
