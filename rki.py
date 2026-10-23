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
of 4th order and 5th order
References:
  [1] Numerical initial value problems in ordinary differential equationsby C. William Gear 
  [2] https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Gauss.E2.80.93Legendre_methods
  [3] https://mathworld.wolfram.com/ButcherTableau.html
'''

from abc import ABC,abstractmethod
from argparse import ArgumentParser
from pathlib import Path
from time import time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np

__version__ = '1.1'
__author__ = 'Simon Crase'

class K_Solver(ABC):
    '''
    This class is used to solve the implicit equations for k
    '''
    def __init__(self,atol=1e-8,max_iterations=200):
        '''
        Parameters:
            max_iterations  Maximum number of iterations for solving implicit equations
            atol            Maximum error we can tolerate when solving implicit equations
        '''
        self.rki = None
        self.atol = atol
        self.max_iterations = max_iterations
   
    @abstractmethod
    def get_k(self,h,k):
        '''
        Estimate the weights, k
        '''

class NaiveSolver(K_Solver):
    def __init__(self,atol=1e-8,max_iterations=200):
        '''
        Parameters:
            max_iterations  Maximum number of iterations for solving implicit equations
            atol            Maximum error we can tolerate when solving implicit equations
        '''
        super().__init__(atol=atol,max_iterations=max_iterations)
        
    def get_k(self,h,y):
        '''
        Estimate the weights, k, using the iterative method
        described in Section 2.5.1 of Gear
        '''        
        diff = 0.0
        k = np.zeros((len(self.rki.b),len(y)))
        for _ in range(self.max_iterations):
            k_next = self.iterate(h, y, k)
            diff = np.linalg.norm(k - k_next)
            if diff < self.atol:
                return k_next
            else:
                k = k_next
        a_max = np.abs(self.rki.a).sum(axis=1).max()  
        L = np.abs(self.rki.dy.get_Jacobian(y)).max()
        raise ImplicitRungeKutta.Failed(
            f'Failed to converge within {self.atol} after {self.max_iterations} iterations: error={diff}\n'
            f'a_max={a_max}\n'
            f'L={L}\n'
            f'1/aL={1/(a_max*L)}'
        )
    
    def iterate(self, h, y, k):
        '''
        Iterate One step in solution of iterative equations for k
        
        Parameters:
            h       Step size
            y       Current value of solution
            k       The k term in Runge Kutta
        '''
        k_next= np.zeros_like(k)
        for i in range(self.rki.s):
            k_next[i,:] = self.rki.dy(y + h*np.inner(self.rki.a[i,:],k.T))
        return k_next   
 
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

    def __init__(self, dy, order, a,b,c,solver=None):
        '''
        Initialize
        
        Parameters:
           dy              Function f for ODE = dy/dx=f(y)
           order           Order of solver
           a               The matrix from the Butcher tableau https://mathworld.wolfram.com/ButcherTableau.html
           b               The b vector from the Butcher tableau
           c               The c vector from the Butcher tableau. This is included for completeness,
                           but it is not used, as dy is independent of time.
        '''
        self.dy = dy
        self.order = order
        self.a = a
        self.b = b
        self.c = c
        self.s = len(b)
        self.solver = solver
        
    def __str__(self):
        '''
        Used to display details of integrator
        '''
        return f'{type(self).__name__}: order={self.order}'

    def attach(self,solver):
        '''
        Attach a solver to integrator
        '''
        self.solver = solver
        solver.rki = self
        
    def step(self, h, y):
        '''
        Compute y after next step
        
        Parameters:
           h    Step size
           y    Current value of y
        '''
        return y + h*np.inner(self.b,self.solver.get_k(h,y).T)       
  

class ImplicitRungeKutta2(ImplicitRungeKutta):
    '''
    4th order Gauss-Legendre
    '''
    def __init__(self, dy, solver=None):
        '''
        Parameters:
           dy              Function f for ODE = dy/dx=f(y)
           max_iterations  Maximum number of iterations for solving implicit equations
           atol            Maximum error we can tolerate when solving implicit equations
        '''
        super().__init__(dy, 4,
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
                ]),
                solver=solver)
   
class ImplicitRungeKutta4(ImplicitRungeKutta):
    '''
    6th order Gauss-Legendre
    '''

    def __init__(self, dy, solver=None):
        '''
        Parameters:
           dy              Function f for ODE = dy/dx=f(y)
           max_iterations  Maximum number of iterations for solving implicit equations
           atol            Maximum error we can tolerate when solving implicit equations
        '''
        super().__init__(dy,  6,
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
                            ]),
                         solver=solver)

def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
    return parser.parse_args()

def run_once(ax,rk,nn = 100):
    rk.attach(NaiveSolver())
    try:
        h = 0.15
        y = np.zeros((nn+1,2))
        y[0,0] = 1
        for i in range(nn):
            y[i+1,:] = rk.step(h,y[i,:])      
        ax.plot(y[:,0], y[:,1])
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
             rk=ImplicitRungeKutta2(lambda y: [y[1], -y[0]]),
             nn=150)
    run_once(ax=fig.add_subplot(1,2,2,adjustable='box',aspect=1.0),
             rk=ImplicitRungeKutta4(lambda y: [y[1], -y[0]]),
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
