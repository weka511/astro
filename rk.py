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

'''Explicit Runge Kutta'''

from abc import ABC,abstractmethod
from argparse import ArgumentParser
from pathlib import Path
from time import strftime,time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np

__version__ = '1.0'
__author__ = 'Simon Crase'

class RungeKutta(ABC):
    '''
    Explicit Runge Kutta
    
    https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Explicit_methods
    
    '''
    def __init__(self, dx, order, a,b,c):
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
        np.testing.assert_array_almost_equal(c,np.sum(a,axis=1),decimal=16)
        self.dx = dx
        self.order = order
        self.a = a
        self.b = b
        self.c = c
        
    def __str__(self):
        '''
        Used to display details of integrator
        '''
        return f'{type(self).__name__}: order={self.order}'
    
    def step(self, h, x):
        '''
        Compute x after next step
        
        Parameters:
           h    Step size
           x    Current value of x
        '''
        k = np.zeros((len(self.b),len(x)))
        k[0,:] = self.dx(x)
        for i in range(1,len(self.b)):
            k[i,:] = self.dx(x + h*sum([self.a[i,j]*k[j,:] for j in range(i)]))
            
        return (x + h*np.inner(self.b,k.T)) 

class Rk4(RungeKutta):
    def __init__(self, dx):
        '''
        Parameters:
           dy              Function f for ODE = dy/dx=f(y)
           max_iterations  Maximum number of iterations for solving implicit equations
           atol            Maximum error we can tolerate when solving implicit equations
        '''
        super().__init__(dx, 4,
                np.array([
                    [0,   0,   0,  0],
                    [1/2, 0,   0,  0],
                    [0,   1/2, 0, 0],
                    [0,   0,   1, 0],
                ]),
                np.array([
                    1/6, 1/3, 1/3, 1/6
                ]),
                np.array([
                    0, 1/2, 1/2, 1
                ]))    

class Nystrom(RungeKutta):
    '''
    Nyström's fifth order method
    https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Nystr%C3%B6m's_fifth-order_method
    '''
    def __init__(self, dx):
        '''
        Parameters:
           dy              Function f for ODE = dy/dx=f(y)
           max_iterations  Maximum number of iterations for solving implicit equations
           atol            Maximum error we can tolerate when solving implicit equations
        '''
        super().__init__(dx, 5,
                np.array([
                    [0,    0,    0,      0,    0, 0],
                    [1/3,  0,    0,      0,    0, 0],
                    [4/25, 6/25, 0,      0,    0, 0],
                    [1/4,  -3,   15/4,  0,    0, 0],
                    [2/27, 10/9, -50/81, 8/81, 0, 0],
                    [2/25, 12/25, 2/15,  8/75, 0, 0],
                ]),
                np.array([
                    23/192, 0, 125/192, 0, -27/64, 125/192
                ]),
                np.array([
                    0, 1/3, 2/5, 1, 2/3, 4/5
                ]))    
        
def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--logs',default='./logs', help=f'Path to log files')
    parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
    return parser.parse_args()


def main():
    '''
    Solve equation for simple harmonic motion, and plot
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
  
    rk4 = Nystrom(lambda y: [y[1], -y[0]])
    h = 0.01
    nn = 1000
    y = np.zeros((nn+1,2))
    y[0,0] = 1
    for i in range(nn):
        y[i+1,:] = rk4.step(h,y[i,:])       
    fig = figure(figsize=(12,12))
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(1,1,1,adjustable='box',aspect=1.0)
    ax1.plot(y[:,0], y[:,1])
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

