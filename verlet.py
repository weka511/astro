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
Velocity Verlet integrator, based on: The leapfrog method and other "symplectic" 
algorithms for integrating Newton's laws of motion, Peter Young, April 14 2009,
https://courses.physics.ucsd.edu/2010/Winter/physics141/Assignments/leapfrog_py.pdf
'''

from argparse import ArgumentParser
from csv import reader
from pathlib import Path
from time import time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np

__version__ = '1.1'
__author__ = 'Simon Crase'

class VelocityVerlet:
    '''
    Velocity Verlet Integrator, as described in Young, Section III
    
    Attributes:
        F 
        F_saved
        is_first_step
        n
    '''
    def __init__(self,F):
        '''
        Parameters:
            F
        '''
        self.F = F
        self.F_half = None
        self.is_first_step = True
        
    def __str__(self):
        '''
        Used to display details of integrator
        '''
        return f'{type(self).__name__}'    
        
    def step(self, h, x):
        '''
        Compute x after next step
        
        Parameters:
           h    Step size
           x    Current value of x
        '''
        n = len(x)//2
        # Calculate v at the mid point of the step
        if self.is_first_step:
            self.F_saved = self.F(x)
            self.is_first_step = False
            
        v_mid_point = x[n:] + 0.5*h*self.F_saved[n:]
 
        # Calculate x at the end of the step       
        x_end = x.copy()
        x_end[0:n] +=  h*v_mid_point
        
        # Calculate x at the end of the step 
        for _ in range(2):
            self.F_saved = self.F(x_end)
            x_end[n:] = v_mid_point + 0.5*h*self.F_saved[n:]
        return x_end
 
    
class SHM:
    '''
    Test using simple harmonc motion
    '''
    def __init__(self,m=1.0,k=1.0):
        self.m = m
        self.k = k
        
    def __len__(self):
        return 1
    
    def __call__(self,x):
        dx = np.zeros_like(x)
        dx[0] = x[1]
        dx[1] = - self.k*x[0]/self.m
        return dx
    
    def get_total_energy(self,y):
        return 0.5 * (self.m*y[1]**2 + self.k * self.m*y[0]**2)
    
    
def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('-N','--N',default=628,type=int)
    parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
    return parser.parse_args()
    
def main():
    '''
    Use Verlet integrator to plot SHM
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    shm = SHM()
    integrator = VelocityVerlet(shm)
    X = np.zeros((args.N+1,2))
    E = np.zeros((args.N+1))
    y = np.array([1,0],dtype=float)
    X[0,:] = y
    E[0] = shm.get_total_energy(y)
    for i in range(args.N):
        y = integrator.step(0.1,y)
        X[i+1,:] = y
        E[i+1] = shm.get_total_energy(y)
    
    fig = figure(figsize=(12,6))
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(1,2,1,adjustable='box',aspect=1.0)
    ax1.scatter(X[:,0],X[:,1],c='xkcd:blue',s=1)
    ax1.scatter(X[0,0],X[0,1],c='xkcd:blue',marker='X',s=25)
    ax2 = fig.add_subplot(1,2,2,adjustable='box',aspect=1.0)
    ax2.plot(E)
    fig.tight_layout(h_pad=2)
    fig.savefig(Path(args.figs)/Path(__file__).stem)    
    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
    if args.show:
        show()
    
if __name__=='__main__':
    main()
