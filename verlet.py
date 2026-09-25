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
Verlet integrator
'''

from argparse import ArgumentParser
from csv import reader
from pathlib import Path
from time import time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np

__version__ = '1.0'
__author__ = 'Simon Crase'

class VelocityVerlet:
    def __init__(self,hamiltonian):
        self.hamiltonian = hamiltonian
        self.F_half = None
        self.first_step = True
        self.n = len(hamiltonian)
        
    def __str__(self):
        '''
        Used to display details of integrator
        '''
        return f'{type(self).__name__}'    
        
    def step(self, h, y):
        '''
        Compute y after next step
        
        Parameters:
           h    Step size
           y    Current value of y
        '''
        q = y[0:self.n]
        p = y[self.n:]
        if self.first_step == True:
            self.F_half = self.hamiltonian.dp(y)
            self.first_step = False
        p_half = p - 0.5*h*self.F_half
        q += h*p_half          # m?
        self.F_half = self.hamiltonian.dp(y)
        p = p_half - 0.5*h*self.F_half
  
        return np.hstack([q,p])
        
class SHM:
    '''
    Test using simple harmonc motion
    '''
    def __init__(self,m=1.0,k=1.0):
        self.m = m
        self.k = k
        
    def __len__(self):
        return 1
        
    def get_total_energy(self,p,q):
        return 0.5 * (p**2/self.m + self.h * self.m*q**2)
    
    def dq(self,y):
        q = y[0:len(self)]
        p = y[len(self):]
        return p/self.m
    
    def dp(self,y):
        q = y[0:len(self)]
        p = y[len(self):]
        return self.k * self.m *q
    
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
    y = np.array([1,0],dtype=float)
    X[0,:] = y
    for i in range(args.N):
        y = integrator.step(0.1,y)
        X[i+1,:] = y
    
    fig = figure(figsize=(8,8))
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(1,1,1,adjustable='box',aspect=1.0)
    ax1.scatter(X[:,0],X[:,1],c='xkcd:blue',s=1)
    ax1.scatter(X[0,0],X[0,1],c='xkcd:blue',marker='X',s=25)
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
