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

'''Restricted 3 body problem after Kotovych and Bowman'''

from argparse import ArgumentParser
from csv import reader
from pathlib import Path
from time import time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np
from integrators import Hamiltonian, KotovychBowman

__version__ = '1.0'
__author__ = 'Simon Crase'

class Restricted3Body(Hamiltonian):
    '''
    Hamiltonian for retricted 3 body problem
    '''
    def __init__(self,mu):
        self.mu = mu
        x = 1
        y = 0
        x_dot = 0
        y_dot = 1
        self.X = np.array([x,y,x_dot,y_dot])
        
    def dx(self):
        '''
        Calculate derivatives 
        '''    
     
    def d_xi(self):
        pass
    
    def transform(self):
        '''
        Transform to the coordinates we will use for integration
        '''
        
    def _invert(self, hamiltonian):
        pass
    
    def get_energy(self):
        '''
        Calculate total energy    equation (12)
        '''
        return 0.5*(self.X[2]**2 + self.X[3]**2) + self.V()
    
    def V(self):
        r1 = np.sqrt((self.X[0]-self.mu)**2 + self.X[1]**2)
        r2 = np.sqrt((self.X[0]+1-self.mu)**2 + self.X[1]**2)
        return -0.5*(self.X[0]**2 + self.X[1]**2) - (1-self.mu)/r1 -self.mu/r2
        
    def dV(self, r, theta, rho, Theta):
        '''
        Calculate derivatives of total energy 
        '''   
    
    def dH(self):
        pass
    
    def inverse_jacobi(self):
        pass
    
def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
    parser.add_argument('--mu', type=float,default=0.01)
    return parser.parse_args()
    
def main():
    '''
    Restricted 3 body problem after Kotovych and Bowman
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    fig = figure(figsize=(12,12))
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(1,1,1,adjustable='box',aspect=1.0)
    hamiltonian = Restricted3Body(args.mu)
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
