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

'''Template for python script'''

from argparse import ArgumentParser
from csv import reader
from pathlib import Path
from time import time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np
from rki import ImplicitRungeKutta4,Driver

__version__ = '1.0'
__author__ = 'Simon Crase'

def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('file_name',help='File name for initial conditions')
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
    parser.add_argument('--data', default='./data', help=f'Path to data files')
    return parser.parse_args()

def read_data(file_name,dim=2):
    '''
    Read initial conditions from a file
    '''
    state = -1
    i = 0
    with open(file_name) as f:
        data_reader = reader(f)
        for row in data_reader:
            values = np.array([float(m) for m in row])
            match(state):
                case -1:
                    masses = values
                    Q = np.zeros((len(masses),dim))
                    P = np.zeros((len(masses),dim))
                    state = 0
                case 0:
                    Q[i,:] = values
                    i += 1
                    if i == len(masses):
                        i = 0
                        state = 1
                case 1:
                    P[i,:] = values
                    i += 1
 
        return Q,P,masses
 
class Hamiltonian:
    r = 0
    theta = 1
    p = 2
    l = 3
    rho = 4
    Theta = 5
    P = 6
    L = 7
    def __init__(self,m,G=1,clone=False,atol=1e-16):
        self.G = G
        self.M = m.sum()        # Section 5
        self.mu = m[0] + m[1]   # Section 5
        self.m = m
        self.g1 = m[0]*m[1]/self.mu    # Reduced mass - Section 5 - just after (28) 
        self.g2 = m[2]*self.mu/self.M  # Reduced mass - Section 5 - just after (28) 
        
    def dH(self,y):  # r theta p l R Theta P L
        dV =self.dV(y)
        return np.array([
            y[Hamiltonian.p]/self.g1,
            y[Hamiltonian.l]/(self.g1*y[Hamiltonian.p]**2),
            y[Hamiltonian.l]**2 / (self.g1*self.r**3) - dV[Hamiltonian.r],
            - dV[Hamiltonian.theta],
            y[Hamiltonian.P]/self.g1,
            y[Hamiltonian.L]/(self.g1*y[Hamiltonan.rho]**2),
            y[Hamiltonian.L]**2 / (self.g1*self.rho**3) - dV[Hamiltonian.rho],
            - dV[Hamiltonian.Theta],            
        ])
    
    def dV(self,y):
        dV = np.zeros((8))
        dV[Hamiltonian.r] = 0
        dv[Hamiltonian.theta] = 0
        dV[Hamiltonian.rho] = 0
        dv[Hamiltonian.Theta] = 0        
        return dV
        
def main():
    '''
    Do whatever...
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    
    R,R_dot,m = read_data(Path(args.data)/args.file_name)
    hamiltonian = Hamiltonian(m)    
    integrator = ImplicitRungeKutta4(lambda y: hamiltonian.dH, 10, 0.000000001)
    
    fig = figure(figsize=(12,12))
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(1,1,1,adjustable='box',aspect=1.0)
    
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
