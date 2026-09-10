#!/usr/bin/env python

# Copyright (C) 2015-2026 Greenweaves Software Limited

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
Hamiltonian for integrating Kepler problem
'''

from argparse import ArgumentParser
from pathlib import Path
from time import time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np
from integrators import Hamiltonian, Integrate2

__version__ = '1.0'
__author__ = 'Simon Crase'


class Kepler(Hamiltonian):
    '''
    The Hamiltonisan for the Kepler problem
    '''
    def __init__(self,x,m,k):
        self.x = x
        self.m = m
        self.k = k
        self.transform()
        
    def create(self,x):
        return Kepler(x,self.m,self.k)
    
    def transform(self):
        r = self.x[0]
        p = self.x[2]
        L = self.x[3]
        self.eta = [-self.k/r, p*p/(2*self.m)+L*L/(2*self.m*r*r),L]
    
    def invert(self,kepler):
        super(Kepler,self).invert(kepler)  #FIXME
        r = -self.k/self.eta[0]
        L = self.eta[2]
        p_squared = 2*self.m*self.eta[1] - L**2/(r*r)        
        p = np.sign(self.x[2])*np.sqrt(p_squared) if p_squared>0 else 0
        self.x[0] = r
        self.x[2] = p
        self.x[3] = L

    def dx(self):
        r = self.x[0]
        L = self.x[3]
        return [self.x[2]/self.m,L/(self.m*r**2), L**2/(self.m*r**3) - self.k / r**2,0]    

    def d_eta(self):
        term = self.k*self.x[2]/(self.m*self.x[0]**2)
        return [term,-term,0]

    def hamiltonian(self):
        return (self.x[2]**2/(2*self.m) + 
               self.x[3]**2/(2*self.m*self.x[0]**2)-self.k/self.x[0])

    def display(self):
        print ("x",self.x)
        print ("eta",self.eta)

def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
    parser.add_argument('-N',type=int,default=10000)
    parser.add_argument('--step',type=float,default=0.001)
    parser.add_argument('-k',type=float,default=1.0)
    parser.add_argument('-r',type=float,default=1.0)
    parser.add_argument('-p',type=float,default=0.0001)
    parser.add_argument('-m',type=float,default=0.001)
    return parser.parse_args()

            
def main():
    '''
    Integrate Kepler problem
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
 
    L = np.sqrt((args.m/args.r**3))
 
    pos = np.zeros((args.N,2))
    kepler = Kepler( np.array([args.r,0,args.p,L]),args.m,args.k)
    hamiltonian = kepler.hamiltonian()
    integrator  = Integrate2(args.step,kepler)
    
    for i in range(args.N):
        integrator.integrate()
        pos[i,:] = kepler.x[0]*np.array([np.cos(kepler.x[1]),np.sin(kepler.x[1])])
    
    fig = figure(figsize=(8, 8))
    ax1 = fig.add_subplot(1, 1, 1)        
    ax1.plot(pos[:,0],pos[:,1])
    ax1.set_title(rf'$\delta H=${kepler.hamiltonian()-hamiltonian:.4e} after {args.N:,} steps of size {args.step}')
  
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