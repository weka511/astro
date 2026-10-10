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
Calculate Equilibrium points L1, L2, and L3, plus Jacobi constants, 
for Lagrange configuration of the 3 body problem
'''

from abc import ABC,abstractmethod
from argparse import ArgumentParser
from pathlib import Path
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
from matplotlib.pyplot import cm
import numpy as np
from scipy.optimize import newton


__version__ = '1.0'
__author__ = 'Simon Crase'

class LagrangianPoint:
    def __init__(self,mu2):
        self.mu1 = 1 - mu2
        self.mu2 = mu2
        
    def get_alpha(self):
        '''
        This function is used to inintialize r:
        Murray & Dermott, (3.75)
        
        Parameters:
            mu2   Smaller of the two main masses
            r2    Distance of smallest mass from mu2
        '''
        return (self.mu2/(3*self.mu1))**(1/3)    
        
    @abstractmethod
    def fn(self,r):
        ...
    
    def solve(self,alpha,bound=0.1,tol=1e-12):
        r = newton(self.fn,alpha-bound,x1=alpha+bound,tol=tol)
        err = self.fn(r)
        return r,err
    
    def get_u(self,r1,r2):
        '''
        Calculate pseudo potential
        
        Murray & Dermott, (3.64)
        
        Parameters:
            mu2   Smaller of the two main masses
            r2    Distance of smallest mass from mu2
        '''  
        return self.mu1*(1/r1 + r1**2/2) + self.mu2*(1/r2 + r2**2/2) - self.mu1*self.mu2/2    
        
class L1(LagrangianPoint):
    def __init__(self,mu2):
        super().__init__(mu2)
        
    def fn(self, r2):
        '''
        We want to  solve  Murray & Dermott (3.74)
        
        Parameters:
            r2  Distance of smallest mass from mu2
        '''
        return (3*r2**3*
                (1 - r2 + r2**2/3) / 
                ((1 + r2 + r2**2) * (1 - r2)**3) 
                - self.mu2/self.mu1)
    
class L2(LagrangianPoint):
    def __init__(self,mu2):
        super().__init__(mu2)
        
    def fn(self,r2):
        '''
        We want to  solve  Murray & Dermott (3.86)
        
        Parameters:
            r2    Distance of smallest mass from mu2
        '''
        return (3*r2**3*
                (1 + r2 + r2**2/3) / 
                ((1 + r2)**2 * (1 - r2**3)) 
                - self.mu2/self.mu1)    
        
class L3(LagrangianPoint):
    def __init__(self,mu2):
        super().__init__(mu2)  
        
    def fn(self,r1):
        '''
        We want to  solve  Murray & Dermott (3.91)
        
        Parameters:
            r1   Distance of smallest mass from mu1
        '''
        return ((1 - r1**3)*(1 + r1)**2/
                (r1**3*(r1**2 + 3*r1 + 3))
                - self.mu2/self.mu1)

    
def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--mu2', default=0.2, type=float,help='Smaller of the two main masses')
    parser.add_argument('--show',default=False,action='store_true')
    parser.add_argument('--limit',default=1.0,type=float)
    parser.add_argument('--step',default=100,type=int)
    return parser.parse_args()

@np.vectorize
def jacobi(x,y,point=None, Cj=0):
    r1 = np.sqrt((x+point.mu2)**2 + y**2)
    r2 = np.sqrt((x-point.mu1)**2 + y**2)
    return point.get_u(r1,r2)

def main():
    args = parse_args()
    
    solver1 = L1(args.mu2)
    alpha = solver1.get_alpha()
    r2_1,err1 = solver1.solve(alpha)
    jacobi1 = 2*solver1.get_u(1-r2_1,r2_1)
    #print (f'L1: r2={r2},err={err},Jacobi={jacobi1}')
    
    solver2 = L2(args.mu2)
    r2_2,err2 = solver2.solve(alpha)
    jacobi2 = 2*solver2.get_u(1-r2_2,r2_2)  
    
    solver3 = L3(args.mu2)
    beta = -(7/12)*solver3.mu2/solver3.mu1 + (7/12)*(solver3.mu2/solver3.mu1)**2 - (13223/20736)*(solver3.mu2/solver3.mu1)**3
    r1_3,err3 = solver3.solve(beta + 1,bound=0.01)       
    jacobi3 = 2*solver3.get_u(r1_3,1-r1_3)
    
    X, Y = np.meshgrid(np.linspace(-args.limit, args.limit + 1/args.step, args.step), 
                       np.linspace(-args.limit, args.limit + 1/args.step, args.step))
    Z = jacobi(X, Y, point=solver1, Cj=jacobi1)
    z0 = int(np.floor(Z.min()))
    z1 = int(np.ceil(Z.max()))
    levels = list(range(z0, z1, 10))
    fig = figure(figsize=(8,8))
    ax = fig.add_subplot(1,1,1)    
    CS3 = ax.contourf(X, Y, Z, levels, cmap=cm.viridis, origin=None)
    CS2 = ax.contour(X, Y, Z, levels=[0], colors='w', origin=None, linewidths=(1,))
    cbar = fig.colorbar(ax.pcolormesh(X, Y, Z), 
                        orientation='vertical', 
                        ticks=None)

    cbar.add_lines(CS2)

    ax.scatter(r2_1,0,marker='+',label=f'L1: r2={r2_1:3f},err={err1:.3g},Jacobi={jacobi1:3f}',c='xkcd:red')
    ax.scatter(r2_2,0,marker='x',label=f'L2: r2={r2_2:3f},err={err2:.3g},Jacobi={jacobi2:3f}',c='xkcd:cyan')
    ax.scatter(r1_3,0,marker='X',label=f'L2: r2={r1_3:3f},err={err3:.3g},Jacobi={jacobi3:3f}',c='xkcd:bright green')
    
    ax.legend(title='Lagrange points')
      
    fig.tight_layout(h_pad=2)
    fig.savefig(Path(args.figs)/Path(__file__).stem)    
    
    if args.show:
        show()
  
main()
