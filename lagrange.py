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
            Smaller of the two main masses
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
            Smaller of the two main masses
        '''
        mu1 = 1 - mu2
        return (3*r2**3*
                (1 + r2 + r2**2/3) / 
                ((1 + r2)**2 * (1 - r2**3)) 
                - mu2/mu1)    
        
class L3(LagrangianPoint):
    def __init__(self,mu2):
        super().__init__(mu2)  
        
    def fn(r1):
        '''
        We want to  solve  Murray & Dermott (3.91)
        
        Parameters:
            Smaller of the two main masses
        '''
        return ((1 - r1**3)*(1 + r1)**2/
                (r1**3*(r1**2 + 3*r1 + 3))
                - self.mu2/self.mu1)

    
def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--mu2', default=0.2, type=float,help='Smaller of the two main masses')
    return parser.parse_args()

def main():
    args = parse_args()
    
    solver1 = L1(args.mu2)
    alpha = solver1.get_alpha()
    r2,err = solver1.solve(alpha)
    jacobi1 = 2*solver1.get_u(1-r2,r2)
    print (f'L1: r2={r2},err={err},Jacobi={jacobi1}')
    
    #start = alpha + alpha**2/3 - alpha**3/9 - 31*alpha**4/81
    #L2 = newton(lambda x:fn_L2(x, args.mu2),start-0.1,x1=start+0.1,tol=1e-12)
    #print (f'L2={L2},fn(r2, mu2)={fn_L2(L2, args.mu2)},Jacobi={2*get_u(L2,args.mu2)}') 

    solver3 = L3(args.mu2)
    beta = -(7/12)*solver3.mu2/solver3.mu1 + (7/12)*(solver3.mu2/solver3.mu1)**2 - (13223/20736)*(solver3.mu2/solver3.mu1)**3
    r1,err = solver1.solve(beta + 1,bound=0.01)
    jacobi3 = 2*solver3.get_u(r1,r1-1)
    z=0
    #L3 = newton(lambda x:fn_L3(x, args.mu2),beta-0.9,x1=beta+0.9,tol=1e-12)
    #print (f'L3={L3},fn(r1, mu2)={fn_L3(L3, args.mu2)},Jacobi={2*get_u(L3,args.mu2)}')     
  
main()
