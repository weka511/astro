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

from argparse import ArgumentParser
import numpy as np
from scipy.optimize import newton

__version__ = '1.0'
__author__ = 'Simon Crase'

def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--mu2', default=0.2, type=float,help='Smaller of the two main masses')
    return parser.parse_args()

def fn_L1(r2, mu2):
    '''
    We want to  solve  Murray & Dermott (3.74)
    
    Parameters:
        Smaller of the two main masses
    '''
    mu1 = 1 - mu2
    return (3*r2**3*
            (1 - r2 + r2**2 / 3) / 
            ((1 + r2 + r2**2) * (1 - r2)**3) 
            - mu2/mu1)

def fn_L2(r2, mu2):
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

def fn_L3(r1, mu2):
    '''
    We want to  solve  Murray & Dermott (3.91)
    
    Parameters:
        Smaller of the two main masses
    '''
    mu1 = 1 - mu2
    return ((1-r1**3)*(1+r1)**2/
            (r1**3*(r1**2+3*r1+3))
            - mu2/mu1)

def get_alpha(mu2):
    '''
    This function is used to inintialize r:
    Murray & Dermott, (3.75)
    
    Parameters:
        mu2   Smaller of the two main masses
        r2    Distance of smallest mass from mu2
    '''
    mu1 = 1 - mu2
    return (mu2/(3*mu1))**(1 / 3)

def get_u(r2, mu2):
    '''
    Calculate pseudo potential
    
    Murray & Dermott, (3.64)
    
    Parameters:
        mu2   Smaller of the two main masses
        r2    Distance of smallest mass from mu2
    '''    
    mu1 = 1 - mu2
    r1 = 1 - r2
    return mu1 * (1/r1 + r1**2/ 2) + mu2 * (1/r2 + r2**2/2) - mu1*mu2/2

def main():
    args = parse_args()
    alpha = get_alpha(args.mu2)
    L1 = newton(lambda x:fn_L1(x, args.mu2),alpha-0.1,x1=alpha+0.1,tol=1e-12)
    print (f'L1={L1},fn(r2, mu2)={fn_L1(L1, args.mu2)},Jacobi={2*get_u(L1,args.mu2)}')
    L2 = newton(lambda x:fn_L2(x, args.mu2),alpha-0.1,x1=alpha+0.1,tol=1e-12)
    print (f'L2={L2},fn(r2, mu2)={fn_L2(L2, args.mu2)},Jacobi={2*get_u(L2,args.mu2)}') 
    beta = -(7/12) *args.mu2/(1-args.mu2) +(7/12) * (args.mu2/(1-args.mu2))**2 - (13223/20736)* (args.mu2/(1-args.mu2))**3
    L3 = newton(lambda x:fn_L3(x, args.mu2),beta-0.9,x1=beta+0.9,tol=1e-12)
    print (f'L3={L3},fn(r1, mu2)={fn_L3(L3, args.mu2)},Jacobi={2*get_u(L3,args.mu2)}')     
  
main()
