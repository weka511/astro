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
Calculate Equilibrium points L1, L2, and L3, of Lagrange configuration of the 3 body problem
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
    parser.add_argument('--mu', default=0.2, type=float)
    return parser.parse_args()

def fn(r, mu):
    '''
    We wanto to solve for fn()==0 - Murray & Dermott (3.740
    '''
    return (3*r**3*(1 - r + r**2 / 3) / 
            ((1 + r + r**2) * (1 - r)**3)
            - mu/(1 - mu))


def get_alpha(mu):
    '''
    Murray & Dermott, (3.75)
    '''
    return (mu/(3*(1 - mu)))**(1 / 3)

def get_u(r, mu):
    '''
    Murray & Dermott, (3.64)
    '''    
    mu1 = 1 - mu
    r1 = 1 - r
    return mu1 * (1/r1 + r1**2/ 2) + mu * (1/r + r**2/2) - mu1*mu/2

def main():
    args = parse_args()
    alpha = get_alpha(args.mu)
    r = newton(lambda x:fn(x, args.mu),alpha-0.1,x1=alpha+0.1,tol=1e-12)
    print (f'r={r},fn(r, mu)={fn(r, args.mu)},u={get_u(r, args.mu)}')
  
main()
