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
import numpy as np
from scipy.optimize import newton

def fn(r, mu):
    '''
    We wanto to solve for fn()==0 - Murray & Dermott (3.740
    '''
    return (3*r**3*(1 - r + r**2 / 3) / 
            ((1 + r + r**2) * (1 - r)**3)
            - mu/(1 - mu))


def solve(r0, r1, mu, f0=None, f1=None, depth=500, eps=1e-12):
    '''
    My equation solver
    '''
    r = 0.5 * (r0 + r1)
    f = fn(r, mu)
    print(depth, r0, r1, r1 - r0, r, f)

    if depth == 0 or r1 - r0 < eps:
        return r
    if f0 == None:
        f0 = fn(r0, mu)
        f1 = fn(r1, mu)

    s0 = np.copysign(1.0, f0)
    s1 = np.copysign(1.0, f1)
    s = np.copysign(1.0, f)
    if s0 == s:
        return solve(r, r1, mu, f, f1, depth=depth - 1)
    if s == s1:
        return solve(r0, r, mu, f0, f, depth=depth - 1)


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
    alpha = get_alpha(0.2)
    mu = 0.2
    r = newton(lambda x:fn(x, mu),alpha-0.1,x1=alpha+0.1,tol=1e-12)
    r1 = solve(alpha - 0.1, alpha + 0.05, 0.2)
    print (f'r={r},fn(r, mu)={fn(r, mu)},u={get_u(r, 0.2)}')
    print (f'r={r1},fn(r, mu)={fn(r1, mu)},u={get_u(r1, 0.2)}')
  

main()
