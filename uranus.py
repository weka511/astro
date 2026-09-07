#!/usr/bin/env python

# Copyright (C) 2019-2026 Simon Crase

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

import argparse
import matplotlib.pyplot as plt
import numpy as np

def generate_orbit(T1=1.413, L=1.546, U=2.101, range_i=range(1,3),rng = np.random.default_rng()):
    Ts = [T1]
    for i in range_i:
        x = rng.uniform()
        ratio = L + x*(U-L)
        Ts.append(ratio * Ts[-1])
    return Ts

def parse_args():
    parser = argparse.ArgumentParser(description='Calculate probability of delta_n<0.1 for problem 1.3.')    
    parser.add_argument('-N',type=int, help='Number of sets of orbital periods',default=100000)
    parser.add_argument('--seed','-s',help='Seed for random number generator',default=None)
    return parser.parse_args()

def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    delta_n = []
    for orbits in  [generate_orbit(rng=rng) for i in range(args.N)]:
        n = [360/T for T in orbits]
        delta_n.append(abs(n[0]-3*n[1] + 2*n[2]))
    bins=10**(np.arange(-2,4,dtype=float))
    plt.xscale('log')
    plt.hist(delta_n,bins=bins)
    plt.title('Probability: {0:.2f}% after {1:,} iterations'.format(100*len([delta for delta in delta_n if delta<0.1])/100000,
                                         args.N))
    plt.show()
    
if __name__=='__main__':
    main()
            