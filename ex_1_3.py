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

'''
Murray & Dermott, Exercise 1.3.
Create a large number of random sets of orbital periods orbital periods for a model 
satellite system, similar to the three inner satellites of Uranus. For each set, 
calculate the mean motion for each satellite and compute |n1-3n2+2n3|. 
Show that the probability that delta_n<0.1 is ~0.3%.
'''

from argparse import ArgumentParser
from pathlib import Path
from matplotlib.pyplot import figure, show
import numpy as np

__version__ = '1.1'
__author__ = 'Simon Crase'

def create_orbital_periods(T1=1.413, L=1.546, U=2.101, k=4, rng=np.random.default_rng()):
    '''
    Create a random set of orbital periods using (1.6) of Murray and Dermott
    
    Parameters:
        T1       Starting value for each set
        L        Lower limit on the ratio of successive orbital periods
        U        Upper limit on the ratio of successive orbital periods
        k        Size of each set
        rng      Generator for random values
    '''
    product = np.zeros((k))
    product[0] = T1
    for i in range(1, k):
        ratio = L + rng.uniform() * (U - L)
        product[i] = ratio * product[i - 1]
    return product

def parse_args():
    parser = ArgumentParser(description='Calculate probability of delta_n<0.1 for problem 1.3.')
    parser.add_argument('-N', type=int, help='Number of sets of orbital periods', default=100000)
    parser.add_argument('--seed', '-s', help='Seed for random number generator', default=None)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots') 
    parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
    return parser.parse_args()

def create_delta_n(N,rng = np.random.default_rng()):
    '''
    Create a large number of random sets of orbital periods orbital periods for a model 
    satellite system, similar to the three inner satellites of Uranus. For each set, 
    calculate the mean motion for each satellite and compute |n1-3n2+2n3|. 
    '''
    delta_n = np.zeros((N))
    weights = np.array([1,-3,2,0])  # Used to calculate: abs(n[0]-3*n[1] + 2*n[3])
    for j,orbital_periods in enumerate([create_orbital_periods(rng=rng) for i in range(N)]):
        n = 360 / orbital_periods
        delta_n[j] = np.abs(np.dot(weights,n))
    return delta_n

def main():
    '''
    Generate orbits, calculate probabailty that delta_n within tolerance, 
    and plot  histogram for delta_n
    '''
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    
    fig = figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xscale('log')
    
    # Draw histogram. This first bin corresponds to the desired interval n<0.1
    n,bins,_ = ax.hist(create_delta_n(args.N,rng=rng), 
                       bins=10**(np.arange(-2, 4, dtype=float)),
                       color='xkcd:blue',
                       label=r'd')
    ax.axvline(bins[1],color='xkcd:red',linestyle=':',label=r'$\delta n=0.1$')
    # Mark area that corresponds to desired values: the weight of 2.25 ia a kludge
    ax.axhline(n[0],color='xkcd:red',xmin=bins[0],xmax=2.25*bins[1],label=r'$\delta n<0.1$')
    ax.set_xlabel(r'$\delta n$')
    ax.set_ylabel('Count')
    ax.legend(loc='upper left')
    ax.set_title(rf'$P[|n_1-3n_2+2n_3|<0.1]\approx$ {100*n[0]/args.N:.2f}% after {args.N:,} iterations')
    for spine in ['top', 'bottom', 'left', 'right']:
        ax.spines[spine].set_color('xkcd:light grey')
        ax.spines[spine].set_linewidth(1)    
    fig.savefig(Path(args.figs) / Path(__file__).stem)
    if args.show:
        show()

if __name__ == '__main__':
    main()
