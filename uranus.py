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
Create a number of random sets of orbital periods orbital periods for a model 
satellite system, similar to the three inner satellites of Uranus. For each set, 
calculate the mean motion for each stellite and check n1-3n2+2n3.
'''

from argparse import ArgumentParser
from pathlib import Path
from matplotlib.pyplot import figure, show
import numpy as np


def create_orbital_periods(T1=1.413, L=1.546, U=2.101, n=4, rng=np.random.default_rng()):
    '''
    Create a number of random sets of orbital periods orbital periods
    '''
    product = np.zeros((n))
    product[0] = T1
    for i in range(1, n):
        x = rng.uniform()
        ratio = L + x * (U - L)
        product[i] = ratio * product[i - 1]
    return product


def parse_args():
    parser = ArgumentParser(description='Calculate probability of delta_n<0.1 for problem 1.3.')
    parser.add_argument('-N', type=int, help='Number of sets of orbital periods', default=100000)
    parser.add_argument('--seed', '-s', help='Seed for random number generator', default=None)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    return parser.parse_args()


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    delta_n = []
    for orbits in [create_orbital_periods(rng=rng) for i in range(args.N)]:
        n = [360 / T for T in orbits]
        delta_n.append(abs(n[0] - 3 * n[1] + 2 * n[2]))
    P = len([delta for delta in delta_n if delta < 0.1]) / args.N
    bins = 10**(np.arange(-2, 4, dtype=float))
    fig = figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xscale('log')
    ax.hist(delta_n, bins=bins)
    ax.set_title(f'Probability: {P:.2f} after {args.N:,} iterations')
    fig.savefig(Path(args.figs) / Path(__file__).stem)
    show()


if __name__ == '__main__':
    main()

