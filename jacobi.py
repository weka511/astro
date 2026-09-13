#!/usr/bin/env python


# Copyright (C) 2016-2026 Greenweaves Software Pty Ltd

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

'''Plot Zero velocity Surfaces for the Jacobi Integral'''

from argparse import ArgumentParser
import math
import sys
from pathlib import Path
from time import time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
from matplotlib.pyplot import cm
import numpy as np

__version__ = '1.0'
__author__ = 'Simon Crase'


def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--show', default=False, action='store_true', help='Used to display figure')
    return parser.parse_args()


@np.vectorize
def jacobi(x, y, n=1, mu2=0.2, Cj=0):
    n2 = n * n
    mu1 = 1 - mu2
    x2 = x * x
    y2 = y * y
    r1 = np.sqrt((x + mu2) * (x + mu2) + y2)
    r2 = np.sqrt((x - mu1) * (x - mu1) + y2)
    return n2 * (x2 + y2) + 2 * (mu1 / r1 + mu2 / r2) - Cj


def bounds(Z):
    minZ = float('inf')
    maxZ = - minZ
    for zz in Z:
        for z in zz:
            if z < minZ:
                minZ = z
            if z > maxZ:
                maxZ = z
    return (math.floor(minZ), math.ceil(maxZ))


def plot_jacobi(fig, n=1, mu2=0.2, Cj=3.9, limit=5, origin='lower'):
    ax = fig.add_subplot(1, 1, 1)
    xlist = np.linspace(-limit, limit + 0.001, 100)
    ylist = np.linspace(-limit, limit + 0.001, 100)
    X, Y = np.meshgrid(xlist, ylist)
    Z = jacobi(X, Y, n, mu2, Cj)
    (z0, z1) = bounds(Z)

    levels = list(range(z0, 0, 10)) + list(range(0, z1 + 1, 10))

    ticks = [z0, 0, z1]
    c = ax.pcolormesh(X, Y, Z)
    cbar = fig.colorbar(c, orientation='vertical', ticks=ticks)
    cbar.ax.set_yticklabels(['min', '0', 'max'])
    CS3 = ax.contourf(X, Y, Z, levels, cmap=cm.jet, origin=origin)
    CS2 = ax.contour(X, Y, Z, levels=[0], colors='w', origin=origin, hold='on', linewidths=(1,))
    cbar.add_lines(CS2)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(r'Zero velocity surfaces for $n={0},\mu_2 = {1},C_j={2}$'.format(n, mu2, Cj))


def main():
    args = parse_args()
    rcParams['text.usetex'] = True
    start = time()
    n = 1
    for (Cj, limit) in zip([3.80465327630637, 3.552, 3.197, 2.84011], [1.5, 2, 2, 1]):
        fig = figure()
        plot_jacobi(fig=fig, Cj=Cj, limit=limit)
        fig.savefig(f'{Path(args.figs) / Path(__file__).stem}-{n}')
        n += 1
    elapsed = time() - start
    minutes = int(elapsed / 60)
    seconds = elapsed - 60 * minutes
    print(f'Elapsed Time {minutes} m {seconds:.2f} s')
    if args.show:
        show()


if __name__ == '__main__':
    main()
