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
from math import floor,ceil
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
    '''
    Evaluate the Jacobi Entegral --Murray & Dermott (3.29)
    '''
    mu1 = 1 - mu2
    r1 = np.sqrt((x + mu2)**2 + y**2)
    r2 = np.sqrt((x - mu1)**2 + y**2)
    return n**2 * (x**2 + y**2) + 2 * (mu1 / r1 + mu2 / r2) - Cj

def plot_jacobi(fig, n=1, mu2=0.2, Cj=3.9, limit=5, origin='lower'):
    X, Y = np.meshgrid(np.linspace(-limit, limit + 0.001, 100), 
                       np.linspace(-limit, limit + 0.001, 100))
    Z = jacobi(X, Y, n, mu2, Cj)
    z0 = floor(Z.min())
    z1 = ceil(Z.max())
    levels = list(range(z0, 0, 10)) + list(range(0, z1 + 1, 10))
    ticks = [z0, 0, z1]
    
    ax = fig.add_subplot(1, 1, 1)
    
    cbar = fig.colorbar(ax.pcolormesh(X, Y, Z), 
                        orientation='vertical', 
                        ticks=ticks)
    cbar.ax.set_yticklabels(['min', '0', 'max'])
    CS3 = ax.contourf(X, Y, Z, levels, cmap=cm.jet, origin=origin)
    CS2 = ax.contour(X, Y, Z, levels=[0], colors='w', origin=origin, linewidths=(1,))
    cbar.add_lines(CS2)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(r'Zero velocity surfaces for $n={0},\mu_2 = {1},C_j={2}$'.format(n, mu2, Cj))

def main():
    args = parse_args()
    rcParams['text.usetex'] = True
    start = time()
    Params = np.array([
                    [3.80465327630637,1.5],
                    [3.552,2],
                    [3.197,2],
                    [2.84011,1]
    ])
    m,_ = Params.shape
    for i in range(m):
        fig = figure()
        plot_jacobi(fig=fig, Cj=Params[i,0], limit=Params[i,1])
        fig.savefig(f'{Path(args.figs) / Path(__file__).stem}-{i+1}')
  
    elapsed = time() - start
    minutes = int(elapsed / 60)
    seconds = elapsed - 60 * minutes
    print(f'Elapsed Time {minutes} m {seconds:.2f} s')
    if args.show:
        show()

if __name__ == '__main__':
    main()
