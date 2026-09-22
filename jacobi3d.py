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

from argparse import ArgumentParser
from pathlib import Path
from time import time
import numpy as np
from matplotlib.pyplot import figure,show,cm
import matplotlib.colors as clrs
from matplotlib import rcParams
import jacobi

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

def create(limit=2, eps=0.001, minZ=-3.9, maxZ=-2.84, steps=1000):
    xlist = np.linspace(-limit, limit + eps, steps)
    ylist = np.linspace(-limit, limit + eps, steps)
    X, Y = np.meshgrid(xlist, ylist)
    Z = -jacobi.jacobi(X, Y)
    Z[Z < minZ] = np.nan
    Z[Z > maxZ] = np.nan
    return X,Y,Z

def plot_3d(file_name,X,Y,Z,limit=2, eps=0.001, minZ=-3.9, maxZ=-2.84):
    '''
    Parameters:
        limit
        eps
        minZ
        maxZ
        steps
    '''
    fig = figure()
    ax = fig.add_subplot(111, projection='3d')    
    surf = ax.plot_surface(X, Y, Z, cmap=cm.jet, norm=clrs.Normalize(vmin=minZ, vmax=maxZ, clip=False))
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    fig.savefig(file_name)

def main():
    args = parse_args()
    rcParams['text.usetex'] = True
    start = time()    
    X,Y,Z = create(limit=1)
    plot_3d(f'{Path(args.figs) / Path(__file__).stem}',X,Y,Z,limit=1)
    elapsed = time() - start
    minutes = int(elapsed / 60)
    seconds = elapsed - 60 * minutes
    print(f'Elapsed Time {minutes} m {seconds:.2f} s')
    if args.show:
        show()
if __name__ == '__main__':
    main()
