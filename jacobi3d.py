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

'''Plot Potential surfaces for the Jacobi Integral'''

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
    parser.add_argument('--limit',default=1,type=float,help='Used to set upper and lower bounds for X and Y')
    parser.add_argument('--minZ',default=-3.9,type=float,help='Used to normalize Z (lower limit)')
    parser.add_argument('--maxZ',default=-2.8,type=float,help='Used to normalize Z (upper limit)')
    parser.add_argument('--steps',default=1000,type=int,help='Number of steps when data generated')
    return parser.parse_args()

def create_surface(limit=2, minZ=-3.9, maxZ=-2.84, steps=1000):
    '''
    Create X, Y, and Z cordinates for surfaces
    
    Parameters:
        limit   Used to set upper and lower bounds for X and Y 
        minZ    Used to clip and normalize Z (lower limit)
        maxZ    Used to clip and normalize Z (upper limit)
        steps   Number of steps when data generated
    '''    
    X, Y = np.meshgrid(np.linspace(-limit, limit, steps), 
                       np.linspace(-limit, limit, steps))
    Z = -jacobi.jacobi(X, Y)
    Z[Z < minZ] = np.nan
    Z[Z > maxZ] = np.nan
    return X,Y,Z

def plot_3d(file_name,X,Y,Z,limit=1,  minZ=-3.9, maxZ=-2.84,heading=''):
    '''
    Plot the surfaces
    
    Parameters:
        file_name   Name of file where figure will be stored
        X           X coordinates for data to be plotted
        Y           Y coordinates for data to be plotted
        Z           Z coordinates for data to be plotted
        limit       Used to set upper and lower bounds for X and Y 
        minZ        Used to normalize Z (lower limit)
        maxZ        Used to normalize Z (upper limit)
    '''
    fig = figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection='3d')      
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    fig.colorbar(ax.plot_surface(X, Y, Z, 
                                 cmap=cm.jet, 
                                 norm=clrs.Normalize(vmin=minZ, 
                                                     vmax=maxZ,
                                                     clip=False)), 
                 shrink=0.5,
                 aspect=5)
    fig.suptitle(heading)
    fig.savefig(file_name)

def main():
    '''
    Calculate and plot Potential surfaces for the Jacobi Integral
    '''
    args = parse_args()
    rcParams['text.usetex'] = True
    start = time()
    
    X,Y,Z = create_surface(limit=args.limit,
                           maxZ=args.maxZ,
                           minZ=args.minZ,
                           steps=args.steps)
    
    plot_3d(f'{Path(args.figs) / Path(__file__).stem}',X,Y,Z,
            limit=args.limit,
            maxZ=args.maxZ,
            minZ=args.minZ,
            heading=f'{ Path(__file__).stem}: limit={args.limit}, maxZ={args.maxZ}, minZ={args.minZ}')
    
    elapsed = time() - start
    minutes = int(elapsed / 60)
    seconds = elapsed - 60 * minutes
    print(f'Elapsed Time {minutes} m {seconds:.2f} s')
    
    if args.show:
        show()
        
if __name__ == '__main__':
    main()
