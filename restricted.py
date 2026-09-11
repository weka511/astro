#!/usr/bin/env python

# Copyright (C) 2015-2017,2026 Greenweaves Software Pty Ltd

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
Restricted 3 body problem 
'''

from argparse import ArgumentParser
from pathlib import Path
from time import time
import numpy as np
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
from rki import ImplicitRungeKutta4, Driver

__version__ = '1.0'
__author__ = 'Simon Crase'


class ForceCalculator:
    def __init__(self,G, M1, M2):
        self.G = G
        self.M1 = M1
        self.M2 = M2
        
    def get_denominator(self,x1, x2, y1, y2):
        '''
        Used to calculate the denominator in the equation for calculating the
        gravitational attraction between two particles
        '''
        return ((x1 - x2)**2 + (y1 - y2)**2)**(3.0/2.0)
    
    
    def __call__(self,x):
        '''
        Calculate force
        '''
        x1 = x[0]
        y1 = x[1]
        x2 = x[2]
        y2 = x[3]
        x3 = x[4]
        y3 = x[5]
        denom12 = self.get_denominator(x1, x2, y1, y2)
        denom31 = self.get_denominator(x3, x1, y3, y1)
        denom23 = self.get_denominator(x2, x3, y2, y3)
        return np.array([
            x[6],    #x1'
            x[7],    #y1'
            x[8],    #x2'
            x[9],    #y2'
            x[10],   #x3'
            x[11],   #y3'
            - self.G * self.M2 * (x1 - x2) / denom12,     #dx1'
            - self.G * self.M2 * (y1 - y2) / denom12,     #dy1'
            - self.G * self.M1 * (x2 - x1) / denom12,     #dx2'
            - self.G * self.M1 * (y2 - y1) / denom12,     #dy2'
            - self.G * self.M1 * (x3 - x1) / denom31 - self.G * self.M2 * (x3 - x2) / denom23, #dx3'
            - self.G * self.M1 * (y3 - y1) / denom31 - self.G * self.M2 * (y3 - y2) / denom23, #dy3'
        ])


def T(x, G, M1, M2):
    '''
    Kinetic energy
    '''
    return 0.5*(M1 * (x[6]**2 + x[7]**2) + M2 * (x[8]**2 + x[9]**2))  # FIXME - not used at present


def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--show', default=False, action='store_true', help='Used to display figure')
    parser.add_argument('-N', type=int, default=1000)
    parser.add_argument('-M1', type=float, default=250000)
    parser.add_argument('-M2', type=float, default=1.0)
    parser.add_argument('-G', type=float, default=1.0)
    parser.add_argument('-L', type=float, default=25.0) 
    parser.add_argument('--epsilon', type=float, default=0.1) 
    return parser.parse_args()

def create_model(M1=250000,M2=1.0,G=1.0,L=25.0,epsilon = 0.1):
    R2 = L * M1 / (M1 + M2)
    x2 = R2
    y2 = 0
    x1 = -L * M2 / (M1 + M2)
    y1 = 0.0
    x3 = 0.5 * (x1 + x2) * (1 + epsilon)
    y3 = L * np.sqrt(3.0) / 2.0
    omega = np.sqrt(G * (M1 + M2) / (L * L * L))
    return np.array([
            x1, y1,
            x2, y2,
            x3, y3,
            0, x1 * omega,
            0, x2 * omega,
            -0.5 * np.sqrt(3) * L * omega, 0.5 * L * omega
        ])

def main():
    rcParams['text.usetex'] = True
    start = time()
    args = parse_args()

    XY = np.zeros((args.N+1, 12))
    XY[0,:]  = create_model(G=args.G,M1=args.M1,M2=args.M2,L=args.L,epsilon=args.epsilon)
    calculator = ForceCalculator(args.G,args.M1,args.M2)
    driver = Driver(ImplicitRungeKutta4(calculator,
                                        max_iterations=200, atol=1e-18), 
                    h_minimum=1.e-8, h=0.5, h_maximum=1.0, atol=1e-12)

    for i in range(args.N):
        XY[i+1, :] = driver.step(XY[i,:])

    fig = figure(figsize=(8, 8))
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.plot(XY[:, 0], XY[:, 1], 'b', label=f'M1({args.M1})')
    ax1.plot(XY[:, 2], XY[:, 3], 'r', label=f'M2({args.M2})')
    ax1.plot(XY[:, 4], XY[:, 5], 'g', label=f'm({args.epsilon})')

    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_title(f'M1={args.M1:5.1f},M2={args.M2:5.0f},L={args.L:5.1f},N={args.N:,}')
    ax1.legend()

    fig.savefig(Path(args.figs) / Path(__file__).stem)
    elapsed = time() - start
    minutes = int(elapsed / 60)
    seconds = elapsed - 60 * minutes
    print(f'Elapsed Time {minutes} m {seconds:.2f} s')
    if args.show:
        show()


if __name__ == '__main__':
    main()
