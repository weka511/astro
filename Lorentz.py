#!/usr/bin/env python

# Copyright (C) 2015-2026 Greenweaves Software Limited

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
Test the ImplicitRungeKutts Integrator by calculating the
evolution of the Lorentz Attractor
'''

from argparse import ArgumentParser
from logging import basicConfig,getLogger,INFO,FileHandler,StreamHandler,Formatter,DEBUG
from pathlib import Path
from time import time,strftime
from warnings import filterwarnings
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from integrators import IntegratorFactory

__version__ = '1.1'
__author__ = 'Simon Crase'


class Lorentz:
    def __init__(self, sigma, rho, beta):
        self.sigma = sigma
        self.rho = rho
        self.beta = beta

    def __call__(self, x):
        return np.array([
            self.sigma * (x[1] - x[0]),
            x[0] * (self.rho - x[2]) - x[1],
            x[0] * x[1] - self.beta * x[2]
        ])

    def __str__(self):
        return 'sigma={0},rho={1}.beta={2}'.format(self.sigma, self.rho, self.beta)

class Logger: 
    '''
    Singleton Logger
    '''
    instance = None
    
    @staticmethod
    def create(path_name:str):
        '''
        Set up console logger and file logger
        
        Parameters:
            path_name    Path name for log file
        '''
        def add_handler(handler,logger,formatter=Formatter('%(message)s'),level:int=INFO):
            handler.setLevel(level)
            handler.setFormatter(formatter)
            logger.addHandler(handler)    
        Logger.instance = getLogger(__name__)
        Logger.instance.setLevel(DEBUG)
        add_handler(FileHandler(path_name),Logger.instance,level=DEBUG)
        add_handler(StreamHandler(),Logger.instance)   
        np.set_printoptions(linewidth=np.nan) # Prevent lines being split when we log numpy arrays
        
def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--show', default=False, action='store_true', help='Used to display figure')
    parser.add_argument('-N', '--N', type=int, default=10000)
    parser.add_argument('--step', type=float, default=0.01)
    parser.add_argument('--integrator',
                        choices=IntegratorFactory.get_choices(),
                        default=IntegratorFactory.get_default(),
                        help='Choice of integrator')
    parser.add_argument('--atol', type=float, default=1.0e-8, help='Tolerance for solving equation in implicit Runge Kutta')
    parser.add_argument('-m', '--max_iterations', type=int, default=250,
                        help='Number of iterations for solving equation in implicit Runge Kutta')
    parser.add_argument('--logs',default='./logs', help=f'Path to log files')
    return parser.parse_args()


def main():
    '''
    Calculate the solution to the Lorentz equations and plot it
    '''
    filterwarnings("error", category=RuntimeWarning)
    rcParams['text.usetex'] = True
    start = time()
    args = parse_args()
    Logger.create(f'{Path(args.logs)/Path(__file__).stem}{strftime('%Y%m%d%H%M%S')}.log')
    lorentz = Lorentz(10, 28, 8.0 / 3.0)
    integrator = IntegratorFactory.create(lorentz, args.integrator, args.atol, args.max_iterations)
    XYZ = np.zeros((args.N + 1, 3))
    XYZ[0, :] = np.array([1, 0, 1])
    try:
        for i in range(args.N):
            XYZ[i + 1, :] = integrator.step(args.step, XYZ[i, :])
    except RuntimeWarning as e:
        Logger.instance.info(f'Runtime warning at i={i}')
        Logger.instance.error(e)
        XYZ = XYZ[0:i,:]

    fig = figure(figsize=(12, 12))
    fig.suptitle(integrator)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(XYZ[:, 0], XYZ[:, 1], XYZ[:, 2], s=1, edgecolors='face')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(lorentz)
    fig.savefig(Path(args.figs) / Path(__file__).stem)
    elapsed = time() - start
    minutes = int(elapsed / 60)
    seconds = elapsed - 60 * minutes
    Logger.instance.info(f'Elapsed Time {minutes} m {seconds:.2f} s')
    if args.show:
        show()


if __name__ == '__main__':
    main()
