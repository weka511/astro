#!/usr/bin/env python

#   Copyright (C) 2026 Simon Crase

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

'''Adams Bashforth Predictor with Adams Moulton Corrector'''

from argparse import ArgumentParser
from csv import reader
from logging import basicConfig, getLogger, INFO, FileHandler, StreamHandler, Formatter, DEBUG
from pathlib import Path
from time import strftime, time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np
from rk import Rk4

__version__ = '1.0'
__author__ = 'Simon Crase'


class Predictor:
    @staticmethod
    def get_coefficients(order:int):
        beta = [
            [1],
            [3, -1],
            [23, -16, 5],
            [55, -59, 37, -9],
            [1901, -2774, 2616, -1274, 251],
            [4277, -7923, 9982, -7298, 2877, -475]
        ]
        sums = [1, 2, 12, 24,720,1440]
        return np.array(beta[order - 1]),sums[order - 1]
        

    def __init__(self, dx, order, starter):
        self.dx = dx
        self.beta,self.divisor = Predictor.get_coefficients(order)
        assert self.divisor == self.beta.sum()
        self.order = order 
        self.state = 0
        self.DX = None
        self.starter = starter

    def __str__(self):
        '''
        Used to display details of integrator
        '''
        return f'{type(self).__name__}: order={self.order}'

    def step(self, h, x):
        '''
        Compute x after next step
        
        Parameters:
           h    Step size
           x    Current value of x
        '''
        if self.state == 0:
            n, = x.shape
            self.X = np.zeros((self.order, n))
            self.X[self.order - 1 - self.state, :] = x
            self.state += 1
            self.X[self.order - 1 - self.state, :] = self.starter.step(h, x)
            self.state += 1
        elif self.state < self.order - 1:
            self.state += 1
            self.X[self.order - 1 - self.state, :] = self.starter.step(h, x)
        else:
            pass
        return self.X[0, :]


class Corrector:
    @staticmethod
    def get_coefficients(order:int):
        beta = [
            [1],
            [1, 1],
            [5, 8, -1],
            [9, 19, -5, 1],
            [251, 646, -264, 106, -19],
            [475, 1427, -798, 482, -173, 27]
        ]
        sums = [1, 2, 12, 24,720,1440]
        return np.array(beta[order - 1]),sums[order - 1]  
    
    def __init__(self, dx, beta, predictor):
        self.dx = dx
        self.beta = beta
        self.order, = beta.shape

    def __str__(self):
        '''
        Used to display details of integrator
        '''
        return f'{type(self).__name__}: order={self.order}'

def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--logs', default='./logs', help=f'Path to log files')
    parser.add_argument('--show', default=False, action='store_true', help='Used to display figure')
    parser.add_argument('--order',default=6,type=int)
    return parser.parse_args()


class Logger:
    '''
    Singleton Logger
    '''
    instance = None

    @staticmethod
    def create(path_name: str):
        '''
        Set up console logger and file logger
        
        Parameters:
            path_name    Path name for log file
        '''
        def add_handler(handler, logger, formatter=Formatter('%(message)s'), level: int = INFO):
            handler.setLevel(level)
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        Logger.instance = getLogger(__name__)
        Logger.instance.setLevel(DEBUG)
        add_handler(FileHandler(path_name), Logger.instance, level=DEBUG)
        add_handler(StreamHandler(), Logger.instance)
        np.set_printoptions(linewidth=np.nan) # Prevent lines being split when we log numpy arrays


def main():
    '''
     Solve equation for simple harmonic motion, and plot
    '''
    rcParams['text.usetex'] = True
    start = time()
    args = parse_args()
    def dx(x): return [x[1], -x[0]]
    starter = Rk4(dx)
    predictor = Predictor(dx,args.order,starter)
    Corrector.get_coefficients(6)
    h = 0.01
    nn = 1000
    y = np.zeros((nn + 1, 2))
    y[0, 0] = 1
    for i in range(nn):
        y[i + 1, :] = predictor.step(h, y[i, :])
    fig = figure(figsize=(12, 12))
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(1, 1, 1, adjustable='box', aspect=1.0)
    ax1.plot(y[:, 0], y[:, 1])
    fig.tight_layout(h_pad=2)
    fig.savefig(Path(args.figs) / Path(__file__).stem)
    elapsed = time() - start
    minutes = int(elapsed / 60)
    seconds = elapsed - 60 * minutes
    print(f'Elapsed Time {minutes} m {seconds:.2f} s')
    if args.show:
        show()


if __name__ == '__main__':
    #try:
    main()
    #except Exception as e:
        #Logger.instance.exception(e,stack_info=True)
