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
from typing import Callable
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np
import numpy.typing as npt
from rk import Rk4

__version__ = '1.0'
__author__ = 'Simon Crase'


class Predictor:
    '''
    This class implements the Adams-Bashforth predictor
    
    Attributes:
        dx
        beta 
        divisor
        order
        state
        DX
        starter
        filled
    '''
    @staticmethod
    def get_coefficients(order:int):
        '''
        Calculate the numerators for the coefficients
        '''
        beta = [
            [1],
            [3, -1],
            [23, -16, 5],
            [55, -59, 37, -9],
            [1901, -2774, 2616, -1274, 251],
            [4277, -7923, 9982, -7298, 2877, -475]
        ]
        return np.array(beta[order-1])
        
    def __init__(self, dx:Callable[[npt.NDArray[float]], npt.NDArray[float]], order:int, starter):
        self.dx = dx
        self.beta = Predictor.get_coefficients(order)
        self.divisor = self.beta.sum()
        self.order = order 
        self.state = -1
        self.DX = np.zeros((order+1,0))
        self.starter = starter
        self.filled = False
        
    def __str__(self):
        '''
        Used to display details of integrator
        '''
        return f'{type(self).__name__}: order={self.order}'

    def step(self, h:float, x:npt.NDArray[float]):
        '''
        Compute x after next step
        
        Parameters:
           h    Step size
           x    Current value of x
        '''
        if self.state == -1:
            n, = x.shape
            self.DX = np.full((self.order+1, n),np.nan)
            self.state += 1
            
        if self.state < self.order:        
            self.DX[self.state, :] = self.dx(x)
            self.state += 1    
            return self.starter.step(h,x)              
        else:
            x_new = x + h* np.dot(self.beta,self.DX[:-1])/self.divisor
            self.DX[-1,:] = self.dx(x_new)
            self.filled = True
            return x_new
        
    def roll(self):
        '''
        Used after a step to make room for the latest fx(x)
        '''
        if self.filled:
            self.DX = np.roll(self.DX,-1,axis=0)

class Corrector:
    '''
    This class implements the Adams-Moulton corrector
    
    Attributes:
        beta 
        divisor
        predictor
        n
    '''
    @staticmethod
    def get_coefficients(order:int):
        '''
        Calculate the numerators for the coefficients
        '''        
        beta = [
            [1],
            [1, 1],
            [5, 8, -1],
            [9, 19, -5, 1],
            [251, 646, -264, 106, -19],
            [475, 1427, -798, 482, -173, 27]
        ]
        return np.array(beta[order-1])
    
    def __init__(self, predictor,n:int=1):
        self.beta = Corrector.get_coefficients(predictor.order)
        self.divisor = self.beta.sum()
        self.predictor = predictor
        self.n = n
        
    def __str__(self):
        '''
        Used to display details of integrator
        '''
        return f'{type(self).__name__}: order={self.order}'
    
    def step(self, h:float, x:npt.NDArray[float]):
        '''
        Compute x after next step
        
        Parameters:
           h    Step size
           x    Current value of x
        '''
        prediction = self.predictor.step(h,x)
        if self.predictor.filled:
            for i in range(self.n):
                x_new = x + h* np.dot(self.beta,self.predictor.DX[:-1])/self.divisor
                self.predictor.DX[-1,:] = self.predictor.dx(x_new)
            self.predictor.roll()
            return x_new
        else:
            return prediction

def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--logs', default='./logs', help=f'Path to log files')
    parser.add_argument('--show', default=False, action='store_true', help='Used to display figure')
    parser.add_argument('--order',default=6,type=int)
    parser.add_argument('--N',default=1,type=int)
    parser.add_argument('--step',default=0.001,type=float)
    parser.add_argument('--n',default=1,type=int)
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

def dx(x): return [x[1], -x[0]]

def main():
    '''
     Solve equation for simple harmonic motion, and plot
    '''
    rcParams['text.usetex'] = True
    start = time()
    args = parse_args()
    starter = Rk4(dx)
    corrector = Corrector(Predictor(dx,args.order,starter),
                          n=args.n)
    h = args.step
    nn = int(2*args.N*np.pi/h)
    y = np.zeros((nn + 1, 2))
    y[0, 0] = 1
    for i in range(nn):
        y[i + 1, :] = corrector.step(h, y[i, :])
       
    fig = figure(figsize=(12, 12))
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(1, 1, 1, adjustable='box', aspect=1.0)
    ax1.scatter(y[0,0],y[0,1],marker='+')
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
