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

'''Plot U'''

from argparse import ArgumentParser
from csv import reader
from logging import basicConfig,getLogger,INFO,FileHandler,StreamHandler,Formatter,DEBUG
from pathlib import Path
from time import strftime,time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np


__version__ = '1.0'
__author__ = 'Simon Crase'

def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--logs',default='./logs', help=f'Path to log files')
    parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
    parser.add_argument('--mu2', default=0.2, type=float,help='Smaller of the two main masses')
    parser.add_argument('--limit',default=1.0,type=float)
    parser.add_argument('--step',default=0.0625,type=int)    
    return parser.parse_args()


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

def get_U(x,y,mu2=0.2):
    mu1 = 1 - mu2
    r1 = np.sqrt((x+mu2)**2 + y**2)                                       # (3.8)
    r2 = np.sqrt((x-mu1)**2 + y**2)                                       # (3.9)
    U = 2*(mu1*(1/r1 + r1**2/2) + mu2*(1/r2 + r2**2/2) - mu1*mu2/2)    # (3.64)
    return U

def main():
    '''
    Compute U and plot it
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    Logger.create(f'{Path(args.logs)/Path(__file__).stem}{strftime('%Y%m%d%H%M%S')}.log')
    jacobi = np.vectorize(lambda x,y:get_U(x,y,mu2=args.mu2))
    
    X = np.arange(-args.limit,args.limit,args.step)
    Y = np.arange(-args.limit,args.limit,args.step)
    X,Y = np.meshgrid(X,Y)
    Z = np.log(jacobi(X,Y))
    
    fig = figure(figsize=(12,12))
    fig.suptitle(Path(__file__).stem)
    
    ax1 = fig.add_subplot(1,1,1,projection='3d')
    ax1.plot_surface(X,Y,Z,cmap='coolwarm',linewidth=0,antialiased=False)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('-2U')
    
    fig.tight_layout(h_pad=2)
    fig.savefig(Path(args.figs)/Path(__file__).stem)    
    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
    if args.show:
        show()
    
if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        Logger.instance.exception(e,stack_info=True)
