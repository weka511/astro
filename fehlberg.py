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

'''Explore use of Fehlberg 78'''

from argparse import ArgumentParser
from csv import reader
from logging import basicConfig,getLogger,INFO,FileHandler,StreamHandler,Formatter,DEBUG
from pathlib import Path
from time import strftime,time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np
from odespy import rkexplicit,abm4
from restricted import Synodic

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
    parser.add_argument('--fn', choices=['abm','rk'],default='rk')
    parser.add_argument('--method',choices=['rk45','rk78','rkf45','rkf78','rkck','rkv'],default='rkf78')
    parser.add_argument('--step',type=float,default=0.01,help='Stepsize for integrator')
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
        
def vdp_func(t, y, p):
    dy = np.zeros((len(y), ))
    dy[0] = y[1]
    dy[1] = p[0]*(1 - y[0]**2) * y[1] - y[0]
    return dy

def main():
    '''
    Evolve orbit using Fehlberg 78 and plot
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    Logger.create(f'{Path(args.logs)/Path(__file__).stem}{strftime('%Y%m%d%H%M%S')}.log')
    synodic = Synodic(0.000953875)
    
    trange = [0, 100.0]
    yinit = [-1.02745, 0,0,0.04032]
    params = [1]
    
    match args.fn:
        case 'abm':
            ts, ys, lte= abm4(lambda t,y,p: synodic(y), trange, yinit, params, args.step,
                                            nord='inf')
            Logger.instance.info(f'lte={lte}')
        case rk:
            ts, ys, yhat, stats = rkexplicit(lambda t,y,p: synodic(y), trange, yinit, params, 
                                             method=args.method, nord='inf',abstol=1e-12)    
    
            Logger.instance.info(stats)
    
    n = ys.shape[1] 
    fig = figure(figsize=(12,12))
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(1,1,1)
    ax1.plot(ys[0,:],ys[1,:])
    match args.fn:
        case 'abm':    
            ax1.set_title(f'{args.fn}, step={args.step}')
        case 'rk':    
            ax1.set_title(f'{args.fn} {args.method}')        
    fig.tight_layout(h_pad=2,pad=2)
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
