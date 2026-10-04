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

'''
Monte Carlo simulation for Murray & Dermott 1.6 (c) and (d)
'''

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
    parser.add_argument('-I','--Max_imax',type=int, default=7)
    return parser.parse_args()

class Logger: 
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

def purge_duplicates(ratios):
    result = set()
    for a,b in ratios:
        factor = 2
        while factor <= a:
            if a%factor == 0 and b%factor == 0:
                a //= factor
                b //= factor
            factor += 1
        result.add((a,b))
                   
    return list(result)

def calculate_ratios(Max_imax):
    ratios = [(1,2)]
    Nrs = []
    i_maxen = []
    Nr_naive = []
    for imax in range(3,Max_imax+1):
        i_maxen.append(imax)
        Nr_naive.append(imax*(imax-1)//2)
        for i in range(1,imax):
            ratios.append((i,imax))
        ratios = sorted(purge_duplicates(ratios),key=lambda ratio:ratio[1]*imax+ratio[0])
        Nrs.append((len(ratios)))
        Logger.instance.debug(ratios)
        
    return i_maxen,Nr_naive,Nrs

def main():
    '''
    Monte Carlo simulation for Murray & Dermott 1.6 (c) and (d)
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    Logger.create(f'{Path(args.logs)/Path(__file__).stem}{strftime('%Y%m%d%H%M%S')}.log')
    
    i_maxen,Nr_naive,Nrs = calculate_ratios(args.Max_imax)
        
    fig = figure(figsize=(12,12))
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(1,1,1,adjustable='box',aspect=1.0)
    ax1.scatter(i_maxen,Nr_naive,marker='x')
    ax1.scatter(i_maxen,Nrs,marker='+')
    fig.tight_layout(h_pad=2)
    fig.savefig(Path(args.figs)/Path(__file__).stem)    
    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
    if args.show:
        show()
    
if __name__=='__main__':
    main()
