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
from matplotlib.ticker import MaxNLocator
import numpy as np
from ratios import Ratios,Logger

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
 
            
def main():
    '''
    Monte Carlo simulation for Murray & Dermott 1.6 (c) and (d)
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    Logger.create(f'{Path(args.logs)/Path(__file__).stem}{strftime('%Y%m%d%H%M%S')}.log')
    
    ratios = Ratios()
    ratios.build(args.Max_imax)
        
    fig = figure(figsize=(12,12))
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(2,2,1)
    ax1.scatter(ratios.i_maxen,ratios.Nr_naive,marker='x',label='Naive')
    ax1.scatter(ratios.i_maxen,ratios.Nrs,marker='+',label='Smart')
    ax1.set_ylabel('$N_r$')
    ax1.set_xlabel('$i_{max}$')
    ax1.set_ylim((0,max(ratios.Nr_naive)+1))
    ax1.set_xlim((0,args.Max_imax+1))    
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.yaxis.set_major_locator(MaxNLocator(integer=True)) 
    ax1.legend(loc='upper left')
    
    ax2 = fig.add_subplot(2,2,2)
    ax2.scatter(ratios.i_maxen,ratios.eps_max,marker='x')
    ax2.set_ylabel(r'$\epsilon_{max}$')
    ax2.set_xlabel('$i_{max}$')
    ax2.set_ylim((0,1.05*max(ratios.eps_max)))
    ax2.set_xlim((0,args.Max_imax+1))    
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
   
    
    fig.savefig(Path(args.figs)/Path(__file__).stem)    
    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
    if args.show:
        show()
    
if __name__=='__main__':
    main()
