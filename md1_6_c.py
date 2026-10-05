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

class Ratios:
    def __init__(self):
        self.ratios = [(1,2)]
        self.Nrs = []
        self.i_maxen = []
        self.Nr_naive = []
        self.primes = []
        self.eps_max = []
        
    def build(self,Max_imax):
        primes = create_primes(Max_imax)
        for imax in range(3,Max_imax+1):
            self.i_maxen.append(imax)
            self.Nr_naive.append(imax*(imax-1)//2)
            self._add_ratios(imax,primes)
            self.eps_max.append(0.5/(imax*(imax-1)))
                 
    def _add_ratios(self,imax,primes):
        for i in range(1,imax):
            self.ratios.append((i,imax))
        factors = factorize(imax,primes)
        if len(factors) > 1:
            self.ratios = sorted(self._purge_duplicates(self.ratios),
                                 key=lambda ratio:ratio[1]*imax+ratio[0])
        self.Nrs.append((len(self.ratios)))
        Logger.instance.info(self.ratios)
            
    def _purge_duplicates(self,ratios):
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

def create_primes(N):
    '''
    Create a list of prime numbers using the sieve of Eratosthenes
    
    Parameters:
        N          Largest condidate to be considered
        
    Returns:
        List of primes up to (and, if relevant), N
    '''
    canditates = list(range(2,N+1))
    i = 0
    while i < len(canditates):
        sieved = canditates[0:i+1]
        for j in range(i+1,len(canditates)):
            if canditates[j] % canditates[i] != 0:
                sieved.append(canditates[j])
        canditates = sieved
        i += 1
    return canditates
            
def factorize(n,primes):
    '''
    Factorize a number into a list of primes
    
    Parameters:
        n       Number to be factorized
        primes  List of primes to be considered
        
    Returns:
      List of factors (repeated as many times as necessary)
    '''
    factors = []
    for p in primes:
        while n%p == 0:
            factors.append(p)
            n //= p
    return factors        
            
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
