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

'''Template for python script'''

from argparse import ArgumentParser
from csv import reader
from pathlib import Path
from sys import float_info
from time import time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np
import pandas as pd

__version__ = '1.0'
__author__ = 'Simon Crase'

        
def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
    parser.add_argument('names',nargs='*')
    parser.add_argument('--data', default='./data', help=f'Path to data files')
    parser.add_argument('--planets', default='planets')
    parser.add_argument('--R_min', default=100, type=float)
    parser.add_argument('--e_max', default=0.15, type=float)
    return parser.parse_args()

def read_and_filter(path,R_min=100,e_max=0.15):
    df = pd.read_csv(path)
    df1 = df[(df['T'] > 0) & (df['R'] > R_min) & (df['e'] < e_max)]
    return df1['T'].to_numpy()
 
def create_ratios(T):
    m, = T.shape
    product = np.zeros((m*(m-1))//2)
    k = 0
    for i in range(m):
        for j in range(i):
            product[k] = T[i]/T[j] if T[i]<T[j] else T[j]/T[i]
            k += 1
    return product

def get_periods(path,e_max=0.15):
    '''
    Retrieve periods for planets. If the column for the period has no data,
    use Kepler's Third Law to calulate periods
    
    Parameters:
        path      Path to planetary data
    '''
    df = pd.read_csv(path)
    df1 = df[df['e'] < e_max]
    try:
        return df1['T']/to_numpy()
    except KeyError:
        return df1['a'].to_numpy()**(3/2)

class SimpleRatios:
    def __init__(self,imax=7):
        self.imax = imax
        self.pairs = sorted([(i1,i2) for i2 in range(1,imax+1) for i1 in range(1,i2)],key=lambda x:x[0]/x[1])
        self.ratios = [i1/i2 for i1,i2 in self.pairs]
        differences = [a-b for a in self.ratios for b in self.ratios if a > b]
        self.eps_max = 0.5*min(differences)
        

    def get_match(self,ratio):
        best = None
        distance = float_info.max
        for i in range(len(self.pairs)):
            if abs(ratio-self.ratios[i]) < distance:
                best = i
                distance = abs(ratio-self.ratios[i])
        return self.pairs[best],distance
                
  

def create_observed_ratios(args):
    ratios = []
    ratios.append(
        create_ratios(
            get_periods((Path(args.data)/args.planets).with_suffix('.csv'),
                        e_max=args.e_max)))

    for name in args.names:
        T = read_and_filter((Path(args.data)/name).with_suffix('.csv'),
                            R_min=args.R_min,
                            e_max=args.e_max)
        if len(T) < 2: continue
        ratios.append(create_ratios(T))
 
    return np.sort(np.hstack(ratios))

def main():
    '''
    Do whatever...
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    simple_ratios = SimpleRatios()
    observed_ratios = create_observed_ratios(args)
    matches = []
    for ratio in observed_ratios:
        pair,distance = simple_ratios.get_match(ratio)
        if distance < simple_ratios.eps_max:
            matches.append(pair)
    matched_ratios = [i1/i2 for i1,i2 in matches]  
    
    colours = [
        'xkcd:purple','xkcd:green','xkcd:blue','xkcd:pink','xkcd:brown','xkcd:red',
        'xkcd:light blue','xkcd:teal','xkcd:orange','xkcd:light green','xkcd:magenta','xkcd:yellow',
        'xkcd:sky blue','xkcd:grey','xkcd:lime green','xkcd:light purple','xkcd:violet','xkcd:dark green',
        'xkcd:turquoise','xkcd:lavender','xkcd:dark blue','xkcd:tan','xkcd:cyan','xkcd:aqua',
        'xkcd:forest green','xkcd:mauve','xkcd:dark purple','xkcd:bright green','xkcd:maroon','xkcd:olive'
    ]    
    fig = figure(figsize=(6,6))
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(2,1,1)
    ax1.bar(observed_ratios,5,
            width=0.1,
            facecolor=colours,
            edgecolor='xkcd:white'
    )
    ax2 = fig.add_subplot(2,1,2)
    ax2.bar(matched_ratios,5,
            width=0.1,
            facecolor=colours,
            edgecolor='xkcd:white'
    )
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
