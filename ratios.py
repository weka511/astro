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
    parser.add_argument('names',nargs='+')
    parser.add_argument('--data', default='./data', help=f'Path to data files')
    return parser.parse_args()

def read_and_filter(path,D_min=100,e_max=0.15):
    df = pd.read_csv(path)
    df1 = df[(df['T'] > 0) & (df['R'] > 0.5*D_min) & (df['e'] < e_max)]
    return df1['T'].to_numpy()
 
def create_ratios(T):
    m, = T.shape
    product = np.zeros((m*(m-1))//2)
    k = 0
    for i in range(m):
        for j in range(i):
            product[k] = T[i]/T[j]
            k += 1
    return product

def main():
    '''
    Do whatever...
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    ratios = []
    for name in args.names:
        T = read_and_filter((Path(args.data)/name).with_suffix('.csv'))
        if len(T) < 2: continue
        ratios.append(create_ratios(T))
    all_ratios= np.hstack(ratios)
    all_ratios = np.sort(all_ratios)
    fig = figure(figsize=(12,12))
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(1,1,1,adjustable='box',aspect=1.0)
    ax1.hist(all_ratios,bins=100)
    #ax1.bar(all_ratios,10,
            #color=[
                #'xkcd:purple','xkcd:green','xkcd:blue', 'xkcd:pink','xkcd:brown','xkcd:red',
                #'xkcd:light blue','xkcd:teal','xkcd:orange', 'xkcd:light green','xkcd:magenta','xkcd:yellow',
                #'xkcd:sky blue','xkcd:grey','xkcd:lime green', 'xkcd:light purple','xkcd:violet','xkcd:dark green',
            #])
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
