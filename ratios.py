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
Exercise 1.6 (e) Use the data in Appendix A to find the periods of all possible pairs of
periods among the planates and the prograde salellites of Mars, Jupiter, Saturn, Uranus, 
and Neptune, with mean radii > 100 km and orbital eccentricites < 0.1. Taking i_max =7, show
that thirty pairs of objects have ratios of orbital periods witin epsilon-max of a permitted
commensurability.
'''

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

class SimpleRatios:
    def __init__(self,imax=7):
        self.imax = imax
        self.pairs = sorted([(i1,i2) for i2 in range(1,imax+1) for i1 in range(1,i2)],key=lambda x:x[0]/x[1])
        self.ratios = [i1/i2 for i1,i2 in self.pairs]
        differences = [a-b for a in self.ratios for b in self.ratios if a > b]
        self.eps_max = 0.5*min(differences)
        
    def get_match(self,ratio):
        best = np.argmin(abs(ratio-self.ratios))
        return self.pairs[best],abs(ratio-self.ratios[best])
    
    def get_sequence(self,i1,i2):
        return [i for i in range(len(self.pairs)) if (i1,i2) == self.pairs[i]][0]
    
class Resonance:
    
    @staticmethod
    def create(names,periods):
        return [Resonance(names[i],names[j],periods[i],periods[j]) 
                          for i in range(len(names)) 
                          for j in range(i+1,len(names))]
    
    def __init__(self,name1,name2,T1,T2):
        self.name1 = name1
        self.name2 = name2
        self.T1 = T1
        self.T2 = T2
        self.i1 = None
        self.i2 = None
        
    def can_match(self,simple_ratios):
        best_pair,distance = simple_ratios.get_match(self.T1/self.T2)
        if distance < simple_ratios.eps_max:
            self.i1,self.i2 = best_pair
            return True
        return False
    
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

 

def get_periods(path,key='Planet',e_max=0.15,R_min=0):
    '''
    Retrieve periods for planets and satellites. If the column for the period has no data,
    use Kepler's Third Law to calculate periods
    
    Parameters:
        path      Path to planetary data
    '''
    df = pd.read_csv(path)
    df1 = df[df['e'] < e_max]
    if R_min > 0:                            #FIXME
        df1 = df1[df1['R'] >= R_min]
    try:
        return df1[key].to_list(),df1['T'].to_numpy()
    except KeyError:
        return df1[key].to_list(),df1['a'].to_numpy()**(3/2)


def main():
    '''
    Use the data in Appendix A to find the periods of all possible pairs of
    periods among the planates and the prograde salellites of Mars, Jupiter, Saturn, Uranus, 
    and Neptune, with mean radii > 100 km and orbital eccentricites < 0.1. Taking i_max =7, show
    that thirty pairs of objects have ratios of orbital periods witin epsilon-max of a permitted
    commensurability.
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    simple_ratios = SimpleRatios()
    names,periods = get_periods((Path(args.data)/args.planets).with_suffix('.csv'),e_max=args.e_max)
       
    resonances = Resonance.create(names,periods)
    
    for primary in args.names:
        names,periods = get_periods((Path(args.data)/primary).with_suffix('.csv'),
                                    key='Satellite',
                                    e_max=args.e_max,R_min=args.R_min)
        
        resonances += Resonance.create(names,periods)
    
    resonances1 = [resonance for resonance in resonances if resonance.can_match(simple_ratios)]
    
    colours = [
        'xkcd:purple','xkcd:green','xkcd:blue','xkcd:pink','xkcd:brown','xkcd:red',
        'xkcd:light blue','xkcd:teal','xkcd:orange','xkcd:light green','xkcd:magenta','xkcd:yellow',
        'xkcd:sky blue','xkcd:grey','xkcd:lime green','xkcd:light purple','xkcd:violet','xkcd:dark green',
        'xkcd:turquoise','xkcd:lavender','xkcd:dark blue','xkcd:tan','xkcd:cyan','xkcd:aqua',
        'xkcd:forest green','xkcd:mauve','xkcd:dark purple','xkcd:bright green','xkcd:maroon','xkcd:olive'
    ]    
    fig = figure(figsize=(12,12))
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(1,1,1)
    for r in resonances1:
        ax1.scatter(r.i1,r.i2,
                    label=f'{r.name1}-{r.name2} {r.i2}:{r.i1}',
                    c=colours[simple_ratios.get_sequence(r.i1,r.i2)])
    ax1.legend(ncols=max(1,len(resonances1)//12))
    
    #ax1.bar(observed_ratios,5,
            #width=0.1,
            #facecolor=colours,
            #edgecolor='xkcd:white'
    #)
    #ax2 = fig.add_subplot(2,1,2)
    #ax2.bar(matched_ratios,5,
            #width=0.1,
            #facecolor=colours,
            #edgecolor='xkcd:white'
    #)
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
