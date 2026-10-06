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
and Neptune, with mean radii > 100 km and orbital eccentricites < 0.15. Taking i_max =7, show
that thirty pairs of objects have ratios of orbital periods witin epsilon-max of a permitted
commensurability.
'''

from argparse import ArgumentParser
from logging import basicConfig,getLogger,INFO,FileHandler,StreamHandler,Formatter,DEBUG
from pathlib import Path
from sys import float_info
from time import strftime,time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Patch
import numpy as np
from ratios import Ratios,SimpleRatios,Resonance,get_periods

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
        
def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
    parser.add_argument('names',nargs='*',help='File name for satellite data')
    parser.add_argument('--data', default='./data', help=f'Path to data files')
    parser.add_argument('--planets',default='planets',help='File name for planetary data')
    parser.add_argument('--R_min',default=100, type=float,help='We require data for satellites to have R greater than this value')
    parser.add_argument('--e_max',default=0.15, type=float,help='We require data for satellites to have eccentricity less than this value')
    parser.add_argument('--imax',default=7,type=int,help='Limit on denominator for acceptable ratios')
    parser.add_argument('--logs',default='./logs', help=f'Path to log files')
    return parser.parse_args()

def create_xkcd_colours():
    '''
    48 colours from https://blog.xkcd.com/2010/05/03/color-survey-results/
    '''
    return [
        'xkcd:purple','xkcd:green','xkcd:blue','xkcd:pink','xkcd:brown','xkcd:red',
        'xkcd:light blue','xkcd:teal','xkcd:orange','xkcd:light green','xkcd:magenta','xkcd:yellow',
        'xkcd:sky blue','xkcd:grey','xkcd:lime green','xkcd:light purple','xkcd:violet','xkcd:dark green',
        'xkcd:turquoise','xkcd:lavender','xkcd:dark blue','xkcd:tan','xkcd:cyan','xkcd:aqua',
        'xkcd:forest green','xkcd:mauve','xkcd:dark purple','xkcd:bright green','xkcd:maroon','xkcd:olive',
        'xkcd:salmon','xkcd:beige','xkcd:royal blue','xkcd:navy','xkcd:lilac','xkcd:black',
        'xkcd:hot pink','xkcd:light brown','xkcd:pale green','xkcd:peach','xkcd:olive green','xkcd:dark pink',
        'xkcd:periwinkle','xkcd:sea green','xkcd:lime','xkcd:indigo','xkcd:mustard','xkcd:light pink'
    ]  

def frac(s):
    ss = [int(x) for x in s.split(':')]
    return 100*ss[1] + ss[0]

def main():
    '''
    Use the data in Appendix A to find the periods of all possible pairs of
    periods among the planates and the prograde salellites of Mars, Jupiter, Saturn, Uranus, 
    and Neptune, with mean radii > 100 km and orbital eccentricites < 0.1. Taking i_max =7, show
    that thirty pairs of objects have ratios of orbital periods within epsilon-max of a permitted
    commensurability.
    '''
    rcParams['text.usetex'] = True
    rcParams['figure.constrained_layout.use'] = True
    start  = time()
    args = parse_args()
    Logger.create(f'{Path(args.logs)/Path(__file__).stem}{strftime('%Y%m%d%H%M%S')}.log')
    
    ratios = Ratios()
    ratios.build(args.imax)
    p_commensurable = ratios.get_p()
    
    simple_ratios = SimpleRatios()
    names,periods = get_periods((Path(args.data)/args.planets).with_suffix('.csv')) 
    resonances = Resonance.create(names,periods)
    
    all_acceptable_resonances = [resonance for resonance in resonances if resonance.can_match(simple_ratios)]
    
    Label_Observations = ['All Planets']
    P_Observations = [ratios.get_P(len(all_acceptable_resonances),len(periods),p_commensurable)]
    Logger.instance.info(f'P ({Label_Observations[-1]}) = {P_Observations[-1]}')
    
    for primary in args.names:
        names,periods = get_periods((Path(args.data)/primary).with_suffix('.csv'),
                                    key='Satellite',
                                    e_max=args.e_max,
                                    R_min=args.R_min)
        
        resonances = Resonance.create(names,periods)
        acceptable_resonances = [resonance for resonance in resonances if resonance.can_match(simple_ratios)]
        if len(acceptable_resonances) > 0:
            Label_Observations.append(primary)
            P_Observations.append(ratios.get_P(len(acceptable_resonances),len(periods),p_commensurable))
            Logger.instance.info(f'P ({Label_Observations[-1]}) = {P_Observations[-1]}')
            all_acceptable_resonances += acceptable_resonances
        else:
            Logger.instance.info(f'No acceptable resonance for {primary}')
        
    unique_ratios = {(resonance.i1,resonance.i2) for resonance in all_acceptable_resonances}
    Logger.instance.debug('Unique resonances')
    for ratio in unique_ratios:
        Logger.instance.debug(f'{ratio[1]}:{ratio[0]}')

    colours = create_xkcd_colours()  

    fig = figure(figsize=(12,12))
    fig.suptitle(Path(__file__).stem)
    ax = fig.add_subplot(2,1,1)
    Logger.instance.debug('All resonances')
    counts = {}
    labels = {}
    for r in all_acceptable_resonances:
        try:
            counts[f'{r.i1}:{r.i2}'] += 1
            labels[f'{r.i1}:{r.i2}'] += f', {r.name1}-{r.name2}'
        except KeyError:
            counts[f'{r.i1}:{r.i2}'] = 1
            labels[f'{r.i1}:{r.i2}'] = f'{r.i1}:{r.i2} {r.name1}-{r.name2}'
 
        Logger.instance.debug(f'{r.name1}-{r.name2} {r.i2}:{r.i1}')
    keys =  sorted(list(counts.keys()),key=frac)  
    ax.bar(keys,[counts[key] for key in keys],
           width=0.5,
           color=colours,label=[labels[key] for key in keys])
    ax.set_xlim((-0.6, len(keys) - 0.4))
    ax.legend(ncols=2,title=f'There are {len(all_acceptable_resonances)} acceptable resonances')
    
    ax.set_xlabel('Resonances')
    ax.set_ylabel('Count')
    ax.set_xmargin(0.5)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    names = ', '.join(name.title() for name in args.names)
    ax.set_title(f'Planets, plus satellites of {names}: $R>${args.R_min}, $e<${args.e_max}, $imax=${args.imax}.')
 
    ax2 = fig.add_subplot(2,1,2)
    ax2.bar(Label_Observations,P_Observations,width=0.5,color=colours[len(keys):],label=Label_Observations)
    ax2.legend()
    ax2.set_xlabel('Primary')
    ax2.set_ylabel('P')
    ax2.set_title('Probability of observed number of commensurabilities')
    
    fig.savefig(Path(args.figs)/Path(__file__).stem)    
    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
    if args.show:
        show()
    
if __name__ == '__main__':
    main()
