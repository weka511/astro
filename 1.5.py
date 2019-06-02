# Copyright (C) 2019 Greenweaves Software Limited

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this software.  If not, see <http://www.gnu.org/licenses/>

def identify_commensurabilities(periods,tolerance=0.001):
     ratios = sorted([(i,j,periods[i]/periods[j]) for j in range(len(periods)) for i in range(j)],
                     key=lambda x:x[2])

     for p in range(1,11):
          commensurability = p/(p+1)
          for i,j,ratio in ratios:
               if abs(ratio-commensurability)<tolerance:
                    print ('i={0},j={1},ratio={2},commensurability={3}:{4}({5}) [{6}]'.format(i,
                                                                                        j,
                                                                                        ratio,
                                                                                        p,
                                                                                        p+1,
                                                                                        commensurability,
                                                                                        abs(commensurability-ratio)))
    
if __name__=='__main__':
     with open('data/1.5.txt') as data:
          periods = []
          for line in data:
               periods.append(float(line.strip()))
          identify_commensurabilities(periods)
            