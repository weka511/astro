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

import random, matplotlib.pyplot as plt, numpy as np

def generate_orbit(T1=1.413, L=1.546, U=2.101, range_i=range(1,3)):
    Ts = [T1]
    for i in range_i:
        x     = random.random()
        ratio = L + x*(U-L)
        Ts.append(ratio * Ts[-1])
    return Ts

if __name__=='__main__':
    N      = 1000000
    random.seed()
    delta_n = []
    for orbits in  [generate_orbit() for i in range(N)]:
        n = [360/T for T in orbits]
        delta_n.append(abs(n[0]-3*n[1] + 2*n[2]))
    bins=10**(np.arange(-2,4,dtype=float))
    plt.xscale('log')
    plt.hist(delta_n,bins=bins)
    plt.title('Probability: {0}%'.format(100*len([delta for delta in delta_n if delta<0.1])/100000))
    plt.show()
            