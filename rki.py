# Copyright (C) 2015 Greenweaves Software Pty Ltd

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

import math

# TODO: more general case

class ImplicitRungeKutta(object):
    def step(self,h,y):
        raise NotImplementedError

class ImplicitRungeKuttaException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
class ImplicitRungeKutta2(ImplicitRungeKutta):
    def __init__(self,dy,max_iterations,iteration_error):
        self.gamma1          = 1.0/2.0
        self.gamma2          = 1.0/2.0
        self.alpha1          = 1.0/2.0-math.sqrt(3.0)/6.0
        self.alpha2          = 1.0/2.0+math.sqrt(3.0)/6.0
        self.beta11          = 1.0/4.0
        self.beta22          = 1.0/4.0
        self.beta12          = 1.0/4.0 - math.sqrt(3.0)/6.0
        self.beta21          = 1.0/4.0 + math.sqrt(3.0)/6.0
        self.dy              = dy
        self.max_iterations  = max_iterations
        self.iteration_error = iteration_error

    def step(self,h,y):
        (k1,k2)=self.iterate(h,y,[0 for yy in y],[0 for yy in y])
        for i in range(self.max_iterations):
            (k1_new,k2_new)=self.iterate(h,y,k1,k2)
            if self.distance(k1,k1_new) < self.iteration_error and self.distance(k2,k2_new) < self.iteration_error:
                return [yy + self.gamma1*kk1 + self.gamma2*kk2 for (yy,kk1,kk2) in zip(y,k1_new,k2_new)]
            else:
                k1,k2 = k1_new,k2_new
        raise ImplicitRungeKuttaException(                     \
            'Failed to Converge within {0} after {1} iterations'.format(  \
                self.iteration_error,                          \
                self.max_iterations)                           \
        )
 
    def iterate(self,h,y,k1,k2):
        y_k1_k2 = zip(y,k1,k2)
        return (
            [h * ff for ff in self.dy([yy+self.beta11*kk1+self.beta12*kk2 for (yy,kk1,kk2) in y_k1_k2])],
            [h * ff for ff in self.dy([yy+self.beta21*kk1+self.beta22*kk2 for (yy,kk1,kk2) in y_k1_k2])]
        )
    
    def distance(self,k,k_new):
        return max([abs(a-b) for (a,b) in zip(k,k_new)])
    

    
if __name__=='__main__':
    import matplotlib.pyplot as plt
    nn=1000
    h=2*math.pi/nn
    rk=ImplicitRungeKutta2(lambda (y): [y[1],-y[0]],10,0.000000001)
    y=[1,0,1]
    u=[]
    v=[]
    for i in range(nn):
        y= rk.step(h,y)
        u.append(y[0])
        v.append(y[1])
    plt.plot(u,v) 
    

        