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

def predict(h,f,x):
    #print x
    #print f(x)
    return [x0 + h * f0 for (x0,f0) in zip(x,f(x))]

def correct(h,f,x,xpred):
    return [x0 + (h/2) * (f0 + fpred) for (x0,f0,fpred) in zip(x,f(x),f(xpred))]

def kepler_pred(x,m,k): # r, theta, p, l
    return [
        x[2]/m, # p/m
        x[3]/(m*x[0]*x[0]), # l/(m*r*r)
        x[3]*x[3]/(m*x[0]*x[0]*x[0]) - k / (x[0]*x[0]),
        0
    ]

def hamiltonian(x,m,k):
    return x[2]*x[2]/(2*m) + x[3]*x[3]/(2*m*x[0]*x[0])-k/x[0]

if __name__=='__main__':
    import matplotlib.pyplot as plt
    u=[]
    v=[]
    w=[]
    z=[]
    h=0.01
    k = 1
    r = 1
    m = 0.001
    L = math.sqrt((m/r*r*r ))
    x=[r,0,0,L]
    #print x
    nn=10000
    step = 100
    for i in range(nn):
        x1 = predict(h,lambda(z) : kepler_pred(z,m,k),x)
        x=correct(h,lambda(z) : kepler_pred(z,m,k),x,x1)
        #print x
        if i%step==0:
            w.append(x1[0]*math.cos(x1[1]))
            z.append(x1[0]*math.sin(x1[1]))            
            u.append(x[0]*math.cos(x[1]))
            v.append(x[0]*math.sin(x[1]))
    plt.plot(w,z,'b',u,v,'r')
    
    print hamiltonian(x,m,k)