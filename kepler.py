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

import math, utilities

def predict(h,f,x):
    #print x
    #print f(x)
    return [x0 + h * f0 for (x0,f0) in zip(x,f(x))]

def eta(x,m,k):
    return [
        -k/x[0],
        x[2]*x[2]/(2*m) + x[3]*x[3]/(2*m*x[0]*x[0]),
        x[3]
    ]

def correct(h,f,eta0,eta1,x0,x1):
    return [e0 + (h/2) * (f0 + f1) for (e0,f0,f1) in zip(eta0,f(eta0,x0),f(eta1,x1))]


    
def inverse(eta,m,k,x):
    r=-k/eta[0]
    L=eta[2]
    psq=2*m*eta[1]-L*L/(r*r)
    if psq>0:
        p=utilities.signum(x[2])*math.sqrt(2*m*eta[1]-L*L/(r*r))
    else:
        p=0
    return [r, x[1], p, L]

def kepler_pred(x,m,k): # r, theta, p, l
    return [
        x[2]/m, # p/m
        x[3]/(m*x[0]*x[0]), # l/(m*r*r)
        x[3]*x[3]/(m*x[0]*x[0]*x[0]) - k / (x[0]*x[0]),
        0
    ]

def kepler_corr(eta,m,k,x):
    term=k*x[2]/(m*x[0]*x[0])
    return [term,-term,0]

def hamiltonian(x,m,k):
    return x[2]*x[2]/(2*m) + x[3]*x[3]/(2*m*x[0]*x[0])-k/x[0]

def ff(eta,x):
    return kepler_corr(eta,m,k,x)

if __name__=='__main__':
    import matplotlib.pyplot as plt
    u=[]
    v=[]
    w=[]
    z=[]
    h=0.001
    k = 1
    r = 1.0000
    p = 0.0001
    m = 0.001
    L = math.sqrt((m/r*r*r ))
    x=[r,0,p,L]
    #print x
    nn=10000
    step = 1
    print hamiltonian(x,m,k)
     
    for i in range(nn):
        x1 = predict(h,lambda(z) : kepler_pred(z,m,k),x)
        eta0=eta(x,m,k)
        eta1=eta(x1,m,k)
        eta2=correct(h,ff,eta0,eta1,x,x1)
        x=inverse(eta2,m,k,x1)
        #print x
        if i%step==0:
            w.append(x1[0]*math.cos(x1[1]))
            z.append(x1[0]*math.sin(x1[1]))            
            u.append(x[0]*math.cos(x[1]))
            v.append(x[0]*math.sin(x[1]))
    plt.plot(w,z,'b',u,v,'r')
    
    print hamiltonian(x,m,k)