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

import math, utilities, integrators

class Kepler(integrators.Hamiltonian):
    def __init__(self,x,m,k):
        self.x = x
        self.m = m
        self.k = k
        self.transform()
        
    def create(self,x):
        return Kepler(x,self.m,self.k)
    
    def transform(self):
        self.eta=[
            -self.k/self.x[0],
            self.x[2]*self.x[2]/(2*self.m) + \
            self.x[3]*self.x[3]/(2*self.m*self.x[0]*self.x[0]),
            x[3]
        ]
    
    def invert(self,kepler):
        r=-self.k/self.eta[0]
        L=self.eta[2]
        psq=2*self.m*self.eta[1]-L*L/(r*r)        
        p=utilities.signum(x[2])*math.sqrt(psq) if psq>0 else 0
        self.x[0] = r
        self.x[1] = kepler.x[1]
        self.x[2] = p
        self.x[3] = L

    def dx(self):
        return [
            self.x[2]/self.m, # p/m
            self.x[3]/(self.m*self.x[0]*self.x[0]), 
            self.x[3]*self.x[3]/(self.m*self.x[0]*self.x[0]*self.x[0]) - 
                self.k / (self.x[0]*self.x[0]),
            0
        ]    

    def d_eta(self):
        term=self.k*self.x[2]/(self.m*self.x[0]*self.x[0])
        return [term,-term,0]

    def hamiltonian(self):
        return self.x[2]*self.x[2]/(2*self.m) + \
               self.x[3]*self.x[3]/(2*self.m*self.x[0]*self.x[0])-self.k/self.x[0]

    def display(self):
        print "x",self.x
        print "eta",self.eta
        




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
    kepler=Kepler(x,m,k)
    nn=100000
    step = 1
    print kepler.hamiltonian()
    integrator = integrators.Integrate2(h,kepler)
    
    for i in range(nn):
        integrator.integrate()
        u.append(kepler.x[0]*math.cos(kepler.x[1]))
        v.append(kepler.x[0]*math.sin(kepler.x[1]))
    plt.plot(u,v)
  
    
    print kepler.hamiltonian()