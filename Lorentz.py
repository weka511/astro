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

import rki

class Lorentz:
    def __init__(self,sigma,rho,beta):
        self.sigma=sigma
        self.rho=rho
        self.beta=beta
    def dx(self,x):
        return [
            self.sigma*(x[1]-x[0]),      \
            x[0]*(self.rho-x[2])-x[1],   \
            x[0]*x[1]-self.beta*x[2]
        ]
    
if __name__=='__main__':
    import matplotlib.pyplot as plt
    nn=1000
    h=0.01
    lorentz=Lorentz(10,28,8.0/3.0)
    rk=rki.ImplicitRungeKutta2(lambda (x): lorentz.dx(x),10,0.000000001)
    y=[1,0,1]
    u=[]
    v=[]
    for i in range(nn):
        y= rk.step(h,y)
        u.append(y[0])
        v.append(y[1])
    plt.plot(u,v) 