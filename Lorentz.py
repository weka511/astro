# Copyright (C) 2015-2017 Greenweaves Software Pty Ltd

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

# This tests the ImplicitRungeKutts Inntegrator by calculating the
# evolution of the Lorentz Attractor



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
    def __str__(self):
        return 'sigma={0},rho={1}.beta={2}'.format(self.sigma,self.rho,self.beta)
       
if __name__=='__main__':
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import rki
    
    nn=10000
    h=0.01
    lorentz=Lorentz(10,28,8.0/3.0)
    rk=rki.ImplicitRungeKutta2(lambda x: lorentz.dx(x),10,0.0000001)
    v=[1,0,1]
    xs=[]
    ys=[]
    zs=[]
    for i in range(nn):
        v= rk.step(h,v)
        xs.append(v[0])
        ys.append(v[1])
        zs.append(v[2])
        
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')    
    ax.scatter(xs, ys, zs)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(lorentz)
    plt.show()    