import numpy as np, matplotlib.pyplot as plt, math as m 


@np.vectorize
def jacobi(x,y,n=1,mu2=0.2,Cj=3.9):
    n2  = n*n
    mu1 = 1-mu2
    x2  = x*x
    y2  = y*y
    r1  = m.sqrt((x+mu2)*(x+mu2) + y2)
    r2  = m.sqrt((x-mu1)*(x-mu1) + y2)
    diff= n2*(x2 + y2) + 2*(mu1/r1 + mu2/r2) - Cj
 
    return diff

def bounds(Z):
    minZ = float('inf')
    maxZ = - minZ
    for zz in Z:
        for z in zz:
            if z < minZ:
                minZ = z
            if z > maxZ:
                maxZ = z
    return (minZ,maxZ)     

plt.figure()
limit=5
xlist = np.linspace(-limit, limit, 100)   
ylist = np.linspace(-limit, limit, 100)
X,Y = np.meshgrid(xlist, ylist)   
Z = jacobi(X,Y)
(z0,z1) = bounds(Z)
print (z0,z1)
#levels=np.linspace(z0,z1,num=21,endpoint=True)
levels=list(range(-2,0,10))+list(range(0,60,10))
c = plt.pcolor(X,Y,Z)
d = plt.colorbar(c,orientation='horizontal')
#q = plt.winter()

CS3=plt.contourf(X, Y, Z, levels, colors=[ 'g'] )
plt.xlabel('x')
plt.ylabel('y')
plt.clabel(CS3)
CS3.cmap.set_under('m')
CS3.cmap.set_over('c')
#plt.colorbar(CS3)
plt.show()