import numpy as np, matplotlib.pyplot as plt, math as m 

origin = 'lower'

@np.vectorize
def jacobi(x,y,n=1,mu2=0.2,Cj=3.9):
    n2  = n*n
    mu1 = 1-mu2
    x2 = x*x
    y2 = y*y
    r1  = m.sqrt((x+mu2)*(x+mu2)+y2)
    r2  = m.sqrt((x-mu1)*(x-mu1)+y2)
    diff= n2*(x2 + y2) + 2*(mu1/r1 + mu2/r2)-Cj
    return clip(diff)

def clip(diff,limit=1):
    if diff>limit:
        return limit
    if diff<-limit:
        return -limit
    return diff

plt.figure()
xlist = np.linspace(-20.0, 20.0, 100)   
ylist = np.linspace(-20.0, 20.0, 100)
X,Y = np.meshgrid(xlist, ylist)   
Z = jacobi(X,Y)

levels=np.linspace(-0.5,0.5,num=21,endpoint=True)
CS3=plt.contour(X, Y, Z, 
            levels,
            colors=[ 'r','g','b'],
            #linestyles = 'solid',
            #origin=origin,
            #extend='both'
            )
plt.clabel(CS3)
CS3.cmap.set_under('m')
CS3.cmap.set_over('c')
#plt.colorbar(CS3)
plt.show()