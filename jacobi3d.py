
from mpl_toolkits.mplot3d import Axes3D
import numpy as np, matplotlib.pyplot as plt, math as m, jacobi,matplotlib.colors as clrs

def plot_3d(limit=2,eps=0.001,minZ=-3.81,maxZ=-2.8):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
       
    xlist = np.linspace(-limit, limit+eps, 100)   
    ylist = np.linspace(-limit, limit+eps, 100)
    X,Y = np.meshgrid(xlist, ylist)   
    Z = -jacobi.jacobi(X,Y)
    Z[Z<minZ]= np.nan
    Z[Z>maxZ]= np.nan
    norm = clrs.Normalize(vmin = minZ, vmax = maxZ, clip = False)
    ax.plot_surface(X, Y, Z,cmap=plt.cm.jet,norm=norm)

if __name__=='__main__':    
    plot_3d()
