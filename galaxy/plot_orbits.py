'''
 Copyright (C) 2017 Greenweaves Software Pty Ltd

 This is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This software is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this software.  If not, see <http://www.gnu.org/licenses/>

  Companion to barnes_hut.cpp - plot images
'''
import os, re, sys, numpy as np, matplotlib.pyplot as plt,mpl_toolkits.mplot3d,random

colours=['r','g','b','m','c','y']
def plot(data,rows=[]):
    def get_coordinates(body,i):
        return [d[body][i] for d in data]
    ax = plt.gcf().add_subplot(111, aspect='equal', projection='3d')
    for body in range(len(data[0])):
        ax.scatter(get_coordinates(body,0)[0],get_coordinates(body,1)[0],get_coordinates(body,2)[0],c=colours[body%len(colours)],
                   label="Body: {0}".format(rows[body]))
        ax.scatter(get_coordinates(body,0),get_coordinates(body,1),get_coordinates(body,2),edgecolor=colours[body%len(colours)],s=1) 
        plt.legend(loc='best')


def extract(config_path = './configs/',rows=[0,1,2,55,100,400]):
    result=[]
    for file_name in os.listdir(config_path):
        m = re.search('body_[0-9]+00.dat',file_name)
        if m:
            positions = np.loadtxt(os.path.join(config_path,m.group(0)))
            result.append([positions[i] for i in rows])
    return result

if __name__=='__main__':
    M=6
    N=1000
    rows=random.sample(range(N),M) 
    data=extract(rows=rows)
    plot(data,rows=rows)
    plt.show()
    
        
