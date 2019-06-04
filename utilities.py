# Copyright (C) 2015-2019 Greenweaves Software Limited

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

#Some useful functions that don't fit anywhere else

import string, time, math,__main__ as main,os

# get_data_file_name
#
# Construct data file name from name of main program

def get_data_file_name(path='data',ext='dat'):
    base = os.path.splitext(os.path.basename(main.__file__))[0]
    return os.path.join(path,'{0}.{1}'.format(base,ext))

# get_planetary_data
#
# Construct map containing elements for planets

def get_planetary_data(data_file_name):
    with open(data_file_name) as data_file:
        data = {}
        for line in data_file:
            parts = line.strip().split(',')  
            data[parts[0]] = [abs(float(part)) for part in parts[1:]]
        return data 
    


def signum(x):
    '''Determine sign of its argument. Returne -1, 0, or +1.'''
    if x<0: return -1
    if x>0: return +1
    return 0

def guarded_sqrt(x):
    '''Calculate a square root of a positive number, otherwise return 0'''
    return math.sqrt(x) if x>0 else 0

def newton_raphson(x,f,df,epsilon,N=50):
    '''
    Solve an equation using the Newton-Raphson method.
    
    ParametersL
       x       Starting value
       f       Function for equation: f(x)=0
       df      Derivative of f
       epsilon Maximum acceptable error
       N       Maximum number of iterations
    '''
    x0 = x
    for i in range(N):
        x1 = x0-f(x0)/df(x0)
        if abs(x1-x0)<epsilon:
            return x1
        else:
            x0 = x1    
    return x0

def get_angle(r):
    abs_theta=0 if r[0]==0 else math.atan(r[1]/r[0])
    return abs_theta+adjust_quadrant(r)

def adjust_quadrant(r):
    if r[0]>=0 and r[1]>=0: return 0
    if r[0]<0 and r[1]>=0: return math.pi/2
    if r[0]<0 and r[1]<0: return path.pi
    return 3*math.pi/2

def get_r(z):
    [x,y]=z
    return math.sqrt(x*x*y*y)

def get_r_velocity(zdot,theta):
    [xdot,ydot]=zdot
    return math.cos(theta)*xdot + math.sin(theta)*ydot

def get_theta_dot(zdot,theta,r):
    [xdot,ydot]=zdot
    return (math.cos(theta)*ydot - math.sin(theta)*xdot)/r
