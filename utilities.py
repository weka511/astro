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

import string, time, math

# slip_zip
#
# given a list, return a new list of triples, with the list "slipped", e.g.
# [1,2,3,4,5] becomes
# [(None, 1, 2), (1, 2, 3), (2, 3, 4), (3, 4, 5), (4, 5, None)]
# This function is uses to iterate through the layers, with each layer being
# sanwiched between the layes immediately above and below.

def slip_zip(x):
    return zip ([None]+x[:-1],x,x[1:]+[None])

# format_latitude
#
# Convert latitude to North/South form
def format_latitude(latitude):
    if latitude>0:
        NS="N"
    elif latitude<0:
        NS="S"
    else: 
        NS=' '
    return (NS,abs(latitude))

colours=[
    'b',   #Blue
    'r',   #Red
    'g',   #Green
    'c',   #Cyan
    'm',   #Magents
    'y',   #Yellow
    'k'    #Black
]

def get_colour(index):
    return colours[index%len(colours)]

# Find extremem, given that y0<y1>y2, or y0>y1<y2
# Fit a parabola to (x0,y0, (x1,y1), and (x2,y2), and find its extremum

def extremum(x0,x1,x2,y0,y1,y2):
    return (0.5 * 
            ((x2-x1)*(x2+x1)*y0 + (x0-x2)*(x0+x2)*y1 + (x1-x0)*(x1+x0)*y2) /
            ((x2-x1)*y0 + (x0-x2)*y1 + (x1-x0)*y2))

def signum(x):
    if x<0: return -1
    if x>0: return +1
    return 0

def guarded_sqrt(x):
    return math.sqrt(x) if x>0 else 0

def newton_raphson(x,f,df,epsilon,N):
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

# Used to record temperaturs in a log
class TemperatureRecord:
    def __init__(self,day,hour,hours_in_day):
        self.temperatures=[]
        self.day=day+hour/float(hours_in_day)
        
    def add(self,temperature):
        self.temperatures.append(temperature)

# Used to store temperature values (abstract class,
# needs to be implemeneted below
class TemperatureLog:
    def add(self,record):
        raise NotImplementedError("TemperatureLog.add(...)")
    def write(line):
        pass
    def get_max_min(self,channel):
        xs,ys=self.extract(channel)
        ys_for_period=[]
        xxx=[]
        ymin=[]
        ymax=[]
        x_previous=-1
        for x,y in zip(xs,ys):
            if x_previous<0:
                x_previous=x
                ys_for_period=[]
            if x-x_previous>=1:
                xxx.append(x_previous)
                ymin.append(min(ys_for_period))
                ymax.append(max(ys_for_period))
                x_previous=x
                x_previous=-1
                ys_for_period=[]
            ys_for_period.append(y)

        return (xxx,ymin,ymax)    
    def close(self):
        pass
        
# Used to store temperature values in an external file     
class ExternalTemperatureLog(TemperatureLog):
    def __init__(self,logfile,sep=' '):
        self.logfile=logfile
        self.sep=sep
        self.skipping=True
        self.start=time.time()
        
    def add(self,record):
        if self.skipping:
            self.logfile.write("START\n")
            self.skipping = False        
        self.logfile.write("{0:f}".format(record.day))
        for temperature in record.temperatures:
            self.logfile.write("{0} {1}".format(self.sep,temperature))
        self.logfile.write('\n')

    def write(self,line):
        self.logfile.write(line+'\n')
        
    def extract(self,channel):
        x=[]
        y=[]
        for line in self.logfile:
            if line[0:3]=="END": return (x,y)
            if self.skipping:
                self.skipping = line[0:5]!="START"
            else:
                parts=line.split()
                x.append(float(parts[0]))
                y.append(float(parts[channel]))
        self.rewind()
        self.skipping = True
        return (x,y)

    def get(self,key):
        for record in self.logfile:
            pair = string.split(record,'=')
            if len(pair)>0 and pair[0].lower()==key.lower():
                result=float(pair[1])
                self.rewind()
                return result
            
    #rewind, in case we want to extract data again
    def rewind(self):
        self.logfile.seek(0)
     
    def close(self):
        self.logfile.write("END, Elapsed time = {0:.1f} seconds\n".format(time.time()-self.start))
        
# Used to store temperature values internally           
class InternalTemperatureLog(TemperatureLog):
    def __init__(self):
        self.history=[]
        
    def add(self,record):
        self.history.append(record)
    
    def extract(self,layer_number,skip=1):
        days=[]
        result=[]
        for record in self.history:
            days.append(record.day)
            result.append(record.temperatures[layer_number])
        return (days,result)

if __name__=="__main__":
    with open('output.txt', 'w') as f:
        log=ExternalTemperatureLog(f)
        r1=TemperatureRecord(0,1,24)
        r1.add(1)
        r1.add(1.345)
        r1.add(99)
        log.add(r1)
        r2=TemperatureRecord(1,1,24)
        r2.add(2)
        r2.add(2.345678)
        r2.add(99)
        log.add(r2) 
        
        for a,b,c in slip_zip([1,2,3,4,5]):
            print (a,b,c)
        