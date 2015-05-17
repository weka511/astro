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

# given a list, return a new list of triples, with the list "slipped", e.g.
# [1,2,3,4,5] becomes
# [(None, 1, 2), (1, 2, 3), (2, 3, 4), (3, 4, 5), (4, 5, None)]
# This function is uses to iterate through the layers, with each layer being
# sanwiched between the layes immediately above and below.


def slip_zip(x):
    return zip ([None]+x[:-1],x,x[1:]+[None])

# Find extremem, given that y0<y1>y2, or y0>y1<y2
# Fit a parabola to (x0,y0, (x1,y1), and (x2,y2), and find its extremum

def extremum(x0,x1,x2,y0,y1,y2):
    return (0.5 * 
            ((x2-x1)*(x2+x1)*y0 + (x0-x2)*(x0+x2)*y1 + (x1-x0)*(x1+x0)*y2) /
            ((x2-x1)*y0 + (x0-x2)*y1 + (x1-x0)*y2))

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
    
# Used to store temperature values in an external file     
class ExternalTemperatureLog(TemperatureLog):
    def __init__(self,logfile,sep=' '):
        self.logfile=logfile
        self.sep=sep
        self.skipping=True
        
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
            if self.skipping:
                self.skipping = line[0:5]!="START"
            else:
                parts=line.split()
                x.append(float(parts[0]))
                y.append(float(parts[channel]))
        self.logfile.seek(0)  #rewind, in case we want to extract data again
        self.skipping = True
        return (x,y)

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
        
        print slip_zip([1,2,3,4,5])
        