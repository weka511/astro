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

def slip_zip(x):
    return zip ([None]+x[:-1],x,x[1:]+[None])


class TemperatureRecord:
    def __init__(self,day,hour,hours_in_day):
        self.temperatures=[]
        self.day=day+hour/float(hours_in_day)
        
    def add(self,temperature):
        self.temperatures.append(temperature)

class TemperatureLog:
     def add(self,record):
         raise NotImplementedError("TemperatureLog.add(...)")
     
class ExternalTemperatureLog(TemperatureLog):
    def __init__(self,logfile,sep=' '):
        self.logfile=logfile
        self.sep=sep
        
    def add(self,record):
        self.logfile.write("{0:f}".format(record.day))
        for temperature in record.temperatures:
            self.logfile.write("{0} {1}".format(self.sep,temperature))
        self.logfile.write('\n')

    def extract(self,channel):
        x=[]
        y=[]
        for line in self.logfile:
            parts=line.split()
            x.append(float(parts[0]))
            y.append(float(parts[channel]))
        self.logfile.seek(0)  #rewind, in case we want to extract data again
        return (x,y)
        
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
        