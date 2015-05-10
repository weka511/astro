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

import utilities, io, os,matplotlib.pyplot as plt

def display(inputfile,figure):
    with open(inputfile, 'r') as f:
        history = utilities.ExternalTemperatureLog(f)
        plt.figure(figure)
        plt.title("Diurnal variation in temperature for {0}".format(inputfile))
        plt.xlabel("Time")
        plt.ylabel("Temperature - Kelvin")
        plt.grid(True)
        x1,y1=history.extract(1)
        x2,y2=history.extract(2) 
        x3,y3=history.extract(3)
        x4,y4=history.extract(4)
        x5,y5=history.extract(5)
        x10,y10=history.extract(10)
        x20,y20=history.extract(20)
        plt.plot(x1,y1,'r-',x2,y2,'g-',x3,y3,'b-',x4,y4,'c-',x5,y5,'m-',x10,y10,'y-',x20,y20,'k-')
        plt.savefig(os.path.splitext(inputfile)[0])
        
if __name__=="__main__":

    figure=0
    for name in [
        'equator.txt',
        '10S.txt',
        '20S.txt',
        '30S.txt',
        '40S.txt',
        '50S.txt',
        '60S.txt',
        '70S.txt',
        '80S.txt',
        '90S.txt',
        '10N.txt',
        '30N.txt',
        '20N.txt',
        '40N.txt',
        '50N.txt',
        '60N.txt',
        '70N.txt',
        '80N.txt',
        '90N.txt',
        ]:
        figure+=1
        display(name,figure)
    
    plt.show()
