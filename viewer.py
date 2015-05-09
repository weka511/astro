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

import utilities

def display(inputfile,figure):
    with open(inputfile, 'r') as f:
        history = utilities.ExternalTemperatureLog(f)
        x,y=history.extract(1)
        plt.figure(figure)
        plt.title("Diurnal variation in temperature for {0}".format(inputfile))
        plt.xlabel("Time")
        plt.ylabel("Temperature - Kelvin")
        plt.grid(True)    
        plt.plot(x,y)    
        
if __name__=="__main__":
    import matplotlib.pyplot as plt
    figure=0
    for name in ['equator.txt',
                 '10S.txt',
                 '20S.txt',
                 '30S.txt',
                 '40S.txt',
                 '50S.txt',
 #                '60S.txt',
                 '70S.txt',
                 '80S.txt',
                 '90S.txt',
                 ]:
        figure+=1
        display(name,figure)
    
    plt.show()
