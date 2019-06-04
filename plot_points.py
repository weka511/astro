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

if __name__=='__main__':
    import tracking, matplotlib.pyplot as plt
    log             = tracking.Tracking('run')
    x1s             = []
    y1s             = []
    x2s             = []
    y2s             = []    
    x3s             = []
    y3s             = []
    sampling_period = 1
    for name in log.get_reader():
        print (name)
        xs = log.get_one_status(name)
        i   = 0
        while 12*i<len(xs):
            if i%sampling_period==0:
                x1s.append(xs[12*i])
                y1s.append(xs[12*i+1])
                x2s.append(xs[12*i+2])
                y2s.append(xs[12*i+3])        
                x3s.append(xs[12*i+4])
                y3s.append(xs[12*i+5])
            i+=1
    plt.plot( x3s,y3s,'ro', x1s,y1s,'bo', x2s,y2s,'go')
    plt.show()