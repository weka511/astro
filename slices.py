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

import math, planet, solar

class Slice:
    def __init__(self,latitude,longitude,depth):
        self.latitude = latitude
        self.longitude = longitude
        self.depth = depth

if __name__=="__main__":
    import matplotlib.pyplot as plt
    
    mars = planet.Mars()
    solar = Solar(mars)
        
    beam_irradience_top=solar.beam_irradience(mars.a)
    print "Mean beam irradience at top of atmosphere = {0:6.2f}".format(beam_irradience_top)    