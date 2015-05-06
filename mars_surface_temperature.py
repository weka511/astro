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

# TODO:

#The daily and yearly position of the sun

#Heat flow into and out of the subsurface

#The effect of condensation and evaporation on the energy balance

#The total amount of volatile at the surface at each point in time

# see http://ccar.colorado.edu/asen5050/projects/projects_2001/benoit/solar_irradiance_on_mars.htm

# DF]Solar Radiation on Mars - NASA Technical Reports Server ...
#ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19890018252.pdf


import math

S = 1371 # Solar constant at the mean Sun-Earth distance of l AU, in N/m2
a = 1.5236915 # the Mars semimajor axis in AU,
e = 0.093377 # Mars eccentricity

def beam_irradience(r):
    return S/(r*r)

def instantanaeous_distance(longitude):
    theta = math.radians(longitude - 248)
    return a*(1-e*e)/(1 + e *math.cos(theta))

def beam_irradience_horizontal(G):
    coz_z=math.sin(latitude)*math.sin(declination) +                  \
        math.cos(latitude)*math.cos(declination)*math.cos(hour_angle)
    return G*cos_z


if __name__=="__main__":
    print "Mean beam irradience at top of atmosphere={0:6.2f}".format(beam_irradience(a))
    