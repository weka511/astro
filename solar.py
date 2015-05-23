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

# Model for solar irradiation, based on Solar Radiation on Mars, 
# Joseph Appelbaum & Dennis Flood, Lewis Research Center, NASA 

import math, planet,utilities

class Solar:
    def __init__(self,planet):
        self.planet=planet
        
#   Beam Irradience in W/m2
#   Appelbaum & Flood equation (1)
    def beam_irradience(self,r):
        return self.planet.S/(r*r)
 
 # Beam irradience on a horizonal surface
    #   Appelbaum & Flood equations (5) & (6)
    def surface_irradience(self,areocentric_longitude,latitude,T):
        cos_zenith_angle = self.planet.cos_zenith_angle(areocentric_longitude,latitude,T)
        beam_irradience = self.beam_irradience(self.planet.instantaneous_distance(areocentric_longitude))
        return max(0,cos_zenith_angle*beam_irradience)
    
if __name__=="__main__":
    import matplotlib.pyplot as plt, unittest
    from scipy.integrate import quad
    
    class TestMarsMethods(unittest.TestCase):
        def setUp(self):
            self.mars = planet.Mars()
            self.solar = Solar(self.mars)          
    
    #   This test is based on figure 3 of Appelbaum & Flood, with two adjustments.
    #   Appelbaum & Flood have the perihelion and aphelion at LS = 249 degress and
    #   68 degress respectively (Figure 3). But in the section preceding equation (3)
    #   the perihelion is stated to be at 248 degrees. The aphelion is 180 degrees
    #   before or after the perihelion, so this has also been shited.
    
        def test_beam_irradience_as_function_areocentric_longitude(self):
            d0=-1
            d1=-1
            for i in range(360):
                d2=self.mars.instantaneous_distance(i)
                if d0>-1 and d1>-1:
                    if d0>d1 and d1<d2:
                        x=utilities.extremum(i-2,i-2,i,d0,d1,d2)
                        d=self.mars.instantaneous_distance(x)
                        irr=self.solar.beam_irradience(d)
                        self.assertAlmostEqual(248,x)
                        self.assertAlmostEqual(718,irr,delta=1)
                    if d0<d1 and d1>d2:
                        x=utilities.extremum(i-2,i-2,i,d0,d1,d2)
                        d=self.mars.instantaneous_distance(x) 
                        irr=self.solar.beam_irradience(d)
                        self.assertAlmostEqual(68,x)
                        self.assertAlmostEqual(493,irr,delta=1)
                d0=d1
                d1=d2   
                 
        def test_beam_irradience(self):
            self.assertAlmostEqual(1371/(1.5236915**2),
                                   self.solar.beam_irradience(self.mars.a),
                                   places=1)

        # The next few tests are based on Table II i Appelbaum & Flood
        
        def test_top_atmosphere(self):
            integral,error=quad(integrand,12,13,args=(self.solar,69))
            self.assertAlmostEqual(488,integral,delta=1.5)
            
        def test_top_atmosphere2(self):
            integral,error=quad(integrand,13,14,args=(self.solar,69))
            self.assertAlmostEqual(460,integral,delta=1)            
 
         
        def test_top_atmosphere3(self):
            integral,error=quad(integrand,14,15,args=(self.solar,249))
            self.assertAlmostEqual(376,integral,delta=1)             
 
        def test_top_atmosphere7(self):
            integral,error=quad(integrand,18,19,args=(self.solar,69))
            self.assertAlmostEqual(25,integral,delta=1)

 
    def integrand(x, solar,ls):
        return solar.surface_irradience(ls,22.3,x)
    
    try:
        unittest.main()
    except SystemExit as inst:
        pass    
    
    
    def generate_irradiance(planet,areocentric_longitude,latitude):
        x=[]
        y=[]
        for T in range(12,22):
            irradiance=solar.surface_irradience(areocentric_longitude,latitude,T)
            x.append(T)
            y.append(irradiance)
        return (x,y)
     
    mars = planet.create("Mars")
    solar = Solar(mars)
    
    beam_irradience_top=solar.beam_irradience(mars.a)
    print "Mean beam irradience at top of atmosphere = {0:6.2f} W/m2".\
          format(beam_irradience_top)
    
    plt.figure(3)
    plt.title("Beam irradience at top of Mars atmosphere")
    plt.xlabel("Areocentric longitude - degrees")
    plt.ylabel("Beam irradience at top of Mars atmosphere - W/m2")
    plt.grid(True) 
    xs=[]
    ys=[]
    d0=-1
    d1=-1

    for i in range(360):
        xs.append(i)
        d2=mars.instantaneous_distance(i)
        ys.append(solar.beam_irradience(d2))
        if d0>-1 and d1>-1:
            if d0>d1 and d1<d2:
                x=utilities.extremum(i-2,i-2,i,d0,d1,d2)
                d=mars.instantaneous_distance(x)
                irr=solar.beam_irradience(d)
                print "perihelion", x,d,irr
            if d0<d1 and d1>d2:
                x=utilities.extremum(i-2,i-2,i,d0,d1,d2)
                d=mars.instantaneous_distance(x) 
                irr=solar.beam_irradience(d)
                print "aphelion", x,d,irr
        d0=d1
        d1=d2
    plt.plot(xs,ys)
    
    plt.figure(4)
    plt.title("Variation of solar declination angle")
    plt.axis([0, 360, -25, 25])
    plt.xlabel("Areocentric longitude - degrees")
    plt.ylabel("Solar Declination Angle - degrees")
    plt.grid(True) 
    x=[]
    y=[]
    for i in range(360):
        x.append(i)
        y.append(math.degrees(math.asin(mars.sin_declination(i))))
    plt.plot(x,y)

    plt.figure(6)
    plt.title("Diurnal Variation of Beam Irradience on a horizontal surface")
    plt.xlabel("Solar Time - Hours")
    plt.ylabel("Beam Irradiance - W/m2")
    (x1,y1)=generate_irradiance(mars,69,22.3)
    (x2,y2)=generate_irradiance(mars,120,22.3)
    (x2,y2)=generate_irradiance(mars,120,22.3)
    (x3,y3)=generate_irradiance(mars,153,22.3)
    (x4,y4)=generate_irradiance(mars,249,22.3)
    (x5,y5)=generate_irradiance(mars,299,22.3)
    plt.plot(x1,y1,"r",x2,y2,"g",x3,y3,"b",x4,y4,"c",x5,y5,"m")
    plt.axis([12, 20, 0, 600])
    plt.grid(True)    
    plt.show()    