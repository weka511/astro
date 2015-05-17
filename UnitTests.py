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

# Tests

import unittest,solar,planet,utilities

class TestMarsMethods(unittest.TestCase):
    def setUp(self):
        self.mars = planet.Mars()
        self.solar = solar.Solar(self.mars)          

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
        
    def test_top_atmosphere(self):
        self.assertAlmostEqual(488,
                               self.solar.surface_irradience(69,22.3,14),
                               places=1)
    #def test_top_atmosphere2(self):
        #self.assertAlmostEqual(460,
                               #self.solar.surface_irradience(69,22.3,13),
                               #places=1) 
    #def test_top_atmosphere7(self):
            #self.assertAlmostEqual(29,
                                   #self.solar.surface_irradience(69,22.3,19),
                                   #places=1)    
        
try:
    unittest.main()
except SystemExit as inst:
    pass
    #if inst.args[0] is True: # raised by sys.exit(True) when tests failed
        #raise