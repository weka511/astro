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

import unittest,solar,planet

class TestMarsMethods(unittest.TestCase):
    def setUp(self):
        self.mars = planet.Mars()
        self.solar = solar.Solar(self.mars)          

    def test_beam_irradience(self):
        self.assertAlmostEqual(1371/(1.5236915**2),
                               self.solar.beam_irradience(self.mars.a),
                               delta=0)
        
    def test_top_atmosphere(self):
        self.assertAlmostEqual(488,
                               self.solar.surface_irradience(69,22.3,14),
                               places=1)
    def test_top_atmosphere2(self):
        self.assertAlmostEqual(460,
                               self.solar.surface_irradience(69,22.3,13),
                               places=1) 
    def test_top_atmosphere7(self):
            self.assertAlmostEqual(29,
                                   self.solar.surface_irradience(69,22.3,19),
                                   places=1)    
        
try:
    unittest.main()
except SystemExit as inst:
    if inst.args[0] is True: # raised by sys.exit(True) when tests failed
        raise