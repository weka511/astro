# astro
Astronomical Calculations inspired by Mike Brown's "Science of the Solar System" course

## Leighton and Murray

I have revisited [Leighton and Murray's Behavior of Carbon Dioside and Other Volatiles on Mars](http://www.mars.asu.edu/christensen/classdocs/Leighton_BehavioCO2_science_66.pdf)


| File | Purpose |
| ------------------------- | ----------------------------|
| leighton.py | Driver for executing model |
| planet.py |  Repository for basic data about planets |
| solar.py |  Model for solar irradiation, based on [Solar Radiation on Mars, Joseph Appelbaum & Dennis Flood, Lewis Research Center, NASA](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19890018252.pdf) |
| thermalmodel.py | Slices of Mars' interior, together with model for heat flow |
| utilities.py | Utility fucntions for log files, zipping lists, choosing colours for plots |
| viewer.py | Used to plot data files from leighton.py |

### TODO

 * The effect of condensation and evaporation on the energy balance
 * The total amount of volatile at the surface at each point in time

