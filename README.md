# astro
Astronomical Calculations inspired by Mike Brown's "Science of the Solar System" course

## Leighton and Murray

I have revisited [Leighton and Murray's Behavior of Carbon Dioside and Other Volatiles on Mars](http://www.mars.asu.edu/christensen/classdocs/Leighton_BehavioCO2_science_66.pdf)


| File | Purpose |
| ------------------------- | ------------------------------------------------------------|
| leighton.py | Driver for executing model |
| planet.py |  Repository for basic data about planets |
| solar.py |  Model for solar irradiation, based on [Solar Radiation on Mars, Joseph Appelbaum & Dennis Flood, Lewis Research Center, NASA](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19890018252.pdf) |
| thermalmodel.py | Slices of Mars' interior, together with model for heat flow |
| utilities.py | Utility functions for log files, zipping lists, choosing colours for plots |
| viewer.py | Used to plot data files from leighton.py |
| physics.py | Repository for physical laws and constants |
| kepler2.py | Responsible for determining distance of planet from Sun at a specified time |
## Orbital calculations

I have started working on some code for orbital calculatione, based on [Kotovich & Bowman: an Exactly Conservative Integrator for the n-body Problem](http://arxiv.org/pdf/physics/0112084)

| File | Purpose |
|--------------------------|---------------------------------------------------------------------|
| kepler.py |  Hamiltonian for integrating Kepler problem |
| integrators.py |  Simple integrator based on Kotovich & Bowman |
| plot_points.py | Display data that has been stored by tracking.py |
| rki.py | Implicit Runge Ketta (symplectic) integrators |
| restricted.py| Restricted 3 body problem |
| Threebody.py | Hamiltonian for 2 Dimensional, but otherwise general,  3 body problem |
| tracking.py | Record results in logfile so they can be played back, and analyses can be restarted |
| Lorentz.py | This tests the ImplicitRungeKutts Integrator by calculating the evolution of the Lorentz Attractor |

## Murray and Dermott

Code based on [Murray and Dermott, Solar System Dynamics](https://www.cambridge.org/core/books/solar-system-dynamics/108745217E4A18190CBA340ED5E477A2)

| File | Purpose |
|--------------------------|---------------------------------------------------------------------|
| jacobi.py |  Zero velocity Surfaces for the Jacobi Integral |
| jacobi3d.py |  Potential surfaces |

