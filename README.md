# astro
Astronomical Calculations inspired by Mike Brown's course [Science of the Solar System](https://www.coursera.org/learn/solar-system/home/info)

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

## Barnes Hut

The [galaxy](https://github.com/weka511/astro/tree/master/galaxy) folder contains a [Barnes Hut](https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation) simulation of a galaxy
