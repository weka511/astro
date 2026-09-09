# Astronomical Calculations inspired by Mike Brown's course [Science of the Solar System](https://www.coursera.org/learn/solar-system/home/info)

## Orbital calculations: [Kotovich & Bowman: an Exactly Conservative Integrator for the n-body Problem](http://arxiv.org/pdf/physics/0112084)

 File | Purpose |
--------------------------|---------------------------------------------------------------------
kepler.py|Hamiltonian for integrating Kepler problem
integrators.py|Simple integrator based on Kotovich & Bowman
plot_points.py|Display data that has been stored by tracking.py
rki.py|Implicit Runge Ketta (symplectic) integrators
restricted.py|Restricted 3 body problem
Threebody.py|Hamiltonian for 2 Dimensional, but otherwise general,  3 body problem
tracking.py|Record results in logfile so they can be played back, and analyses can be restarted
Lorentz.py|This tests the ImplicitRungeKutta Integrator by calculating the evolution of the Lorentz Attractor

## Examples from [Murray and Dermott, Solar System Dynamics](https://www.cambridge.org/core/books/solar-system-dynamics/108745217E4A18190CBA340ED5E477A2)

#|File| Description 
----|--------------------|------------------------------------------------------------
1|-| Structure of the Solar System
1.3|uranus.py|Murray & Dermott, Exercise 1.3. Create a number of random sets of orbital periods orbital periods for a model satellite system, similar to the three inner satellites of Uranus. For each set, calculate the mean motion for each satellite and check n1-3n2+2n3.
1.4|saturn.py|Murray & Dermott, Exercise 1.4. Taking the orbital periods listed in Table A.9 but excluding Epimetheus, Telesto, Calypso and Helen, use the criteria gSiven in Section 1.7 to show there are 28 ratios of mean motions to consider in the Saturn System.
1.5|commensurability.py|Murray & Dermott, Exercise 1.5. Identify commensurability and estimate probability of value occurring by chance.
2||The Two Body Problem
2.2|earth_mars.py|Murray and Dermott, Exercise 2.2: times of orbital conjunction between earth & Mars.
||orbital.py|Orbital calculations to support *earth_mars.py*
3||The Restricted Three Body Problem
||jacobi.py|Zero velocity Surfaces for the Jacobi Integral
||jacobi3d.py|Potential surfaces

## Docs directory

 File | Purpose 
--------------------------|---------------------------------------------------------------------
astro.bib|Bibliography
orbit_space.tex|Derivations of equations for *earth_mars.py*

## Data files

File | Description
---------------------------|---------------------------------------------------------------------
commensurability.csv|Orbital periods for exercise 1.4
earth_mars.csv|Data files used by *earth_mars.py*
saturn.csv|Orbital periods for exercise 1.5

