# Astronomical Calculations inspired by Mike Brown's course [Science of the Solar System](https://www.coursera.org/learn/solar-system/home/info)

## Orbital calculations: [Kotovich & Bowman: an Exactly Conservative Integrator for the n-body Problem](http://arxiv.org/pdf/physics/0112084)

| File | Purpose |
|--------------------------|---------------------------------------------------------------------|
|kepler.py|Hamiltonian for integrating Kepler problem|
|integrators.py|Simple integrator based on Kotovich & Bowman|
|plot_points.py|Display data that has been stored by tracking.py|
|rki.py|Implicit Runge Ketta (symplectic) integrators|
|restricted.py|Restricted 3 body problem|
|Threebody.py|Hamiltonian for 2 Dimensional, but otherwise general,  3 body problem|
|tracking.py|Record results in logfile so they can be played back, and analyses can be restarted|
|Lorentz.py|This tests the ImplicitRungeKutta Integrator by calculating the evolution of the Lorentz Attractor|

## Examples from [Murray and Dermott, Solar System Dynamics](https://www.cambridge.org/core/books/solar-system-dynamics/108745217E4A18190CBA340ED5E477A2)

|#|File| Description |
|----|--------------------|------------------------------------------------------------|
|1|-| Structure of the Solar System|
|1.3|uranus.py|Satellites of Uranus|
|1.4|saturn.py|Satellites of Saturn|
|1.5|commensurability.py|Identify commensurability and estimate probability of this particular value occurring by chance.|
|2||The Two Body Problem|
||orbit_space.tex|Calculations to support values for *orbital.py*|
|2.2|earth_mars.py|Earth Mars Distance at opposition|
||orbital.py|Orbital calculations|
|3||The Restricted Three Body Problem|
||jacobi.py|Zero velocity Surfaces for the Jacobi Integral|
||jacobi3d.py|Potential surfaces|

## Miscellaneous

| File | Purpose |
|--------------------------|---------------------------------------------------------------------|
|astro.bib|Biblography|
