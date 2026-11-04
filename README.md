# Astronomical Calculations

## Examples from [Murray and Dermott, Solar System Dynamics](https://www.cambridge.org/core/books/solar-system-dynamics/108745217E4A18190CBA340ED5E477A2)

#|Ex|File| Description 
----|-----|--------------------|------------------------------------------------------------
1||-| **Structure of the Solar System**
1.2|1.1|ex_1_1.py|Murray & Dermott, Exercise 1.1. Calculate major axes using Kepler's model (Platonic solids). Fit best possible permuation of Platonic Solids, neglecting eccentricty.
-|1.2|ex_1_2.py|Murray & Dermott, Exercise 1.2. Calculate perihelion and aphelion distances using Kepler's model (Platonic solids), allowing for  eccentricity
-||platonic.py|Common code needed by Murray & Dermott, Exercises 1 and 2
1.5|1.3|ex_1_3.py|Murray & Dermott, Exercise 1.3. Create a large number of random sets of orbital periods for a model satellite system, similar to the three inner satellites of Uranus. For each set, calculate the mean motion for each satellite and create hostogram for |n1-3n2+2n3}.
1.7|1.4|ex_1_4.py|Taking the orbital periods listed in Table A.9 but excluding Epimetheus, Telesto, Calypso and Helen, use the criteria given in Section 1.7 to show there are 28 ratios of mean motions to consider in the Saturn System.
-|1.5|ex_1_5.py|We are given the periods of six planets orbiting a star. By considering the fifteen possible ratios and the ten first order commensurabilities, identify the pairs of planets such that ratios are within 0.0001 of the commensurabilities. Estimate probability of this occurring by chance if the periods were randomly distributed. 
-|1.6|ex_1_6.py|Exercise 1.6 (e) Use the data in Appendix A to find the periods of all possible pairs of periods among the planets and the prograde satellites of Mars, Jupiter, Saturn, Uranus, and Neptune, with mean radii > 100 km and orbital eccentricities < 0.15. Taking i_max =7, show that thirty pairs of objects have ratios of orbital periods within epsilon-max of a permitted commensurability.
-||ratios.py|Functions and classes to support exercises from Section 1.7
2|||**The Two Body Problem**
2.4|2.2|ex_2_2.py|Determine average times or orbital conjunction between Earth & Mars. Show that the minimum distance varies by factor of almost 2. Determine orbital motions over the period 1982-2002. 
-||orbital.py|Orbital calculations to support Exercise 2.2
3|||**The Restricted Three Body Problem**
3.2||restricted.py|Integrate equations (3.16) & (3.17) for particle near Lagrange Point L4
3.3||jacobi.py|Zero velocity curves for the Jacobi Integral
-||jacobi3d.py|Potential surfaces for the Jacobi Integral
3.6||lagrange.py|Calculate Equilibrium points L1, L2, and L3, plus Jacobi constants, for Lagrange configuration of the 3 body problem
-||fig_3_7.py|Plot the zero-velocity curves for the three critical values of the Jacobi Constant that pass through the points L1, L2, and L3
-||fig_3_8.py|Plot the three-dimensional surface defined by Cj=2U and the Lagrangian equilibrium points.
3.7||fig_3_10.py|Figure 3-10 Plot the radial distances of the five Lagrangian points from the centre of mass as a function of log mu2
3.9||ex.sh|Batch file to generate figure 3.16(a), 3.16(b),3.17(a), and 3.17(b)
A|||Appendix A, **Solar System Data**
-||julian.py|Conversion between Calendar dates and Julian dates
-||md_data.py|Accessor function for data from Appendix A, plus code to validate data

## ODE integrators 

 File | Purpose |
--------------------------|---------------------------------------------------------------------
adams.py|Adams-Bashforth-Moulton Predictor/Corrector
euler.py|Simple Euler method for integrations ODE
fehlberg.py|Explore use of Fehlberg's method
integrators.py|Used to set up an integrator
rk4.py|Explicit Runge Kutta integrators
rki.py|Implicit Runge Kutta (symplectic) integrators
template.py|Template for python script
Lorentz.py|This tests an Integrator by calculating the evolution of the Lorentz Attractor
verlet.py|Integrator using Velocity Verlet algorithm
