# galaxy

## Barnes Hut

[Barnes Hut](https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation) simulation of a galaxy

| File | Header | Purpose |
|-----------|--------------|---------------------------------------------------------------------|
| |catch.hpp | [Catch2](https://github.com/catchorg/Catch2) Unit testing framework |
| galaxy.cpp |galaxy.h| Simulate Galaxy using Barnes Hut|
| Makefile || Build Barnes Hut simulation |
| make_3d.py | |Plot output from galaxy |
| make_img.py || Plot output from galaxy |
|| particles.h | Represents the particles whose motion is being simulated|
| plot_orbits.py || Select a few stars at random and plot their orbits |
| treecode.cpp | treecode.h | Barbes-Hur Tree|
| utils.cpp | utils.h | Assorted utility functions|
| verlet.cpp | verlet.h | Verlet integrator|
