# Bioregeneration-IX-Modeling

This repo has code associated with the following journal article: 

"Aponte-Morales, V. E., Payne, K. A., Cunningham, J. A., & Ergas, S. J. (2018). Bioregeneration of 
chabazite during nitrification of centrate from anaerobically digested livestock waste: experimental 
and modeling studies. Environmental science & technology, 52(7), 4090-4098."

The Matlab script IX_model implements a homogeneous surface diffusion model (HSDM)
which predicts the kinetics of IX between NH4+ and Na+
At the surface of the chabazite (solid-water interface), assume
an "ion exchange isotherm" applies:
q = Q*K*[NH4+] / {[Na+] + K*[NH4+]}
Best-fit values of Q and K come from equilibrium isotherm
experiments performed by Veronica Aponte-Morales and analyzed by
Karl Payne.

The solution of the differential equations is a finite difference method

