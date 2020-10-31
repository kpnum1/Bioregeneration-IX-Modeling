# Bioregeneration-IX-Modeling

This repo has code associated with the following journal article: 

"Aponte-Morales, V. E., Payne, K. A., Cunningham, J. A., & Ergas, S. J. (2018). Bioregeneration of 
chabazite during nitrification of centrate from anaerobically digested livestock waste: experimental 
and modeling studies. Environmental science & technology, 52(7), 4090-4098."


**Bioregen_model.m** 

Implements a homogeneous surface diffusion model (HSDM)
which predicts the kinetics of IX between NH4+ and Na+
The biological process is a two-step nitrification process
described by an Andrew's model to account for inhibition 

The solution of the differential equations is a finite difference method
& the plots are for a bioreactor amended with chabazite
Figure 1: NH4+ and Na+  concentrations with time
Figure 2: NO2- and NO3- concentrations with time

**Isotherm_and_Removal.m**

Plots the IX isotherm & plots
and calculates the removal of NH4+ as a 
function of chabazite dose


