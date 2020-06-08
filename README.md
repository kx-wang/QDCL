# QDCL

This code calculates the wavefunctions and simulates for the transport properties in THz 
Quantum Dot Cascade Lasers (QDCLs). The effects of longitudinal-optical (LO) phonons and longitudinal-acoustic (LA) phonons on the transport are included. The transport model used is a time-local density matrix model, within the second-order Born-Markov approximation, also known as the Interaction approach or Redfield equations.

# Instructions: 
1. Calculate the wavefunctions. 

   -run script_one_module.m, script_two_modules.m, to calculate the wavefunctions in the z (growth) direction 
   
   -run script_Radial.m, to calculate the wavefunctions in the radial direction 
   
   -sort the quantum states according to their eigen energy by running script_sort_e.m
   
2. Run script_main.m to calculate solve the density matrix equations, calculate the gain, and current density.
