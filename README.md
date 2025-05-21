[![DOI](https://zenodo.org/badge/603239582.svg)](https://doi.org/10.5281/zenodo.15485446)

# SkelMuscleCa
MATLAB code for ODE-based model of action potential generation, calcium handling, and crossbridge cycling in skeletal muscle fibers, associated with the manuscript, "Systems modeling reveals that store-operated calcium entry modulates force and fatigue during exercise".

Files summary:

- SkelMuscleCa_dydt : Main ODE function to solve for Ca2+ in skeletal muscle
- SkelMuscleCa_dydt: Function defining all 51 ODEs in the model
- SkelMuscleCa_EstimateDriver: Driver script for parameter estimation.
- SkelMuscleObj: Function for computing objective value of parameter estimation.
- SkelMuslceCa_sweepSOCE: Script for sweeping across the parameter range for Figs 7, S5, and S6 
- SkelMuscleCa_TestPlots: Script for generating all the plots shown in the manuscript.
- SkelMuscleCa_Morris: Script for Morris sensitivity analysis.
- SkelMuscleCa_SAOutput.m: Function for computing quantities of interest for sensitivity analysis.

Dependencies:

Our sensitivity analysis relies on UQLab, which can be [downloaded for free](https://uqlab.com/download).
Our optimization uses the Optimization Toolbox and Global optimization Toolbox in MATLAB.
To run sensitivity analysis and optimization in parallel, the Parallel Computing Toolbox is required.
