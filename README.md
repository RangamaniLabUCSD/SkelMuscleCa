# SkelMuscleCa
MATLAB code for ODE-based model of calcium in skeletal muscle. The model was assembled using VCell software, and the resulting ODE model was exported to MATLAB, then altered further.

Files summary:

- SkelMuscleCa.vcml: VCell model used for development
- SkelMuscleCa_dydt : Main ODE function to solve for Ca2+ in skeletal muscle.
- SkelMuscleCa_Result: Script for generating calcium and force output data.
- SkelMuscleCa_dydtEst: ODE function to solve for Ca2+ during parameter estimation.
- SkelMuscleCa_EstimateDriver: Driver script for parameter estimation.
- SkelMuslceCa_paramEst: Function for computing objective value of parameter estimation.
- SkelMuslceCa_sweepScript: Script for sweeping across the parameter range for estimation. 
- SkelMuscleCa_Force: Script fot fitting Ca2+ to Hill model equation of Ca-Force relationship.
- SkelMuscleCa_Plots: Script for generating all the plots.
- SkelMuscleCa_Morris: Script for Morris sensitivity analysis.
- SkelMuscleCa_SAOutput.m: Function for computing quantities of interest for sensitivity analysis.

Dependencies:

Our sensitivity analysis relies on UQLab, which can be [downloaded for free](https://uqlab.com/download).
Our optimization uses the Optimization Toolbox and Global optimization Toolbox in MATLAB.
To run sensitivity analysis and optimization in parallel, the Parallel Computing Toolbox is required.
