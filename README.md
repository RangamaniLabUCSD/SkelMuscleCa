# SkelMuscleCa
MATLAB code for ODE-based model of calcium in skeletal muscle. The model was assembled using VCell software, and the resulting ODE model was exported to MATLAB, then altered further.

Files summary:
- SkelMuscleCa.vcml: VCell model used for development
- SkelMuscleCa_dydt : Main ODE function to solve for Ca2+ in skeletal muscle.
- SkelMuscleCa_Result: Script for generating calcium and force output data. Dependencies: Parallel Computing Toolbox
- SkelMuscleCa_dydtEst: ODE function to solve for Ca2+ during parameter estimation.
- SkelMuscleCa_EstimateDriver: Driver script for parameter estimation. Dependencies: Global Optimization Toolbox, Optimization Toolbox and Parallel Computing Toolbox
- SkelMuslceCa_paramEst: Function for computing objective value of parameter estimation.Dependencies: Global Optimization Toolbox, Optimization Toolbox and Parallel Computing Toolbox
- SkelMuslceCa_sweepScript: Script for sweeping across the parameter range for estimation. 
- SkelMuscleCa_Force: Script fot fitting Ca2+ to Hill model equation of Ca-Force relationship. Dependencies: Optimization Toolbox
- SkelMuscleCa_Plots: Script for generating all the plots.
- SkelMuscleCa_Morris: Script for Morris sensitivity analysis. Dependencies: UQLab, Parallel Computing Toolbox
- SkelMuscleCa_SAOutput.m: Function for computing quantities of interest for sensitivity analysis. Dependencies: UQLab, Parallel Computing Toolbox
