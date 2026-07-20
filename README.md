[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15485446.svg)](https://doi.org/10.5281/zenodo.15485446)

# SkelMuscleCa
MATLAB code for an ODE-based model of calcium dynamics and force production in skeletal muscle fibers
This repository is associated with the manuscript, "Systems modeling reveals that store-operated calcium entry modulates force and fatigue during exercise", currently accessible as a preprint [on bioRxiv](https://www.biorxiv.org/content/10.1101/2025.05.22.655415v1.full).

## Installation and Dependencies:

After downloading the files from this repository, the user should ensure the following dependencies are satisfied:
- Sensitivity analysis relies on UQLab, which can be [downloaded for free](https://uqlab.com/download). The downloaded code must be added to the path in MATLAB to run Morris-associated files.
- Optimization uses the Optimization Toolbox and Global optimization Toolbox in MATLAB.
- To run sensitivity analysis and optimization in parallel, the Parallel Computing Toolbox is required.

## Summary of repository

All plots shown in main and supplementary figures of the manuscript are reproducible by running the associated sections of the script, "SkelMuscleCa_MakeFigs.m".
The function of all other supporting files are described briefly below:
- SkelMuscleCa_dydt: Main ODE function used to solve for state variables over time in the mitochondria-free model using ode15s.
- SkelMuscleCa_dydt_withMito: Main ODE function used to solve for all state variables over time for the model including mitochondria using ode15s.
- SkelMuscleCa_fitMitoParam: File used to fit mitochondrial related parameters to data from [Wei et al 2015](https://pubmed.ncbi.nlm.nih.gov/25963147/)
- SkelMuscleObj: Function that computes model predictions for 3 experiments and then computes the objective function as described in our manuscript. Also returns quantities of interest (QOIs) needed for Morris sensitivity analysis.
- SkelMuscleCa_EstimateDriver: Driver script for particle swarm optimization (PSO)-based parameter estimation.
- SkelMuscleCa_paramEst: Function called by the estimate driver script to initialize PSO. Parallelization options should be configured within this file.
- SkelMuscleCa_Morris: Driver script for Morris sensitivity analysis (requires UQLab)
- SkelMuscleCa_SAOutput.m: Function called by Morris driver script to return QOIs for different parameter sets.
- SkelMuscleCa_sweepSOCE: Script for sweeping across SOCE conductance and calcium sensitivity, as tested in Fig 8 and Figs S7-S8.

An additional function, "prettyGraph.m" is used throughout for custom graph settings and can be adjusted by individual users as desired.

Furthermore, several .mat files are included in the "Data" folder of the repository:
- Exptdata.mat: Contains cell vector "Expt" that contains all experimental data used in optimization.
- p0Struct.mat: Contains structure "p0Struct" with fields "data" and "names", giving the values and names associated with all parameters.
- yinit0.mat: Contains vector "yinit0" with starting guesses for initial values of each state variable in the model.
- MorrisResults.mat: Contains uq_analysis object "MorrisAnalysis" from Morris sensitivity analysis
- pMito.mat: Contains vector "pMito" with calibrated mitochondrial parameters
- pSol_allParam_withMito.mat: Contains vector "pSol" with final solution of normalized parameters from PSO (model with mitochondria).
- pSol_VOnly_withMito.mat: Contains vector "pVec" with estimate of parameters to which voltage was sensitive (from "voltage-only" fit, model with mitochondria)
- pSol_allParam.mat_noMito: Contains vector "pSol" with final solution of normalized parameters from PSO (model without mitochondria).
- pSol_VOnly_noMito.mat: Contains vector "pVec" with estimate of parameters to which voltage was sensitive (from "voltage-only" fit, model without mitochondria)

Furthermore, the data from [Wei et al 2015](https://pubmed.ncbi.nlm.nih.gov/25963147/) extracted using PlotDigitizer and saved as .csv files can be found in the following folders within "Data":
- csv_bound_free: ratio of bound-to-free calcium in associated experiments
- csv_ca_out: extra-mitochondrial calcium in associated experiments
- csv_mito_ca: mitochondrial free calcium measurements in associated experiments
- csv_voltage: mitochondrial membrane voltage in associated experiments

## Contributing guidelines
Please feel free to submit any issues to this Github repository or suggest modifications by submitting a pull request.