# HOIs-steepen-BEF-relationships
Code and data used for the MS " HOIs among trees steepen the BEF relationship

System information:
Model fitting was submitted to EVE cluster (VERSION="9.5 (Blue Onyx)" ID="rocky")
GCC version 12.2.0; 
OpenMPI version 4.1.4;
rstan (Version 2.21.8, GitRev: 2e1f913d3ca3);
R version 4.2.2

All the analyses and simulations were conduct in
R version 4.3.3

## Descriptions of code and data files
1. Code/model_array.R: Code to fit the model candidates with specs (i.e. iteration, acceptance).
2. Code/model.stan: A folder contains all the model formulations tested in the study (18 in total).
3. Code/array_job.sh: For submitting array jobs to HPC
4. Code/IS-div-analysis.R: Code to analyze the results of the pairwise-nonlinear model (mod.10 in the model.stan folder), and produce Fig.2 in the MS. 
5. Code/sim_growth.R: Code to simulate tree growth over 7 years with and withoug net HOI effects, and produce Fig.3 A and B in the MS. The simulation output sim_output.csv (due to large size) can be found at https://zenodo.org/records/16911388.
6. Code/SEM_IS.R: Code to extract interactions from each plot and SEM analysis.
7. Code/functions.R: Functions required for simulating tree growth and extracting tree-tree interaction.
8. Data/params.csv: Parameter files for running the array job on HPC. Files of 1,2, and 8 are needed to run the models on HPC.
9. Data/stan_dat_1.RData: Data for fitting the model.
10. Data/pars_pairwise_nl.csv: Parameters estimated by the pairwise-nonlinear model. For script 4.
11. Data/year_7_growth.csv: The growth data of the seventh year. For script 5.
12. Data/BEF_plot_map.csv: The configuration of each plot. For script 6.


