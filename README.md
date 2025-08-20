# HOIs-steepen-BEF-relationships
Code and data used for the MS " HOIs among trees steepen the BEF relationship

## Descriptions of code and data files
1. Code/model_array.r: Code to fit the model candidates with specs (i.e. iteration, acceptance).
2. Code/model.stan: A folder contains all the model formulations tested in the study (18 in total).
3. Code/params.csv: Parameter files for running the array job on HPC. Files of 1,2,3,--- are need to run the models on HPC.
4. Code/IS-div-analysis.r: Code to analyze the results of the pairwise-nonlinear model (mod.10 in the model.stan folder), and produce Fig.2 in the MS. 
5. Code/sim_growth.r: Code to simulate tree growth over 7 years with and withoug net HOI effects, and produce Fig.3 A and B in the MS.
6. Code/SEM_IS.r: Code to extract interactions from each plot and SEM analysis.
7. Code/functions.r: Functions required for simulating tree growth and extracting tree-tree interaction.
8. Data/stan_dat_1.RData: Data for fitting the model.
9. Data/pars_pairwise_nl.csv: Parameters estimated by the pairwise-nonlinear model. For script 4.
10. Data/year_7_growth.csv: The growth data of the seventh year. For script 5.
11. Data/BEF_plot_map.csv: The configuration of each plot. For script 6.


