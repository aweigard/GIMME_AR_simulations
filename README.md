# GIMME_AR_simulations

This repository contains the code used to implement the simulation/recovery study conducted in the manuscript "The Influence of Autoregressive Relation Strength and Search Strategy on Directionality Recovery in GIMME" by Weigard, Lane, Gates and Beltz (2020 preprint).

It contains a number of R scripts which cover the different stages of the simulation, model fitting, data extraction, recovery analysis, and plotting processes. The exact simulated data and GIMME model fits that were used in the study are available at https://osf.io/zvfuh/.

- simulate_data.R contains code used to simulate the data for all simulation conditions
- The files GIMME_AR_Balanced.R, GIMME_MS_Balanced.R, GIMME_AR_LagGreater.R, GIMME_MS_LagGreater.R, GIMME_AR_ConGreater.R, and GIMME_MS_ConGreater.R contain code used to fit GIMME-AR and GIMME-MS models to data from all simulation conditions using the "gimme" R package (https://cran.r-project.org/web/packages/gimme/index.html). They are meant to be run on a high-performance computing cluster where multiple cores can be used for parallel processing.
- The analyze_recovery_Balanced.R, analyze_recovery_LagGreater.R and analyze_recovery_ConGreater.R contain code used to extract relevent data from GIMME model fits, and compute recovery statistics using functions defined in the calc_functions.R file.
- The save_ms_summary_fit file can be used to extract relevant information about the number of multiple solutions produced by GIMME-MS for each fit to simulated data.
- final_analyses_and_plots.R generates the tables and plots that are used to report recovery values in the manuscript.

If there are any difficulties in reproducing the code in this repository please reach out to me at asweigar@med.umich.edu.
