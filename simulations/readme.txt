This folder contains simulation code and results for Adjustment Sets Selection for Estimating
Optimal Treatment Rules under Confounding by Galanter, Shortreed, and Moodie.

Contents:

truth.RDS
 true optimal rule values with tailoring variables 1 and 2 for each simulation setting

code
 - causal_ball_script_tailoring.R
    Code implementing causal ball from GitHub repository, as was available Jan 1, 2023, referenced in the supplement to Tang et al. (2023) , modified to always select tailoring variable
 - dips_tailoring.R
    Code implementing DiPS from supplement to Cheng et al. (2020), modified to always select tailoring variable
 - GLiDeR_Rfunctions_tailoring.R
    Code implementing GLiDeR from supplement to Koch et al (2018), modified to always select tailoring variable
 - CBPSmod
   Modified version of CBPS package version 0.23 to always select tailoring variable
 - simulation_functions.R
    functions used to run all the simulations, sources the files above
 - pre_specified_simulation_code.R
    code for running the simulations with pre-specified adjustment sets, sources simulation_functions.R
 - selection_simulation_code.R
    code for running the simulations with selection, sources simulation_functions.R
 - true_values.R
    code for generating the true optimal rule values with tailoring variables 1 and 2
 - pre_specified_results.R
    code to get figures and text documents of latex tables for pre specified simulations
 - selection_results.R
    code to get figures and text documents of latex tables for selection simulations

pre specified results
 - text files with latex result tables
 - pdf files with rule value estimation boxplots
 - jpg files with boxplots of propensity score distribution

selection results
 - text files with latex result tables
 - pdf files with rule value estimation boxplots
 - jpg files with boxplots of propensity score distribution
 - pdf files with variabple selection proportion scatterplots

 
