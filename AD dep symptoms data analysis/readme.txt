This folder contains the data analysis code for Adjustment Set Selection for Estimating Optimal Treatment Rules under Confounding by Galanter, Shortreed and Moodie.

Contents:

- causal_ball_script_tailoring.R
   Code implementing causal ball from GitHub repository, as was available Jan 1, 2023, referenced in the supplement to Tang et al. (2023) , modified to always select tailoring variable
- dips_tailoring.R
   Code implementing DiPS from supplement to Cheng et al. (2020), modified to always select tailoring variable
- GLiDeR_Rfunctions_tailoring.R
   Code implementing GLiDeR from supplement to Koch et al (2018), modified to always select tailoring variable
- CBPSmod
   Modified version of CBPS package version 0.23 to always select tailoring variable
- analysis_code.R
   Code to conduct the data analysis, sources the files above
- analysis_summary
   Code to summarize the analysis results into figures
