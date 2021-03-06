# COMBINEDglucoseRASpodocytesACEInhibPKPD
combination of glucoseRASpodocytes and ACEInhibPKPD

A mathematical model of glucose-stimulated local renin-angiotensin system (RAS) in podocytes combined with a PKPD model of ACE inhibition of the RAS

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2550842.svg)](https://doi.org/10.5281/zenodo.2550842)

## Overview
The code for the model is provided as is, without guarantees and without support. The corresponding manuscript describing the model is currently under review for publication.

## Combined Glucose RAS in Podocytes ACE Inhibition in PKPD Model
### Authors
Minu R. Pilvankar, Hui Ling Yong, and Ashlee N. Ford Versypt, 
School of Chemical Engineering,
Oklahoma State University.
Corresponding author: A. N. Ford Versypt, ashleefv@okstate.edu

## Related Publication for Model Details
[1] Pilvankar, Minu R., Michele A. Higgins, and Ashlee N. Ford Versypt. "Mathematical Model for Glucose Dependence of the Local Renin–Angiotensin System in Podocytes." Bulletin of Mathematical Biology, 80, no. 4 (2018): 880-905.

[2] Pilvankar, Minu R., Hui L. Yong, and Ashlee N. Ford Versypt. "A Glucose-Dependent Pharmacokinetic/Pharmacodynamic Model of ACE Inhibition in Kidney Cells." Processes, 7(3), 131, 2019 https://doi.org/10.3390/pr7030131

### Main files
* MAIN.m.
   Runs the model manually either in the Editor or by passing arguments in the
   command line. Can specify plot_mode to produce and save publication-quality plots.

### Dependent & supplemental files

* combinedRAS_ACE_PKPD.m.
   The central file that uses case-specific parameters and dosing information 
   as input to the PKPD model of ACE inhibitor dose impact on AngII 
   concentration in podocytes. Output includes the vectors of time and the concentrations of 
   the drug diacid form, AGT, AngI, AngII, and Renin and the percent of inhibition 
   (not in this order). Note: do not run this file directly. Instead let it be called by call_combinedRAS_ACE_PKPD.m
   (described subsequently) that correctly defines all the input values.

* .mat files contain parameters described in the paper.
   
* call_combinedRAS_ACE_PKPD_.m.
   Used by modifiedrun_PKPD_without_GUI.m. 
   Sets up parameters before calling combinedRAS_ACE_PKPD.m for a scalar 
   value of tfinal_dosing. Output is all of the output from 
   combinedRAS_ACE_PKPD.m with options for plotting and saving figures
   through plot_mode string.
   
* .mat files. 
   Needed to run the programs to pass data, parameters, and calculated values.
   
* NormalCurvefitGlucoseDynamics.m.
 Polyfit to representative normal and diabetic subject glucose dynamic profiles.

### Generating plots for paper [2]
* script_generate_plots.m
Runs the local and systemic RAS model files to generate plots in the paper. Takes in figure number as an input (3-10) corresponding th the figure numbers in the paper.

* .mat files. 
   Needed to run the programs to pass data, parameters, and calculated values.

(c) Minu R. Pilvankar, Hui Ling Yong, and Ashlee N. Ford Versypt, 2018
