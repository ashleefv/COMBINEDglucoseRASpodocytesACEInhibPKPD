# COMBINEDglucoseRASpodocytesACEInhibPKPD
A mathematical model of glucose-stimulated local renin-angiotensin system (RAS) in podocytes combined with a PKPD model of ACE inhibition of the RAS

[![DOI](https://zenodo.org/badge/94033856.svg)](https://zenodo.org/badge/latestdoi/94033856)

## Overview
The code for the model is provided as is, without guarantees and without support. The corresponding manuscript describing the model is currently under review for publication.

## Glucose RAS in Podocytes Model
### Authors
Minu R. Pilvankar, Michele A. Higgins, and Ashlee N. Ford Versypt, 
School of Chemical Engineering,
Oklahoma State University.
Corresponding author: A. N. Ford Versypt, ashleefv@okstate.edu

## Related Publication for Model Details
M. R. Pilvankar, M. A. Higgins, A. N. Ford Versypt, Mathematical Model for Glucose Dependence of the Local Renin-Angiotensin System in Podocytes, submitted 2017.

### Main files

* param_estimation_Approach1and2.m. Runs the model to estimate the parameters using Approach 1 and Appraoch 2 (described in the manuscript). The simulation results compare the change in concentration of Angiotensin II (ANG II) with increasing glucose for Approach 1 and different scenarios of Approach 2. The results are compared using root mean squared error (RMSE) using a data set from literature.
* param_estimation_Approach3.m. Runs the model to estimate the parameters using Approach 3 (described in the manuscript). The simulation results compare the change in concentration of ANG II with increasing glucose for Approach 3.
* Sensanalysis.m. Runs local sensitivity analysis on ANG II with respect to all the parameters. The resulting plots show sensitivity of all the parameters at normal and high glucose states.
* globalsensitivity.m. Runs global sensitivity analysis on ANG II with respect to all the parameters. The resulting plots show two sensitivity indices S_i (ﬁrst-order index) and S_{Ti} (total-order index which accounts for higher-order and non-linear interactions between the parameters) for ANG II with respect to each input parameter.

### Supplementary files
 
* glucoseRASssA12.m. Runs through param_estimation_Approach1and2.m and returns ANG II concentration by solving the set of equations. It is used to pass different set of coefficients for each case of Approach 1 and 2.
* glucoseRASssA3.m. Runs through param_estimation_Approach1and2.m and returns ANG II concentration by solving the set of equations. It is used to pass different set of coefficients for Approach 3.
* glucoseRASssSens.m. Runs through param_estimation_Approach1and2.m and returns ANG II concentration by solving the set of equations. It is used to pass the parameter values at NG and HG with a devaition to do the sensitivity analysis.

 Needed to run the programs to pass data, parameters, and calculated values.
* scenario0.m
* scenario1.m
* scenario2.m
* scenario3.m
* scenario4.m
* scenarioM0.m.
* Approach3.m
* NG.mat
* NGHG.mat

 For running globalsensitivity.m
* efast_sd.m
* efast_ttest.m
* parameterdist.m
* Parameter_settings_eFAST.m
* SETFREQ.m
  
    
* export_fig.m
   MATLAB package to nicely export figures.
   Download: https://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig
   Tips for usage: https://github.com/altmany/export_fig/blob/master/README.md

(c) Ashlee N. Ford Versypt, 2017
