# letermovir_cea
Cost-Effectiveness Model of Letermovir Prophylaxis vs Preemptive Therapy for Cytomegalovirus 

This repository contains a cost-effectiveness model assessing letermovir prophylaxis as compared to
pre-emptive therapy alone for cytomegalovirus in allogeneic hematopoietic stem cell transplant recipients. 

Notes on the model:

-Calculated inputs sheet: This sheet contains manual calculations for some of the costs incorporated into the model. Cells shaded in light blue can be altered by 
the user. References for the default inputs used in this model can be found in the accompanying model manuscript. 

-Parameters sheet: Contains input values to be used for the model. Column B contains the 'live values' that are currently being used by the model. Column C contains 
values used for the probabilstic sensitivity analysis (PSA). Column D contains values used for the base case/deterministic analysis. To generate results for the PSA, 
click the "Run PSA" button. For the accompanying CEA curve and ICER plane, click the "Generate CEA Curve" button. 

- Decision Tree sheet: Contains results for the decision tree portion of the model. The user should not alter any cells on this sheet. 

- Markov Model_b sheet: Contains results for the Markov model portion of the model, in addition to final results. The user should not alter any cells on this sheet. 

- OWSA: Contains results of the one-way sensitivity analyses. Note: The OWSA must be done manually (e.g. inputting each value into the model at a time to generate results). 
The results are not automatically generated. 

- PSA: Contains results of the PSA; the accompanying cost-effectiveness acceptability curve and incremental cost effective plane can be viewed at the far right of the sheet. 

- Survival Extrapolation: Contains results related to extrapolation of a previously published survival curve of allo-HCT recipients. Note: the initial parameters for the distributions 
(e.g. Columns B-F) were generated using R. Sensitvity analyses related to changing the curves using cholesky decomposition matrices can also be viewed here. 
