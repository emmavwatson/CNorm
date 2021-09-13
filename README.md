# CNorm
Custom copy number tumor cohort normalization for isolating effects of individual CNAs
#
First, call the ChrArm_CNA_freq_TCGA_puritycorrection.r function to generate purity-corrected copy number calls and CNA summaries.
#
Then use the output to call Arm_Level_CNorm_diffexp.r function to analyze differential expression associated with a CNA.
#
This function will perform a cohort subsampling to normalize for co-occuring CNAs when possible. 
To simply generate differential expression data without cohort subsampling and CNA normalization, just set cohort_subsample_minFrac = 0.
#
![github_figure_CNorm](https://user-images.githubusercontent.com/90458376/132859486-4df16b17-62ab-44ec-a572-65b4dc78b804.png)
