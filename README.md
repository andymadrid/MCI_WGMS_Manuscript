# MCI_WGMS_Manuscript
Code/scripts for the analysis and figure generation for the manuscript comparing cognitively unimpaired (CU), mild cognitive impairment (MCI), and Alzheimer's disease (AD)

This repository is broken down into a few different sections. One being "Analysis" where scripts take raw methylation values, filter low performing CpGs, impute missing values, generate methylation/coverage matrices, perform differential methylation analysis, compute estimate methylation differences, and (finally) correct for local FDR.

Another section ("Figures") has code used to make the vairous figures from the manuscript (e.g., volcano plot) and to generate data for supplementary tables.

In each section there will be near-identical copies of the scripts that were used for each of the different comparions that were made in the manuscript. These include MCI:CU, AD:MCI, AD:CU, CU:MCI:AD. I'm sure there was a nice, clean way to combine all of the different iterations of the scripts into one concise, cohesive one, but available time and energy led me down this path instead. 

-A.M.
