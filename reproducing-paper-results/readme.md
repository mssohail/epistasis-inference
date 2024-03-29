
This repository contains codes and data required to reproduce results in "Inferring epistasis from genetic time-series data" doi:

# System and software requirements

The codes were written on MATLAB R2017b running Windows 10 64-bit. 

Requires MATLAB 2017b or later (may work on earlier versions but not tested). Also
requires Statistics and Machine Learning Toolbox (v11.2) or later.

# Downloading and running codes

1. Download the MPL codes from github
2. Download data files from [Zenodo](https://doi.org/10.5281/zenodo.6605512).
3. Unzip folders Analysis, Data and RandPermFiles within the folder containing
   epistasis-inference codes.
4. **Edit the dirName.txt file** within the MPL codes folder and **enter complete paths 
   of the Data, the Analysis, and the Figure folders.**
   
   example:
   - dirNameData=D:\epistasis-inference\Data\
   - dirAnalysisName=D:\epistasis-inference\Analysis\
   - dirFiguresName=D:\epistasis-inference\Figures\


# Steps to obtain results presented in the paper

## 1. Results shown in Figure 1 

   Estimating epistasis and selection coefficient under scenarios 
   (a) positive epistasis
   (b) negative epistasis
   (c) sign epistasis
   (d) no epistasis

### Steps: 
1. Run AnalysisMPL_Epi_2sites_v2_allSets_8cases.m: This code estimates the selection coefficients and epistasis terms for the two-locus system analyzed in Figure 1 of the paper. In Line 30, set analyzeDataForFigure = 1
    
2. Run PlotMPL_Epi_FF1_8cases.m: This code plots Figure 1

## 2. Results shown in Figure 2

   Show 2 instances of the data analyzed in Step 1.1. 
   
### Steps:

1. Run PlotMPL_Epi_FF2_TrajExample_new.m: This code plots 1 iteration each of two cases analyzed in Step 1.1.

## 3. Results shown in Figure 3

   Effect of data variability on estimation performance of MPL
   
### Steps:

1. Run AnalysisMPL_Epi_2sites_v2_allSets_8cases.m: This code estimates selection 
   coefficients and epistasis terms for a two-locus system analyzed in Figure 1 
   of the paper. 
   In Line 30, set analyzeDataForFigure = 3
   
2. Run PlotMPL_Epi_FF3_BoxPlotEstnSums_2.m: This code plots Figure 3

## 4. Results shown in Figure 4

   Performance of MPL as a function fo data variability controlled by varying the 
   number of genotypes in initial population
   
### Steps:   

1. Run AnalysisMPL_Epi_5sites_WithR_reg2_v2.m: This code estimates selection coefficients and epistasis terms for a five-locus system analyzed in Figure 4 of the paper. In Line 28, set numStrainsInInitialPop to 5, 10 and 20 to run for each data variability case one by one.
2. Run PlotMPL_Epi_5sites.m: This code takes the analysis of the previous step and prepares it for plotting. In Line 37, set numStrainsInInitialPop to 5, 10 and 20 to run for each data variability case one by one.
3. Run PlotMPL_Epi_FF4_NumStrains_4.m: This code plots Figure 4
     
## 5. Results shown in Figure 5

   Performance of MPL as a function fo data variability controlled by varying the 
   number of replicates used for inference

### Steps:   

1. Run AnalysisMPL_Epi_5sites_rep.m: This code estimates selection
   coefficients and epistasis terms for a five-locus system analyzed in Figure 5 of the paper.
   
   For running the case of 3 replicates
   
	 - In Line 27, set numReps = 3
	 - In line 28, set numItr1 = 999
	 
   For running the case of 5 replicates
   
	 - In Line 27, set numReps = 5	 
	 - In line 28, set numItr1 = 1000
	 
2. Run PlotMPL_Epi_5sites_rep_new.m: This code takes the analysis of the previous step and prepares it for plotting.
   
   For running the case of 3 replicates
   
	 - In Line 48, set thisItrEnd = 3	 
	 - In line 28, set lastItr = 999

   For running the case of 5 replicates
   
	 - In Line 48, set thisItrEnd = 5
	 - In line 28, set lastItr = 1000
	 
3. Run PlotMPL_Epi_FF5_NumRep_4.m: This code plots Figure 5	 
	

## 6. Results shown in Figure 6

   Compute performance of MPL by varying number of samples n_s and time sampling 
   step \delta_t
   
### Steps:

1. Run AnalysisMPL_Epi_5sites_HeatMap_reg2_v2.m: This code estimates selection
     coefficients and epistasis terms by varying n_s and \delta_t
2. Run PlotMPL_Epi_5sites_Heatmap_ng_dT.m: This code takes the analysis of the
     previous step and prepares it for plotting.
3. Run AnalysisMPL_Epi_5sites_HeatMap_reg2_v2_Tused.m
4. Run PlotMPL_Epi_5sites_Heatmap_Tused.m
3. Run PlotMPL_Epi_FF6_MPLE_heatmap_withTused.m: This code plots Figure 6	 

## 7. Results shown in Figure 7

   Comparision of MPL and MPL (no epistasis) over fitness landscape with varying 
   fraction of non-zero epistasis terms.

### Steps:
1. Run AnalysisMPL_Epi_5sites_reg2_v2_loop_Rev1_Sims.m (this code analyzes the various FL scenarios in Figure 7A.
2. Run AnalysisMPL_Epi_5sites_rep_loop.m 
    Change the variable 'thisSetAll' to thisSetAll = [1068141 1068144 1068241 1068244 1068143 1068142 1068243 1068242]
                        'numStrainsInInitialPop' to numStrainsInInitialPop = 5;
3. Run PlotMPL_Epi_5sites_AUROC_Rev1_Sims.m
4. Run PlotMPL_Epi_5sites_rep_new_AUROC_Rev1_Sims_loop.m
5. Run AnalysisMPL_Epi_5sites_rep_loop.m with thisSetAll = []
6. Run PlotMPL_Epi_FF7_modelComp_5site_lowDiv_rep.m:  This code plots Figure 7

## 8. Results shown in Supplementary Figure 1
   
   Effect of data variability on estimation performance of MPL
   
### Steps:

1. Run PlotMPL_Epi_FFS_boxplots_supp.m: This code plots Supplementary Figure 1

## 9. Results shown in Supplementary Figure 2
   
   Fraction of acessible, partiallly accessible and non accessible 
   selection coefficients and epistrasis terms as a function of varying data 
   variability analyzed in Figure 4.

### Steps:

1. Run PlotMPL_Epi_FFS2_fracAccess.m: This code plots Supplementary Figure 2

## 10. Results shown in Supplementary Figure 3
   
   Performance of MPL on sparse fitness landscapses analyzed in Figure 7

### Steps:

1. Run PlotMPL_Epi_FFS3_Landscapes_new.m: This code plots Supplementary Figure 3


## 11. Results shown in Supplementary Figure 5
   
   Performance comparison of MPL and MPL (no epistasis) on low variability data 
	Analyzed in Figure 4

### Steps:

1. Run PlotMPL_Epi_FFS5_ModelComp_LowVarib.m: This code plots Supplementary Figure 5
