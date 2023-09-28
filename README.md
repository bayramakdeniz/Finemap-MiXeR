# Finemap-MiXeR: A variational Bayesian approach for finemapping genomic data



This repo is to run Finemap Mixer 

Required inputs:

1) z-scores of the study 

2) A matrix (correlation matrix of SNPs weighted by Heterozygozity, ($ a_{ij}=\sum_{i=1}^{M} \sqrt{N \hat{H}_{i}} \hat{r}_{j i} $)
3)
4)   for details see[here](https://www.biorxiv.org/content/10.1101/2022.11.30.518509v2.full.pdf))

Note: Instead of weighted correlation matrix A, correlation matrix itself can also be used but not guaranteed for 

# =================================
# How to run the analysis

Step by Step guide how to run Finemap Mixer


## Step 1. 

Download and unzip the repo. This will be our working directory in Matlab.

## Step 2. 
This step is optional and you can skip if you do not want to create files for Susie and Finemap. For this case you need to set "data_creator=0"

Create empthy folders under working directory as finemap/1causal_hapgen/ .... finemap/5causal_hapgen/  and             susie/1causal_hapgen/ ..... susie/5causal_hapgen/ folders. These 10 folders will be used to run the same analysis for SuSie and Finemap

## Step 3.

Run hapgenRun_clean.m

##
Once the experiments have been done, it is possible to run the same analysis using exactly same input data with Finemap and SuSie

For Susie, you may directly run ``susieFromMixerData.R`` or ``susieFromRealData.R``

ToDo: Add More details and Finemapv1.4

# How to run the analysis with Singularity container

ToDo: Update

It may not feasible to run this analysis directly on a HPC due to Matlab Licence issue. In those cases we can use Matlab Runtime if you have compiled your Matlab code as an application as expained in detail [here](https://github.com/comorment/matlabruntime). To do this;

## Step 1 Download ``matlabruntime.sif``,  placed in matlabruntime folder from [here](https://drive.google.com/drive/folders/1mfxZJ-7A-4lDlCkarUCxEf2hBIxQGO69?usp=sharing)

## Step 2. Download the compiled application of our analysis from here [here](https://drive.google.com/file/d/1CCLCR2FJX-feHxo9uHHz-6YDpecVGSCJ/view?usp=sharing) and unzip it. This folder will be our working directory

## Step 3. Move matlabruntime.sif and plink data you want to analyze to the working directory. In this case the name of the plink file is ``chr21``

## Step 4 Run the container as  `singularity exec -B $PWD:/execute  matlabruntime.sif /execute/hapgenruntime 'chr21' `
