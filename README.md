# Finemap-MiXeR: A variational Bayesian approach for finemapping genomic data



This repo is to run Finemap Mixer 

Required inputs:

1) z-scores of the study 

2) A matrix (correlation matrix of SNPs weighted by Heterozygozity, ( $$ a_{ij}=\sum_{i=1}^{M} \sqrt{N \hat{H}_{i}} \hat{r}_{j i} $$).   For details see[here](https://www.biorxiv.org/content/10.1101/2022.11.30.518509v2.full.pdf))

Note: Instead of the weighted correlation matrix A, correlation matrix itself (R) can also be used but is not guaranteed for the optimal results.

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

# How to run the Finemap MiXeR with Singularity container


 
 If you do not have MATLAB licence, you can still run Finemap-MiXeR using Matlab Runtime and Singularity Container. Our method has been compiled using Matlab 2018B compiler  and  we can run this compiled application using Matlab Runtime that is hosted in a singularity container. All this process has already done and if you follow the following steps you can run Finemap-MiXeR without using MATLAB.

## Step 1 Download ``matlabruntime.sif``,  container from [here](https://drive.google.com/file/d/1tJ14nauquF_GZg10gOB1Mj2EIsXXYbLC/view?usp=drive_link)

## Step 2. Download the compiled application of our analysis from here [here](https://drive.google.com/file/d/1CCLCR2FJX-feHxo9uHHz-6YDpecVGSCJ/view?usp=sharing) and unzip it. This folder will be our working directory

## Step 4 Run the container as  `singularity exec -B $PWD:/execute  matlabruntime.sif /execute/hapgenruntime 'chr21' `
