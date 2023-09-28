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

Download plink files generated using hapgen2 tool from [here](https://github.com/comorment/containers/tree/main/reference/hapgen) and move these files to your working directory

## Step 3.

To run the analysis using this genomic data (with various scenarios) run


# How to run the Finemap MiXeR with Singularity container

 If you do not have MATLAB licence, you can still run Finemap-MiXeR using Matlab Runtime and Singularity Container. Our method has been compiled using Matlab 2018B compiler  and  we can run this compiled application using Matlab Runtime that is hosted in a singularity container. All this process has already done and if you follow the following steps you can run Finemap-MiXeR without using MATLAB. Assuming that singularity container is available on your machine, you can also follow the steps below to run the analysis

## Step 1
Download ``matlabruntime.sif``,  container from [here](https://drive.google.com/file/d/1tJ14nauquF_GZg10gOB1Mj2EIsXXYbLC/view?usp=drive_link)

## Step 2. 
Download the compiled application of our analysis from here [here](https://drive.google.com/file/d/1SoLpSclxm5NsGSEz2nWfN0CawuicWid7/view?usp=sharing).

## Step 3. 
Prepare z scores as a vector in .mat format (glmt.mat) and A matrix (a.mat)

## Step 4 

Move container, application and corresponding input (.mat) files in the same directory, you can run the analysis as

Run the container as  `singularity exec -B $PWD:/execute  matlabruntime.sif /execute/hapgenruntime 'a.mat' 'glmt.mat' `

Then the output is created as "ResultsFinemapMixer.mat" and it includes the probability of being causal (probs.mat), effect size of each SNP (effect_size.mat) and its variance (variance_effect.mat)
