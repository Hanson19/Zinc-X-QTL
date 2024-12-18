# Zinc-X-QTL
Data and Code for Zinc Toxicity X-QTL mapping and phenotyping.

All scripts beging with an introduction which includes a brief description of what the script does, what packages are required, and what data files are required to run the script. The scripts are then divided into sections with names that briefly explain what each section does. If using RStudio you can see each section's name by clicking the document outline button at the top right corner of the scripts panel. 

These scripts can be devided into two main sections:

X-QTL Experiment and Analysis and Phenotyping Analysis
## X-QTL Experiment
These scripts will go through X-QTL analysis including estimating the selective pressure of 25mM ZnCl2, identifying QTL using DSPR founder haplotype frequencies, and making figures included in paper and supplemental matierals.

Scripts Include:

XQTL Replicate Collection Data.R : Estimates the selective pressure of 25mM ZnCl2 and estimating the number of females tested per replicate of XQTL mapping.

Zinc Analysis.R : Identifies QTL in genetic distance using the DSPR founder haplotype frequencies, and creates QTL Plot and founder haplotype frequencies changes plot. (This analysis must be run before being able to run any of the remaining scripts)

XQTL Duplicate Control Analysis.R : Compares duplicate control pools collected from the replicates. 

XQTL Founder Frequency Visualization.R : Creates supplemental figures that plots the actual DSPR founder haplotype frequencies, not just the change.

XQTL Location Homeostasis Proteins.R : Plots known zinc homeotasis proteins (ZnTs, Zips, MTs, and MTF-1) agaisnt backdrop of our identified QTL. 

## Phenotyping Analysis
These scripts analyze phenotyping data from X-QTL Subpopulations, Semi-Inbred Lines, Phenotyping Cohort, and RNAi KD of candidate genes. There are three main phenotypes that are measured: Zinc naive female adult resistance to zinc, developmental resistance to heavy metals measuring the percentages of eggs that develope to adulthood, developmental resistance to heavy metals measuring how long for eggs to reach adulthood. 

X-QTL Subpopulations:

XQTL Subpopulations Adult Resistance.R : Female resistance on 100mM ZnCl2

XQTL Subpopulations Developmental Resistance.R : Resistance on 25mM ZnCl2 and on water. Plots figures 3A and 4A.

Semi-Inbred Line Analysis Scripts:

Semi-Inbred Line Analysis.R : Developmental resistance of semi-inbred lines on a variety of heavy metals

Phenotyping Cohorts:

Phenotyping Cohorts Adult Resistance.R : Female resistance on 100mM ZnCl2

Sup_Fig_5_Adult_Resistance.R: Plotting adult zinc resistance from the X-QTL Subpopulations and Phenotyping Cohorts. Need tables made in X-QTL and Phenotyping Cohorts Adult Resistance.R

PP-II Developmental Analysis.R : Developmental resistance of flies on a variety of heavy metals.

RNAi Candidate Genes:

RNAi Analysis.R : Developmental resistance of candidate gene RNAi KDs on 10mM ZnCl2. 

