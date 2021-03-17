# R-scripts_MA

All skripts were built under R for Windows 3.6.1 and can be run directly without any changes as long as the working directory is set to the folder containing the skript to be run. 

## Dependencies
* FactoMineR v.2.4
* permute v.0.9-5
* uwot v.0.1.10
* vegan v.2.5-7

## PCA+CI 
This folder contains the script for a principal component analysis on the exemplary timeseries data set in “rawdata.csv” where 95% confidence intervals are added around the centroids of each week and cell population applied on both a correlation (cor) and covariance (cov) matrix and the centroids of one treatment are connected from week to week resulting in a time trajectory. Exemplary result graphs are included in the png-files “2021-02-24.pca.cor.T2.png” and “2021-02-24.pca.cov.T2.png”.
## dummydata.generator
This folder contains a script to add supplementary data to a timeseries data file prior to a principle response curve analysis. By opening and running "dummydata.generator.v1.R" the csv-file “rawdata.csv” is supplemented with rows of artificial data points resulting in the same number of measurements for each week and cell population and an ID-column is added. The new data set is then t-tested against the original data set per week and treatment to ensure that there is no significant difference in between original and new data. The new table is saved as csv-file for further analysis (cf. „2021-02-24_pcr.data+dummy.csv”) and the results of the t-test as a txt-file (cf. “2021-02-24 _res.t.test.raw+dummy.txt”).
## PRC
In this folder scripts for a PRC on the raw data without statistics as well as a PRC on a supplemented dataset with statistics can be found. Both scripts can be run on the included exemplary data file “pcr.data+dummy.csv”. The script on the raw data “prc_raw.v3.R” creates only a graph as found in “2021-02-24_prc_raw.cor.png” and the results of the PRC in “2021-02-24_prc_raw.cor.res.txt”. The script on the supplemented data set “prc.v3.R” calculates the PRC with statistics ten times as the ID-vector is permutated ten times, creates a graph of each permutation (cf. “2021-02-24_prc.cor.perm<I>x</I>.png”) and summarizes the results of each PRC and the permutation test for significance in the txt-file “2021-02-24_prc.cor.res.txt”. 
## UMAP
The code for the UMAPs can be found in this folder. The script “umap.2d.R” can be run on the included file “Example data.csv” and performs a 2D UMAP, calculates a spline regression between biomarker value  and UMAP-axes and plots the results once by differently colored cell populations and once by the heatmap coloring of each variable. Results can be compared to the exemplary result files „2021-02-24 _Example UMAP.png” and “2021-02-24_Heatmap+var.vec.x.png”. The “umap.3d.v1.R”-Script calculates a 3D UMAP but does not include the heatmap-coloring and plots the UMAP by coloring the different cell populations only (cf. “2021-02-24_3d.umap+var.vec.html”).
