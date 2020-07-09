# LME-dust-drought
Supporting code for analyzes performed for the manuscript Dust-Drought Nexus in the Southwestern United States: A Proxy-Model Comparison Approach by S. H. Arcusa, N. P. McKay, C. M. Carrillo and T. R. Ault 

# Data source

LME data can be downloaded from http://www.cesm.ucar.edu/projects/community-projects/LME/data-sets.html
Paleo data can be downloaded from https://www.ncdc.noaa.gov/paleo/study/22454 (Living Blended Drought Atlas) and https://www.ncdc.noaa.gov/paleo/study/27078 (Blue Lake).

# Required data & pre-processing

For the model:
1) Data for each variable should be stored for each model run in folders called runXXX, where XXX refers to the model run from 002 to 011, for the years 850-2005
2) Variables to download include: SOILICE, SOILIQ, DSTFLXT, DSTDEP, TLAI, TSAI, FSNO, SOILWATER_10CM, U10
3) Follow method from manuscript to prepare the variables. Eg. calculate fv from TLAI and TSAI, calculate fm, and cube of U10
4) Surface data files are also necessary including: "surfdata_1.9x2.5_simyr1350_c131018.nc"

For the paleo:
5) For Blue Lake, the data must be converted to LiPD file from http://lipd.net/playground.
6) For the LBDA, the region of interest should be selected.

# HPC computing

Many of the analyses require HPC. The scripts are not currently written for this to work automatically. 

# Order of analysis

1) Fig1.R
2) Fig5.R
3) Part_1a.R*
4) Part_1b.R**
4) Part_2.R
5) Part_3.R
6) Part_4.R

Notes:

 * Script part_1.R needs to be run for each model ensemble member. Currently manual. One member analysis takes ~10hrs on HPC. The outputs must be manually saved and used in Part_1b. 
 ** Output of each variable matrix should be saved as .csv to be used in Part_2 following file name of coef_SOILWATER.ann.011.csv for example. Change coef to lmg for the variance explained as a fraction.
