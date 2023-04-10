# arkansas_size_spectra
This repo contains code for the time-series analysis of macroinvertebrate size spectra in the Upper Arkansas River, CO

Site AR1 is upstream and is used as a reference. Although it was subject to some mining inputs from the Leadville Mine Drainage Tunnel (LMDT) Until ~ 1993. 

Site AR3 is immediately downstream of California Gulch, an EPA super fund site. It has received extensive restoration activities since the mid 1990's. 

Raw data files take up too much storage space and are not being stored on github. Script 01 stitches and cleans all the raw CSV files and saves the output as an r data file (.RDS). 

Script 02_est_dw... takes the stitched and cleaned data, attaches taxon-specific length weight coefficients and estimates dry weight. 

Script 03_lambda... estimates the exponent of a power law describing the decline in abundance with increasing body size. 