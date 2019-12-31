# lts
Code for analyzing lower tropospheric stability. This package provides the code required to reproduce the results of my paper. It is set up on the assumption that the code will be run on NCAR'S Cheyenne computer.

Data used:
* CESM Large Ensemble (CESM1-CAM5)
* CESM2-CAM6
* CESM2-WACCM
* ERA-Interim
* ERA5

## Components
# Data preparation
Data access is through the Glade server via the `collections` and `rda` directories. Filenames, variable names, and data types vary by data source. I use the variable names from the Community Earth System Model as the starting point, and use helper functions internally to translate into the names used in each source. 

`get_filenames`

# Analysis
I focus on the relationships between lower tropospheric stability, sea ice, and low-level clouds. 

`taylor_diagrams`  
`monthly_distributions`  
`internal_variability`  
`interannual_variability`  

