"""Master script to walk through the data analysis."""
from lts.utilities import *
from lts.computation import *

#### Parameters ####
data_save_loc = '~/data/'
comp_save_loc = '~/computed/'
image_save_loc = '~/images/'

variables = ['ICEFRAC', 'LANDFRAC', 'TREFHT', 'T850', 'CLDLOW', 'CLDTOT', 'TGCLDLWP']

#### Grouping Data ####
get_data(variables, source='CESMLE', saveloc=data_save_loc)
get_data(variables, source='CESM2-CAM6', saveloc=data_save_loc)
get_data(variables, source='CESM2-WACCM', saveloc=data_save_loc)
get_data(variables, source='ERA-I', saveloc=data_save_loc)
get_data(variables, source='ERA5', saveloc=data_save_loc)

#### Computations ####
compute_time_series(variables, dataloc=data_save_loc, saveloc=comp_save_loc)
compute_ensemble_variability(variables, dataloc=data_save_loc, saveloc=comp_save_loc)
compute_interannual_variability(variables, dataloc=data_save_loc, saveloc=comp_save_loc)

#### Visualizations ####
plot_annual_cycle(['LTS', 'TREFHT', 'T850', 'CLDLOW', 'CLDLOW'], dataloc=data_save_loc, saveloc=image_save_loc)
plot_correlation(x='ICEFRAC', y='LTS', dataloc=data_save_loc, saveloc=image_save_loc)
plot_correlation(x='ICEFRAC', y='TREFHT', dataloc=data_save_loc, saveloc=image_save_loc)
plot_correlation(x='ICEFRAC', y='T850', dataloc=data_save_loc, saveloc=image_save_loc)
plot_taylor(variables, dataloc=data_save_loc, saveloc=image_save_loc)
                 
map_annual_cycle(['LTS', 'TREFHT', 'T850', 'CLDLOW', 'CLDLOW'], dataloc=data_save_loc, saveloc=image_save_loc)
map_ensemble_variability(['LTS', 'TREFHT', 'T850', 'CLDLOW', 'CLDLOW'], dataloc=comp_save_loc, saveloc=image_save_loc)
map_interannual_variability(['LTS', 'TREFHT', 'T850', 'CLDLOW', 'CLDLOW'], dataloc=comp_save_loc, saveloc=image_save_loc)
