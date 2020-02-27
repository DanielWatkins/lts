"""Apply xESMF to regrid sea ice to the atmosphere grid.
Requires the parameter object to be sent to it somehow. Temporary 
text files are probably the way to go for now.

Depending on workflow, could make a commandline function with argparse.
"""

# TODO: see if it works to reshape and make a new dataset for cases where there
# are missing lat/lons

import numpy as np
import xarray as xr
import xesmf as xe
import pandas as pd


info = pd.read_csv('data_for_gridding.tmp')
ref_file = info.loc[info.name == 'reference','info'].values[0]
save_loc = info.loc[info.name == 'save_location', 'info'].values[0]
min_lat = int(info.loc[info.name == 'minimum_latitude', 'info'].values[0])
varname = info.loc[info.name == 'varname','info'].values[0]
reference_dataset = xr.open_dataset(ref_file)
new_lat = np.round(reference_dataset['lat'], 4)
new_lon = np.round(reference_dataset['lon'], 4)
reference_dataset.close()

ds_out = xr.Dataset({'lat': (['lat'], new_lat[new_lat >= min_lat]),
                     'lon': (['lon'], new_lon),
                    }
                   )

for targ_file in info.loc[info.name == 'target', 'info']:
    target_dataset = xr.open_dataset(targ_file)[varname]
    fname = targ_file.split('/')[-1].split('.')[0] + '_int.nc'
    regridder = xe.Regridder(target_dataset, ds_out, 'bilinear', reuse_weights=True)
    ds_regridded = regridder(target_dataset)
    ds_regridded.to_netcdf(save_loc + fname)
    target_dataset.close()
    
regridder.clean_weight_file()