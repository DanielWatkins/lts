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

# Get parameters
info = pd.read_csv('data_for_gridding.tmp')
ref_file = info.loc[info.name == 'reference','info'].values[0]
save_loc = info.loc[info.name == 'save_location', 'info'].values[0]
min_lat = int(info.loc[info.name == 'minimum_latitude', 'info'].values[0])
varname = info.loc[info.name == 'varname','info'].values[0]
variable = info.loc[info.name == 'variable','info'].values[0]
source = info.loc[info.name == 'source','info'].values[0]
begin_time = info.loc[info.name == 'begin_time','info'].values[0]
end_time = info.loc[info.name == 'end_time','info'].values[0]

# Get reference values
reference_dataset = xr.open_dataset(ref_file)
new_lat = np.round(reference_dataset['lat'], 4)
new_lon = np.round(reference_dataset['lon'], 4)
reference_dataset.close()

ds_out = xr.Dataset({'lat': (['lat'], new_lat[new_lat >= min_lat]),
                     'lon': (['lon'], new_lon),
                    }
                   )

def make_filename(filename):
    parts = filename.split('_')
    for i in range(len(parts)):
        if len(parts[i]) >= 6: # Might need to do ==6 for CMIP5 and ==8 for CMIP6
            if parts[i][0] == 'r':
                ensemble = parts[i]
                break
                
    if i == len(parts):
        ensemble = ''
        
    return save_loc + '_'.join([source, ensemble, variable]) + '.nc'


def clean_dataset(ds):
    """Ensures that the dataset has time series, lat/lon values, etc that
    will play well with others: i.e., rounding lat lon so that no doubling 
    happens, making time series follow a calendar, etc."""
    import xarray as xr
    import numpy as np
    
        
    data = ds[varname].values
    lats = np.round(ds['lat'].values, 6)
    lons = np.round(ds['lon'].values, 6)
    times = pd.date_range(params.begin_time, freq='MS', periods=len(ds['time']))
    
    return xr.Dataset({variable: (('time', 'lat', 'lon'), data)}, 
               coords={'time': times, 'lat': lats, 'lon': lons})



for targ_file in info.loc[info.name == 'target', 'info']:
    target_dataset = xr.open_dataset(targ_file)[varname]
    fname = targ_file.split('/')[-1]
    regridder = xe.Regridder(target_dataset, ds_out, 'bilinear', reuse_weights=True)
    ds_regridded = regridder(target_dataset)
    ds_regridded = clean_dataset(ds_regridded)
    ds_regridded.sel(time=slice(begin_time, end_time)).to_netcdf(make_filename(targ_file))
    target_dataset.close()
    
regridder.clean_weight_file()



# Doesn't account for data split between multiple files