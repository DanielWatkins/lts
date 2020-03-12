"""A collection of functions for analyzing lower tropospheric stability.



(c) Daniel Watkins 2020
"""

import numpy as np
import xarray as xr
import pandas as pd
import lts.utilities.utilities as util

# TODO
# Set up filenames to pull just first ensemble
# Set up models method to display all cmip5 and cmip6 model names
# Check whether naming conventions will work for all ensemble members
# Set up quantile function for area weighting


class Parameters:
    """Bundle of parameters used throughout analysis.
    source: One of 'cesm-le', 'cesm2-cam6', 'cesm2-waccm', 'era-i', 'era5'.
    lts_type: One of 't850-t2m' (default), 't925-t2m', 't850-t1000'.
    save_location: String, ending with /, default '~/data/'.
    begin_time: Start date for subsetting. Default '1979-01-01'.
    end_time: End date for subsetting. Default '2018-12-31'.
    variables: List of variables
    frequency: 'monthly' (only option right now)
    minimum_latitude: Default 60
    maximum_latitude: Default 90
    minimum_longitude: Default 0
    maximum_longitude: Default 360
    
    ensemble_members: 'all', 'first', 'random'
    """

    def __init__(self, source):
        self.source = source
        self.lts_type = 't850-t2m'
        self.save_location = '~/data/'
        self.begin_time = '1979-01-01'
        self.end_time = '2018-12-31'
        self.variables = ['air_temperature', '2m_temperature']
        self.pressure_levels = [1000, 925, 850]
        self.ensemble_members = 'all'
        self.frequency = 'monthly'
        self.minimum_latitude = 60
        self.maximum_latitude = 90
        self.minimum_longitude = 0
        self.maximum_longitude = 360
        
    def display(self):
        for x in self.__dict__:
            print(x + ':', self.__dict__[x])
        
    def filenames(self):
        """Returns a dictionary with keys for each variable, and
        a each entry is a list of the matching files on glade."""
        fn_dict = {}
        for varname in self.variables:
            fn_dict[varname] = util.get_filenames(
                varname, self)
        return fn_dict
    
    def ensembles(self):
        """Returns a list of the ensemble members available for source."""
        fnames = self.filenames()
        for key in fnames:
            break
        fnames = fnames[key]
        
        def ensemble_finder(fname):
            """Pulls out the ensemble from the file name"""
            components = fname.split('/')[-1].split('_')
            for comp in components:
                if comp[0] == 'r':
                    if comp[1].isdigit():
                        return comp
        return list(np.unique([ensemble_finder(f) for f in fnames]))
        
        
        
#### TO DO ######
    
    def list_available_variables(self):
        """Returns a list of available variables for the source"""
        
        return
    
    def copy():
        """Fully copies parameters object"""
        return
    
    def cmip5_list():
        """Returns a list of all the CMIP5 models available."""
    
        return

    def cmip6_list():
        """Returns a list of all the CMIP6 models available."""
        return

    def reanalysis_list():
        """Returns a list of all the reanalysis models available."""
        return


def get_data(params):
    """Wrapper for load_dataset. Checks to see if required variables
    are included based on the desired lts_type.
    """
    
    filenames = params.filenames()
    
#    if params.lts_type == 't850-t2m':
#        req_vars = ['air_temperature', '2m_temperature']
#        for var in req_vars:
#            if var not in params.variables:
#                params.variables.append(var)
#        for plev in [850]:
#            if plev not in params.pressure_levels:
#                params.pressure_levels.append(var)
#                
#    elif params.lts_type == 't850-t1000':
#        req_vars = ['air_temperature']
#        for var in req_vars:
#            if var not in params.variables:
#                params.variables.append(var)
#        for plev in [1000, 850]:
#            if plev not in params.pressure_levels:
#                params.pressure_levels.append(var)
#              
#                
#    elif params.lts_type == 't925-t1000':
#        req_vars = ['air_temperature']
#        for var in req_vars:
#            if var not in params.variables:
#                params.variables.append(var)
#        for plev in [1000, 925]:
#            if plev not in params.pressure_levels:
#                params.pressure_levels.append(var)
#    else:
#        print('Other definitions of LTS type not supported yet')
        
        
    
    for variable in params.variables:
        print('Loading variable ' + variable)
        if variable in ['sea_ice_concentration', 'sea_ice_thickness']:
            subset_latlon=False
        else:
            subset_latlon=True

        ds_dict = util.load_dataset(variable, params, subset_latlon=subset_latlon)
        
        # either here or within load_dataset, remake dataset to have matching names
        
        for ens in ds_dict:
            
            if variable in ['air_temperature', 'eastward_wind', 'northward_wind']:
                for plevel in params.pressure_levels:
                    fn = '_'.join([params.source, ens, variable, str(plevel)])
                    ds = ds_dict[ens].sel(plevel=plevel*100).drop('plevel')
                    ds = ds.sortby(ds.time)
                    ds.to_netcdf(
                params.save_location + fn + '.nc')
            else:        
                fn = '_'.join([params.source, ens, variable])
                ds = ds_dict[ens].sortby(ds_dict[ens].time)
                ds.to_netcdf(
                    params.save_location + fn + '.nc')
        del ds_dict
    return
    

def load_data(variable, params, plevel=None):
    """Imports downloaded data into workspace."""
    ddict = {}
    for ensemble in params.ensembles():
        fname = util.filename_loaded(variable, params, ensemble, plevel)
        with xr.open_dataset(fname) as ds:
            ddict[ensemble] = ds
    return ddict
    
    
    
    
    
