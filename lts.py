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
    source: One of 'cesm-le' (default), 'cesm2-cam6', 'cesm2-waccm', 'era-i', 'era5'.
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
        print('source: ', self.source)
        print('lts_type: ', self.lts_type)
        print('save_location: ', self.save_location)
        print('begin_time: ', self.begin_time)
        print('end_time: ', self.end_time)
        print('variables: ', self.variables)
        print('pressure_levels: ', self.pressure_levels)
        print('ensemble_members: ', self.ensemble_members)
        print('frequency: ', self.frequency)
        print('maximum_latitude: ', self.minimum_latitude)
        print('minimum_latitude: ', self.maximum_latitude)
        print('maximum_longitude: ', self.maximum_longitude)
        print('minimum_longitude: ', self.minimum_longitude)
        
    def filenames(self):
        fn_dict = {}
        for varname in self.variables:
            fn_dict[varname] = util.get_filenames(
                varname, self)
        return fn_dict


def get_data(params):
    """Wrapper for load_dataset. Checks to see if required variables
    are included based on the desired lts_type.
    """
    
    filenames = params.filenames()
    
    if params.lts_type == 't850-t2m':
        req_vars = ['air_temperature', '2m_temperature']
        for var in req_vars:
            if var not in params.variables:
                params.variables.append(var)
        for plev in [850]:
            if plev not in params.pressure_levels:
                params.pressure_levels.append(var)
                
    elif params.lts_type == 't850-t1000':
        req_vars = ['air_temperature']
        for var in req_vars:
            if var not in params.variables:
                params.variables.append(var)
        for plev in [1000, 850]:
            if plev not in params.pressure_levels:
                params.pressure_levels.append(var)
              
                
    elif params.lts_type == 't925-t1000':
        req_vars = ['air_temperature']
        for var in req_vars:
            if var not in params.variables:
                params.variables.append(var)
        for plev in [1000, 925]:
            if plev not in params.pressure_levels:
                params.pressure_levels.append(var)
    else:
        print('Other definitions of LTS type not supported yet')
        
        
    
    for variable in params.variables:
        print('Loading variable ' + variable)
        if variable in ['sea_ice_concentration', 'sea_ice_thickness']:
            subset_latlon=False
        else:
            subset_latlon=True

        ds_dict = util.load_dataset(variable, params, subset_latlon=subset_latlon)
        
        # either here or within load_dataset, remake dataset to have matching names
        
        for ens in ds_dict:
            ds_dict[ens].to_netcdf(
                params.save_location + params.source + '_' + ens + '_' + variable + '.nc')
        del ds_dict
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



    
    
    
    
    
    
    
    
    
