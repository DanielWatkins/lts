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
        print('Source: ', self.source)
        print('LTS Type: ', self.lts_type)
        print('Save location: ', self.save_location)
        print('Begin time: ', self.begin_time)
        print('End time: ', self.end_time)
        print('Variables: ', self.variables)
        print('Ensemble members: ', self.ensemble_members)
        print('Frequency: ', self.frequency)
        print('Max latitude: ', self.minimum_latitude)
        print('Min latitude: ', self.maximum_latitude)
        print('Max longitude: ', self.maximum_longitude)
        print('Min longitude: ', self.minimum_longitude)
        
    def filenames(self):
        # Update to allow only the first ensemble member to be called
        fn_dict = {}
        for varname in self.variables:
            fn_dict[varname] = util.get_filenames(
                varname,
                self.source,
                self.frequency,
                self.ensemble_members)
        return fn_dict

def get_data(params):
    """Based on information in params object, pull data
    from glade/collections, subset it, and interpolate to
    a pressure level if needed."""
    
    filenames = params.filenames()
    
    if params.lts_type == 't850-t2m':
        req_vars = ['air_temperature', '2m_temperature',
                    'sea_ice_concentration']
        for var in req_vars:
            if var not in params.variables:
                params.variables += var
    elif params.lts_type == 't925-t1000':
        req_vars = ['air_temperature',
                    'sea_ice_concentration']
        for var in req_vars:
            if var not in params.variables:
                params.variables += var            
    else:
        print('Other definitions of LTS type not supported yet')
        
    
        
        
    # Handling ensembles:
    # Handling data split into multiple files:
        
    
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



    
    
    
    
    
    
    
    
    
