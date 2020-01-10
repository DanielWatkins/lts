"""A collection of functions for analyzing lower tropospheric stability.



(c) Daniel Watkins 2020
"""

import numpy as np
import xarray as xr
import pandas as pd
import lts.utilities.varnames as ltsv
import lts.utilities.filenames as ltsf
# import lts.utilities.util as ltsu

class parameters(source):
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
    """

    def __init__(self, source):
        self.source = source
        self.lts_type = 't850-t2m'
        self.save_location = '~/data/'
        self.begin_time = '1979-01-01'
        self.end_time = '2018-12-31'
        self.variables = ['lts', 't850', 't2m']
        self.ensemble_members = 'all'
        self.frequency = 'monthly'
        self.minimum_latitude = 60
        self.maximum_latitude = 90
        self.minimum_longitude = 0
        self.maximum_longitude = 360
        
        

def get_data(params):
    """Based on information in params object, pull data
    from glade/collections, subset it, and interpolate to
    a pressure level if needed."""