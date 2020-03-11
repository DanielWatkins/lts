"""Helper functions (as you may have guessed). The 
main functions here are for generating lists of
filenames from glade, dealing with different ways of
organizing the model datasets (such as splitting large
files into different years), and navigating differing 
file names.

"""
import numpy as np
import os
import pandas as pd


#### Generating Paths to Model and Reanalysis Data ####


def get_filenames(variable, params):
    """Returns list of file names matching spefications in params

    """
    frequency = params.frequency
    ensemble_members = params.ensemble_members
    source = params.source
    varname = get_varnames(variable, params)
    def get_cesmle_file_names(varname, params):
        """Reads filenames for cesm large ensemble from glade
        and returns as a master list."""

        domain = 'atm' # Can add a step here when other variables become available
        model = 'cam'
        
        if varname == 'ICEFRAC':
            model = 'cice'

        assert freq == 'monthly', 'Frequency must be monthly, for now.'

        file_loc = '/'.join(['/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE',
                             domain, 'proc/tseries', freq, varname])
        if domain=='atm':
            if freq=='daily':
                fcode = 'h1'
            elif freq == 'monthly':
                fcode = 'h0'
        elif domain == 'ice':
            if freq == 'monthly':
                fcode = 'h'
            if freq == 'daily':
                fcode = 'h1'
                file_loc += '_d'
            varname += '_nh'

        ensembles = [str(i).zfill(3) for i in range(1,36)]+[str(i) 
                                           for i in range(101,108)]
        
        if ensemble_members == 'first':
            ensembles = [ensembles[0]]
        elif ensemble_members == 'random':
            ensembles = [np.random.choice(ensembles, 1)]

        
        flist = os.listdir(file_loc)
        files = []
        for ens in ensembles:
            prefix = ['.'.join(['b.e11.B20TRC5CNBDRD.f09_g16', 
                                ens, model, fcode, varname]),
                  '.'.join(['b.e11.BRCP85C5CNBDRD.f09_g16', 
                            ens, model, fcode, varname])]
            files += [f for f in flist if f.startswith(prefix[0])] 
            files += [f for f in flist if f.startswith(prefix[1])]
        flist =  [file_loc + '/' + f for f in files]
        flist.sort()
        
        return flist
    
    def get_cmip6_file_names(varname, params):
        """Reads filenames for cmip6 from glade
        and returns as a master list."""
        
        # To extended to full cmip6, will need reference function here
        if params.source == 'cesm2-cam6':
            model_group = 'NCAR'
            model = 'CESM2'
        elif params.source == 'cesm2-waccm':
            model_group = 'NCAR'
            model = 'CESM2-WACCM'
        else:
            print('Need to adapt to full CMIP6 list')
            

  
        assert params.frequency == 'monthly', 'Frequency must be monthly, for now.'

        if varname in ['siconca', 'sithick']:
            dfreq = 'SImon'
        elif varname in ['sftlf']:
            dfreq = 'fx'
        else:
            dfreq = 'Amon'
            
        file_loc = '/'.join(['/glade/collections/cmip/CMIP6/CMIP',
                             model_group, model, 'historical'])
        
        ensembles = os.listdir(file_loc) 
        if ensemble_members == 'first':
            ensembles = [ensembles[0]]
        elif ensemble_members == 'random':
            ensembles = [np.random.choice(ensembles, 1)]
            
        flist = []
        for ens in ensembles:
            floc = '/'.join([file_loc, ens, dfreq,
                                         varname,'gn','latest'])
            
            files = os.listdir(floc)
            
            flist += [floc + '/' + f for f in files]# if f.endswith('.nc')]
        flist.sort()
        return flist

    def get_cmip5_file_names(varname, params):
        """Reads filenames for cmip5 from glade
        and returns as a master list."""
        
        # To extended to full cmip6, will need reference function here
        model_info = cmip5_models()[params.source]
        model_group = model_info['model_group']
        model = model_info['model_name']

  
        assert params.frequency == 'monthly', 'Frequency must be monthly, for now.'
        freq = 'mon'
        
        if varname == 'sftlf':
            freq = 'fx'
            domain = 'atmos'
            dfreq = 'fx'
        elif varname in ['sic', 'sit']:
            domain = 'seaIce'
            dfreq = 'OImon'
        else:
            domain = 'atmos'
            dfreq = 'Amon'
            
        file_loc = '/'.join(['/glade/collections/cmip/cmip5/output1',
                             model_group, model, 'historical', freq, domain, dfreq])
        
        
        ensembles = os.listdir(file_loc) 
        if ensemble_members == 'first':
            ensembles = [ensembles[0]]
        elif ensemble_members == 'random':
            ensembles = [np.random.choice(ensembles, 1)]
            
        flist = []
        for ens in ensembles:
            floc = '/'.join([file_loc, ens, 'latest', varname])
            
            files = os.listdir(floc)
            
            flist += [floc + '/' + f for f in files]# if f.endswith('.nc')]
        flist.sort()
        return flist   
    
    
    ###################### MAIN PART ##############################
    if source=='cesm-le':
        return get_cesmle_file_names(varname, params)
    
    elif source in cmip5_models():
        return get_cmip5_file_names(varname, params)
        
    elif source in cmip6_models():
        return get_cmip6_file_names(varname, params)

    elif source=='era-i':
        return get_erai_file_names(varname, params)
    
    elif source=='era5':
        return get_era5_file_names(varname, params)
    

# def get_cesm2_file_names(ens, variable, freq, domain):

# def get_erai_file_names(variable, freq, domain):

# def get_era5_file_names(variable, freq, domain):

def get_varnames(variable, params):
    """Given a variable name and parameters, provides the name of 
    the corresponding variable in the source data if it exists.
    """
    
    cesmle = {'latitude': 'LAT',
              'longitude': 'LON',
              'land_mask': 'LANDFRAC',
              'sea_ice_concentration': 'ICEFRAC',
              'sea_ice_thickness': 'SIT',
              'total_cloud_cover': 'CLDTOT',
              'low_cloud_fraction': 'CLDLOW',
              '2m_temperature': 'TREFHT',
              'air_temperature': 'T',
              'sensible_heat_flux': 'SHFLX',
              'latent_heat_flux': 'LHFLX',
              '10m_wind_speed': 'U10',
              'sea_level_pressure': 'PSL',
              'surface_downward_longwave': np.nan,
              'surface_upward_longwave': np.nan,
              'surface_downward_shortwave': np.nan,
              'surface_upward_shortwave': np.nan,
              'net_surface_longwave': 'FLNS',
              'net_surface_shortwave': 'FSNS',
              'liquid_water_path': 'TGCLDLWP', 
              'total_water_path': 'ICLDTWP',
              'snow_depth_water_equivalent': 'SNOWHICE',
              'ice area snapshot': 'aisnap',
              'grid cell mean ice thickness': 'hi',
              'grid cell mean snow thickness': 'hs'
             }
              

    cmip6 =  {'10m_wind_speed': 'sfcWind',
              '2m_temperature': 'tas',
              'air_temperature': 'ta',
              'condensed_ice_path': 'clivi',
              'condensed_water_path': 'clwvi',
              'eastward_wind': 'ua',
              'ice_water_path': 'wlivi',
              'land_area_fraction': 'sftlf',
              'latent_heat_flux': 'hfss',
              'latitude': 'lat',
              'longitude': 'lon',
              'low_cloud_fraction': np.nan,
              'mass_fraction_cloud_ice': 'cli',
              'mass_fraction_cloud_liquid': 'clw',
              'net_surface_longwave': np.nan,
              'net_surface_shortwave': np.nan,
              'northward_wind': 'va',
              'precipitable_water': 'pwr',
              'pressure_level': 'plev',
              'relative_humidity': 'hur',
              'sea_ice_concentration': 'siconca',
              'sea_ice_thickness': 'sithick',
              'sea_level_pressure': 'psl',
              'sensible_heat_flux': 'hfls',
              'specific_humidity': 'hus',
              'surface_downward_longwave': 'rlds',
              'surface_downward_shortwave': 'rsds',
              'surface_pressure': 'ps',
              'surface_relative_humidity': 'hurs',
              'surface_specific_humidity': 'huss',
              'surface_upward_longwave': 'rlus',
              'surface_upward_shortwave': 'rsus',
              'toa_longwave_all_sky': 'rlut',
              'toa_longwave_clear_sky': 'rlutcs',
              'total_cloud_cover': 'clt',
              'vertical_velocity': 'wap'}
    
    cmip5 =  {'10m_wind_speed': 'sfcWind',
              '2m_temperature': 'tas',
              'air_temperature': 'ta',
              'condensed_ice_path': 'clivi',
              'condensed_water_path': 'clwvi',
              'eastward_wind': 'ua',
              'ice_water_path': 'wlivi',
              'land_area_fraction': 'sftlf',
              'latent_heat_flux': 'hfss',
              'latitude': 'lat',
              'longitude': 'lon',
              'low_cloud_fraction': np.nan,
              'mass_fraction_cloud_ice': 'cli',
              'mass_fraction_cloud_liquid': 'clw',
              'net_surface_longwave': np.nan,
              'net_surface_shortwave': np.nan,
              'northward_wind': 'va',
              'precipitable_water': 'pwr',
              'pressure_level': 'plev',
              'relative_humidity': 'hur',
              'sea_ice_concentration': 'sic',
              'sea_ice_thickness': 'sit',
              'sea_level_pressure': 'psl',
              'sensible_heat_flux': 'hfls',
              'specific_humidity': 'hus',
              'surface_downward_longwave': 'rlds',
              'surface_downward_shortwave': 'rsds',
              'surface_pressure': 'ps',
              'surface_relative_humidity': 'hurs',
              'surface_specific_humidity': 'huss',
              'surface_upward_longwave': 'rlus',
              'surface_upward_shortwave': 'rsus',
              'toa_longwave_all_sky': 'rlut',
              'toa_longwave_clear_sky': 'rlutcs',
              'total_cloud_cover': 'clt',
              'vertical_velocity': 'wap'}
              
    erai =   {'latitude':np.nan,
              'longitude':np.nan,
              'land_mask': np.nan,
              'sea_ice_concentration': 'CI_GDS4_SFC_S123',
              'total_cloud_cover': 'TCC_GDS4_SFC_S123',
              'low_cloud_fraction': 'LCC_GDS4_SFC_S123',
              '2m_temperature': '2T_GDS4_SFC_S123',
              'air_temperature': 'T_GDS4_ISBL_S123',
              'sensible_heat_flux': 'SSHF_GDS4_SFC_120',
              'latent_heat_flux': 'SLHF_GDS4_SFC_120',
              '10m_wind_speed': '10SI_GDS4_SFC_S123',
              'sea_level_pressure': 'MSL_GDS4_SFC_S123',
              'surface_downward_longwave':np.nan,
              'surface_upward_longwave':np.nan,
              'surface_downward_shortwave':np.nan,
              'surface_upward_shortwave':np.nan,
              'net_surface_longwave': 'STR_GDS4_SFC_120',
              'net_surface_shortwave': 'SSR_GDS4_SFC_120'
             }   
              
    era5 =   {'latitude':np.nan,
              'longitude':np.nan,
              'land_mask': np.nan,
              'sea_ice_concentration':np.nan,
              'total_cloud_cover':np.nan,
              'low_cloud_fraction':np.nan,
              '2m_temperature':np.nan,
              'sensible_heat_flux':np.nan,
              'latent_heat_flux':np.nan,
              '10m_wind_speed':np.nan,
              'sea_level_pressure':np.nan,
              'surface_downward_longwave':np.nan,
              'surface_upward_longwave':np.nan,
              'surface_downward_shortwave':np.nan,
              'surface_upward_shortwave':np.nan,
              'net_surface_longwave':np.nan,
              'net_surface_shortwave':np.nan}     
              
    sources = {'cesm-le': cesmle,
               'cesm2-cam6': cmip6,
               'cesm2-waccm': cmip6,
               'cesm1-cam5': cmip5,
               'cesm1-waccm': cmip5,
               'era-i': erai,
               'era5': era5}
    
    source = params.source
    if source not in sources.keys():
        print('Source must be one of the following: \n----------------------- \n' + '\n'.join(list(sources.keys())))
    elif variable not in sources[source].keys():
        print('Variable must be one of the following: \n----------------------- \n'  + '\n'.join(list(sources[source].keys())))
    else:
        return sources[source][variable]

    
    
def cmip5_models():
    """Returns a dictionary with shortnames as keys, 
    and entries for the model group and model name."""
    
    return {'cesm1-cam5': {'model_group': 'NSF-DOE-NCAR', 'model_name': 'CESM1-CAM5'},
            'cesm1-waccm': {'model_group': 'NSF-DOE-NCAR', 'model_name': 'CESM1-WACCM'}}

def cmip6_models():
    """Returns a dictionary with shortnames as keys, 
    and entries for the model group and model name."""
    
    return {'cesm2-cam6': {'model_group': 'NSF-DOE-NCAR', 'model_name': 'CESM1-CAM5'},
            'cesm2-waccm': {'model_group': 'NSF-DOE-NCAR', 'model_name': 'CESM1-WACCM'}}

    
def weights(lats, lons, area=True):
    """Compute area for each lat/lon grid box. Assumes that
    the polar grid points are offset (so values are centered halfway
    between 90 and the next highest lat), and the rest of the grid points
    represent grid centers. (Not sure how this plays with the interp to 1deg.)
    
    
    
    """
    import numpy as np
    dlat = np.diff(lats)[0]/2
    nlon = len(lons)
    lats = lats.copy()
    lats[0] = lats[0] + dlat
    lats[-1] = lats[-1] - dlat
    R = 6356 # earth radius in km
    sin_lat1 = np.sin((lats + dlat) * np.pi / 180.)
    sin_lat2 = np.sin((lats - dlat) * np.pi / 180.)
    sin_lat = np.abs(sin_lat1 - sin_lat2)
    weights = (np.zeros((len(lats), nlon)) + np.array([sin_lat / (np.sum(sin_lat)*nlon)]).T)
    # Not sure what the np.sum(sin_lat) is for
    if area:
        grid_area = R**2 * weights
        return grid_area
    else:
        return weights
    
def make_tempfile(variable, params):
    """make a temporary file for regridder to read. Requires
    2m_temperature to be used as the ref file. Currently expects
    the lat/lon to be 'lat' and 'lon' respectively. Keeps the file
    name except adds _interpolated.nc at the end."""
    
    
    # TODO: add save_name to the thing, so that the naming convention is followed.
    
    varname = get_varnames(variable, params)
    fnames = params.filenames()
    reffile = fnames['2m_temperature'][0]
    targfiles = fnames[variable]
    saveloc = params.save_location
    min_lat = params.minimum_latitude
    begin_time = params.begin_time
    end_time = params.end_time
    source = params.source
    info = [source, variable, saveloc, reffile, varname, min_lat, 
            begin_time, end_time] + targfiles
    names = ['source', 'variable', 'save_location','reference', 
             'varname', 'minimum_latitude',
            'begin_time', 'end_time'] + ['target' 
                               for i in range(len(targfiles))]
    
    df = pd.DataFrame({'name':names, 'info':info})
    df.to_csv('data_for_gridding.tmp', index=False)
    
def load_dataset(variable, params, subset_time=True, subset_latlon=True):
    """Loads dataset for <variable> based on time and lat/lon subsets from params.
    Currently can only subset latlon if the grid is rectilinear. Returns dictionary with 
    a key for each ensemble member.
    
    Todo: To avoid duplicating coordinates, a new dataset is created with truncated lat lon values
    (4 decimal places). 
    """
    import warnings 
    import xarray as xr
    import pandas as pd
    
    def ensemble_finder(fname):
        """Pulls out the ensemble from the file name"""
        components = fname.split('/')[-1].split('_')
        for comp in components:
            if comp[0] == 'r':
                if comp[1].isdigit():
                    return comp
        
    def check_dates(fname, params):
        """Checks dates, assuming that the time range is given right before the .nc.
        Only takes year into account"""
        daterng = fname.split('_')[-1].split('.')[0].split('-')
        begin = int(daterng[0][0:4])
        end = int(daterng[1][0:4])
        begin_p = int(params.begin_time[0:4])
        end_p = int(params.end_time[0:4])

        if (begin <= begin_p ) & (end >= begin_p):

            return True
        elif (begin >= begin_p) & (begin <= end_p):
            return True
        else:
            return False
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fnames = params.filenames()[variable]
        ensembles = [ensemble_finder(f) for f in fnames]
        fn_dict = {}
        if params.ensemble_members == 'first':
            ens = ensembles[0]
#             if len(ensembles) > 1:
#                 ens = ensembles[0]
#             else:
#                 ens = ensembles[0]
            fnames = [f for f in fnames if (ensemble_finder(f)==ens) and check_dates(f, params)]
            fn_dict[ens] = fnames
        elif params.ensemble_members == 'all':
            for ens in ensembles:
                fnames_ens = [f for f in fnames if (ensemble_finder(f)==ens) and check_dates(f, params)]
                fn_dict[ens] = fnames_ens
        # later: random option?
        else:
            print('Bad selection of ensembles. Must be either first or all.')
            return
        # todo:
        # group fnames by ensemble
        # loop through each group
        #return(fnames, ensembles)

        ds_dict = {}
        for ens in fn_dict:
            ds_list = []
            for fname in fn_dict[ens]:
                dsx = clean_dataset(xr.open_dataset(fname), variable, params)
                
                if subset_time and subset_latlon:
                    ds = dsx.sel(
                            time=slice(params.begin_time, params.end_time),
                            lat=slice(params.minimum_latitude, params.maximum_latitude),
                            lon=slice(params.minimum_longitude, params.maximum_longitude))
                    ds.load()
                elif subset_time:
                    ds = dsx.sel(
                            time=slice(params.begin_time, params.end_time))
                    ds.load()
                elif subset_latlon:
                    ds = dsx.sel(
                            lat=slice(params.minimum_latitude, params.maximum_latitude),
                            lon=slice(params.minimum_longitude, params.maximum_longitude))
                    ds.load()
                ds_list.append(ds)
                dsx.close()
                if len(ds.time) == len(pd.date_range(params.begin_time, params.end_time, freq='1MS')):
                    break
                    del ds
                del ds
            ds_dict[ens] = xr.concat(ds_list, dim='time')
        return ds_dict
    

def clean_dataset(ds, variable, params):
    """Ensures that the dataset has time series, lat/lon values, etc that
    will play well with others: i.e., rounding lat lon so that no doubling 
    happens, making time series follow a calendar, etc."""
    import xarray as xr
    import numpy as np
    
        
    data = ds[get_varnames(variable, params)].values
    lats = np.round(ds[get_varnames('latitude', params)].values, 6)
    lons = np.round(ds[get_varnames('longitude', params)].values, 6)
    times = pd.date_range(params.begin_time, freq='MS', periods=len(ds['time']))
    
    if variable in ['air_temperature', 'eastward_wind', 'northward_wind']:
        plevels = ds.coords[get_varnames('pressure_level', params)].values
        return xr.Dataset({variable: (('time', 'plevel', 'lat', 'lon'), data)}, 
                   coords={'time': times, 'plevel': plevels, 
                           'lat': lats, 'lon': lons})

    return xr.Dataset({variable: (('time', 'lat', 'lon'), data)}, 
               coords={'time': times, 'lat': lats, 'lon': lons})

def filename_loaded(variable, params, ensemble, plevel=None):
    if plevel is not None:
        return params.save_location + '_'.join([params.source, ensemble, variable, str(plevel)]) + '.nc'
    else:
        return params.save_location + '_'.join([params.source, ensemble, variable]) + '.nc'

def compute_lts(params):
    """Make a new dataset for lts based on the lts definition."""
    import xarray as xr
    for ensemble in params.ensembles():
        if params.lts_type == 't850-t2m':
            ds_top = xr.open_dataset(filename_loaded('air_temperature', 
                                                     params, ensemble, 850))
            
            ds_bottom = xr.open_dataset(filename_loaded('2m_temperature', params, ensemble))
        ds_top = ds_top.sortby(ds_top.time)
        ds_bottom = ds_bottom.sortby(ds_bottom.time)
        lts = ds_top['air_temperature'].values - ds_bottom['2m_temperature'].values
        ds_lts = xr.Dataset({'lts': (('time', 'lat', 'lon'), lts.squeeze())}, 
                   coords={'time': ds_bottom.time,
                           'lat': ds_bottom['lat'].values,
                           'lon': ds_bottom['lon'].values})
        ds_top.close()
        ds_bottom.close()
        ds_lts.to_netcdf(filename_loaded('lts', params, ensemble)) # Add lts_type option here?


    
    
    
    
    
    
    
