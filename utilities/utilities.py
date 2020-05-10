"""Helper functions (as you may have guessed). The 
main functions here are for generating lists of
filenames from glade, dealing with different ways of
organizing the model datasets (such as splitting large
files into different years), and navigating differing 
file names.

I'd like to break this apart into a few files so that 
the file is more single-purpose.
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
    if ensemble_members == 'random':
        ensemble_size = params.ensemble_size
        
    source = params.source
    varname = get_varnames(variable, params)
    try:
        np.isnan(varname)
        print('Undefined variable ' + variable + 'for ' + source)
        return
    except:
        pass
    
    def get_cesmle_file_names(varname, params):
        """Reads filenames for cesm large ensemble from glade
        and returns as a master list."""

        domain = 'atm' # Can add a step here when other variables become available
        model = 'cam'
        
#         if varname == 'ICEFRAC':
#             model = 'cice'
        freq = frequency
        freq = 'monthly' # Don't have daily data yet
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
            ensembles = np.random.choice(ensembles, ensemble_size)
            if ensemble_size == 1:
                ensembles = list(ensembles)
        elif ensemble_members == 'all':
            pass
        else:
            print('I don\'t know what you mean, using all ensemble members.')
            
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
            if ensemble_size == 1:
                ensembles = list(np.random.choice(ensembles, ensemble_size))
            if ensemble_size > len(ensembles):
                print('Ensemble size greater than number of ' +
                      'available ensemble members, using all ensemble members')
                pass
            else:
                ensembles = np.random.choice(ensembles, ensemble_size)
        elif ensemble_members == 'all':
            pass
        else:
            print('I don\'t know what you mean, using all ensemble members instead.')

                 
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
            if ensemble_size == 1:
                ensembles = list(np.random.choice(ensembles, ensemble_size))
            if ensemble_size > len(ensembles):
                print('Ensemble size greater than number of ' +
                      'available ensemble members, using all ensemble members')
                pass
            else:
                ensembles = np.random.choice(ensembles, ensemble_size)
        elif ensemble_members == 'all':
            pass
        else:
            print('I don\'t know what you mean, using all ensemble members instead.')


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
    
    cesmle = {'latitude': 'lat',
              'longitude': 'lon',
              'land_area_fraction': 'LANDFRAC',
              'sea_ice_concentration': 'ICEFRAC',
              'sea_ice_thickness': 'SIT',
              'total_cloud_cover': 'CLDTOT',
              'low_cloud_fraction': 'CLDLOW',
              '2m_temperature': 'TREFHT',
              'air_temperature': 'T',
              'pressure_level': 'lev',
              'sensible_heat_flux': 'SHFLX',
              'latent_heat_flux': 'LHFLX',
              '10m_wind_speed': 'U10',
              'sea_level_pressure': 'PSL',
              'surface_pressure': 'PS',
              'surface_downward_longwave': 'FLDS',
              'surface_upward_longwave': np.nan,
              'surface_downward_shortwave': 'FSDS',
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
              'latent_heat_flux': 'hfls',
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
              'sensible_heat_flux': 'hfss',
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
              'latent_heat_flux': 'hfls',
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
              'sensible_heat_flux': 'hfss',
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
    
    
def ensemble_finder(fname, params):
    """Pulls out the ensemble from the file name"""
    def cmip_ensemble_finder(fname):
        """Pulls out the ensemble from the file name"""
        components = fname.split('/')[-1].split('_')
        for comp in components:
            if comp[0] == 'r':
                if comp[1].isdigit():
                    return comp

    def cesmle_ensemble_finder(fname):
        """Pulls out the ensemble from the file name"""
        ensemble = fname.split('/')[-1].split('.')[4]
        return ensemble
    if params.source == 'cesm-le':
        return cesmle_ensemble_finder(fname)
    else:
        return cmip_ensemble_finder(fname)

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
    
    def check_dates(fname, params):
        """Checks dates, assuming that the time range is given right before the .nc.
        Only takes year into account"""
        
        if params.source == 'cesm-le':
            daterng = fname.split('_')[-1].split('.')[-2].split('-')
        else:
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
        ensembles = [ensemble_finder(f, params) for f in fnames]
        fn_dict = {}
        if params.ensemble_members == 'first':
            ens = ensembles[0]
            fnames = [f for f in fnames if (ensemble_finder(f, params)==ens) and check_dates(f, params)]
            fn_dict[ens] = fnames
        elif params.ensemble_members == 'all':
            for ens in ensembles:
                fnames_ens = [f for f in fnames if (ensemble_finder(f, params)==ens) and check_dates(f, params)]
                fn_dict[ens] = fnames_ens
        elif params.ensemble_members == 'random':
            for ens in ensembles:
                fnames_ens = [f for f in fnames if (ensemble_finder(f, params)==ens) and check_dates(f, params)]
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
            print('Loading ensemble member ' + ens)
            ds_list = []
            for fname in fn_dict[ens]:
                if not ((params.source == 'cesm-le') & (variable == 'air_temperature')): 
                    # add other 3d variables somehow
                    # For the cmip5 and cmip6 data, we subset and "clean" the dataset
                    print('cleaning ' + fname)
                    dsx = clean_dataset(xr.open_dataset(fname), variable, params)

                    if subset_time and subset_latlon:
                        ds = dsx.sel(
                                time=slice(params.begin_time, params.end_time),
                                lat=slice(params.minimum_latitude, 
                                          params.maximum_latitude),
                                lon=slice(params.minimum_longitude, 
                                          params.maximum_longitude))
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
                else:
                    # Only does the 850 layer right now
                    
                    varname = get_varnames(variable, params)
                    suffix = '.'.join(fname.split('.')[-2:])
                    fname_p = '/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/PS/b.e11.B20TRC5CNBDRD.f09_g16.' + ens + '.cam.h0.PS.' + suffix
                    ds = vinth2p_cesm(in_file_var=fname,
                                      in_file_pvar=fname_p,
                                      min_lat=params.minimum_latitude,
                                      max_lat=params.maximum_latitude,
                                      p_levels=[850],
                                      time_range=(params.begin_time, params.end_time),
                                      newvar_name=variable)
                    ds_list.append(ds)
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
    """Convenience function to grab the filename of the cleaned and moved dataset."""
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


def get_landmask(params):
    """Loads and saves the land fraction dataset."""
    import xarray as xr
    import numpy as np
    
    fname = get_filenames('land_area_fraction', params)[0]
    with xr.open_dataset(fname) as lf:
        ds_land = xr.Dataset({'land_area_fraction': (('lat', 'lon'), 
                             lf[get_varnames('land_area_fraction', params)])}, 
               coords={'lat': np.round(lf['lat'].values, 6),
                       'lon': np.round(lf['lon'].values, 6)})

        ds_land.sel(lat=slice(params.minimum_latitude, params.maximum_latitude),
                    lon=slice(params.minimum_longitude,
                             params.maximum_longitude)).to_netcdf(
                             filename_loaded('land_area_fraction', params, '')) # Add lts_type 
    

# Weight the months for a representatitive histogram.
def get_month_weights(dataset, timevar, altseasons=False):
    """Mostly copied directly from xarray's help docs. 
    If altseasons, then group November in with winter."""
    import xarray as xr
    import numpy as np
    
    dpm = {'noleap': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
           '365_day': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
           'standard': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
           'gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
           'proleptic_gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
           'all_leap': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
           '366_day': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
           '360_day': [0, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]}

    def leap_year(year, calendar='standard'):
        """Determine if year is a leap year"""
        leap = False
        if ((calendar in ['standard', 'gregorian',
            'proleptic_gregorian', 'julian']) and
            (year % 4 == 0)):
            leap = True
            if ((calendar == 'proleptic_gregorian') and
                (year % 100 == 0) and
                (year % 400 != 0)):
                leap = False
            elif ((calendar in ['standard', 'gregorian']) and
                     (year % 100 == 0) and (year % 400 != 0) and
                     (year < 1583)):
                leap = False
        return leap

    def get_dpm(time, calendar='standard'):
        """
        return a array of days per month corresponding to the months provided in `months`
        """
        month_length = np.zeros(len(time), dtype=np.int)

        cal_days = dpm[calendar]

        for i, (month, year) in enumerate(zip(time.month, time.year)):
            month_length[i] = cal_days[month]
            if leap_year(year, calendar=calendar):
                month_length[i] += 1
        return month_length

    month_length = xr.DataArray(get_dpm(dataset[timevar].to_index(), calendar='noleap'),
                                coords=[dataset[timevar]], name='month_length')

    if altseasons:
        altseasons_list = []
        for m in month_length.time.dt.month:
            if m in [11, 12, 1, 2]:
                altseasons_list.append('NDJF')
            elif m in [3, 4, 5]:
                altseasons_list.append('MAM')
            elif m in [6, 7, 8]:
                altseasons_list.append('JJA')
            else:
                altseasons_list.append('SO')

        test = xr.DataArray(altseasons_list, coords=[month_length.time], 
                            dims={'time':month_length.time.time})
        test.name = 'season'
        month_length = xr.concat(month_length, test)
        weights = month_length.groupby('season') / month_length.groupby('season').sum()
        return weights
    
    else:
        weights = month_length.groupby(timevar + '.season') / month_length.groupby(timevar + '.season').sum()
        return weights
    
def compute_histogram(data, bins, lat, lon, surface, land, ice):
    """Compute weighted histogram values. Data is an array and has 
    dimensions lat, lon. surface is a string: land, ocean_ice,
    ocean, ice. land is the array with landfrac, and ice is an array with
    the same dimensions as data."""
    import numpy as np
    
    if surface == 'land':
        fraction = land > 0.99
        fraction = np.broadcast_to(fraction, data.shape)
    elif surface == 'ocean_ice':
        fraction = (1 - land) > 0.99
        fraction = np.broadcast_to(fraction, data.shape)
    elif surface == 'ocean':
        fraction = (1 - (ice + land)) > 0.99
    elif surface == 'ice':
        fraction = ice > 0.15
    elif surface == 'high_sic':
        fraction = ice > 0.85
    elif surface == 'low_sic':
        fraction = (ice > 0.15) & (ice <= 0.85)
    elif surface == 'north70':
        fraction = np.zeros(data.size)
        fraction[lat >= 70, :] = True
    else:
        print('Bad surface')
        return
    
    w = np.broadcast_to(weights(lat, lon, np.min(lat)), data.shape)
    w = w * fraction
    w = w/w.sum()
    x, y = np.histogram(np.ravel(data), weights=np.ravel(w), density=True, bins=bins)
    return x


def get_histograms(params_list, variable='lts', plevel=None, period='seasons', minimum_latitude=65, 
                   surfaces=['ice', 'ocean_ice'], bins = np.arange(-20,30,1)):
    """Computes the histograms for a variable over each ice surface for the chosen period.
    Ice surface options include ice, ocean_ice (tested), or land, north70 (untested). Period can be
    seasons, altseasons (november grouped with winter), or months. Bins default is set for lts."""
    import warnings
    import pandas as pd
    import xarray as xr
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") 
        # This is so we don't have to see the "mean of empty slice" error a million times.
        
        bin_centers = (bins[:-1] + bins[1:])/2
        
        if period == 'seasons':
            periods = ['DJF', 'MAM', 'JJA', 'SON']
            altseasons = False
        elif period == 'altseasons':
            periods = ['NDJF', 'MAM', 'JJA', 'SO']
            altseasons = True
        else:
            periods = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
            numbers = np.arange(1, 13)
            altseasons = False
        results = {sfc: {timestep: {} for timestep in periods} for sfc in surfaces}

        for params in params_list:
            print(params.source)
            temp_results = {sfc: {sn: {} for sn in periods} for sfc in surfaces}
            ensembles = params.ensembles()
            for surface in surfaces:
                df = pd.DataFrame(columns=bin_centers,
                                  index=ensembles)
                for timestep in periods:
                    for ens in ensembles:
                        with xr.open_dataset(filename_loaded(variable, params, ens, plevel)) as ds_lts:
                            with xr.open_dataset(filename_loaded('sea_ice_concentration', params, ens)) as ds_sic:
                                ds_sic = ds_sic.sel(lat=slice(minimum_latitude, params.maximum_latitude))
                                if not ('sea_ice_concentration' in ds_sic):
                                    ds_sic.rename({'sic': 'sea_ice_concentration'}, inplace=True)
                                # TODO: fix the cesm1 sea ice thing so the name is right
                                ds_land = xr.open_dataset(filename_loaded('land_area_fraction', params, ''))
                                ds_lts = ds_lts.sel(lat=slice(minimum_latitude, params.maximum_latitude))
                                ds_land = ds_land.sel(lat=slice(minimum_latitude, params.maximum_latitude))

                                mw = get_month_weights(ds_lts, 'time', altseasons=altseasons)
                                if period == 'months':
                                    #return numbers, timestep, mw.time.dt.month, periods
                                    time_selection = mw.time.dt.month == numbers[np.array(periods) == timestep]
                                else:
                                    time_selection = mw.season == timestep
                                df.loc[ens,:] = compute_histogram(
                                        data=ds_lts[variable].loc[time_selection, :, :].values,
                                        bins=bins,
                                        lat=ds_lts['lat'].values,
                                        lon=ds_lts['lon'].values,
                                        surface=surface,
                                        land=ds_land['land_area_fraction'],
                                        ice=ds_sic['sea_ice_concentration'].loc[time_selection, :, :].values)

                    temp_results[surface][timestep][params.source] = df.copy()
                    df_cat = pd.concat(temp_results[surface][timestep], axis=0)
                    if len(results[surface][timestep]) != 0:
                        results[surface][timestep] = pd.concat([results[surface][timestep], df_cat.copy()], axis=0)
                    else:
                        results[surface][timestep] = df_cat
        return results

def vinth2p_cesm(in_file_var, in_file_pvar, **kwargs):
    """Interpolate CESM data from hybrid coordinates to standard
    pressure coordinates. This code is a thin wrapper for the vinth2p function
    from Ngl, which itself is based on the vinth2p function in NCL.
    
    in_file_var is the filename of the variable to be interpolated.
    in_file_press is the filename of the pressure data (either PS or PSL).
    out_file is the filename where the data should be written.
    
    Optional inputs:
    min_lat: minimum latitude to include in the output 
    max_lat: maximum latitude to include in the output
    p_levels: levels to interpolate to
    time_range: date range in form (string date, string date)
    newvar_name: name for variable
    """

    import netCDF4 as nc
    import Ngl
    import numpy as np
    import xarray as xr
    
    ds_var = xr.open_dataset(in_file_var)
    ds_p = xr.open_dataset(in_file_pvar)
    # Get the names of the variables
    pvar = in_file_pvar.split('.')[-3]
    var = in_file_var.split('.')[-3]
    
    # Handle optional arguments
    if 'min_lat' in kwargs:
        min_lat = kwargs['min_lat']
    else:
        min_lat = -90.
        
    if 'max_lat' in kwargs:
        max_lat = kwargs['max_lat']
    else:
        max_lat = 90.0
        
    if 'min_lon' in kwargs:
        min_lon = kwargs['min_lon']
    else:
        min_lon = 0.
        
    if 'max_lon' in kwargs:
        max_lon = kwargs['max_lon']
    else:
        max_lon = 360.0
        
    if 'time_range' in kwargs:
        begin_time, end_time = kwargs['time_range']
    else:
        begin_time = ds_var.time[0].values
        end_time = ds_var.time[-1].values
    
    if 'p_levels' in kwargs:
        p_levels = kwargs['p_levels']
    else:
        # Use WMO mandatory pressure levels
        p_levels = np.array([1000, 925, 850, 700, 500, 400, 
                             300, 250, 200, 150, 100, 70, 
                             50, 30, 20, 10, 7, 5, 3, 2, 1])
    
    
    var_interp = Ngl.vinth2p(
                    datai=ds_var.sel(
                        time=slice(begin_time, end_time),
                        lat=slice(min_lat, max_lat),
                        lon=slice(min_lon, max_lon)
                        ).variables[var].values, 
                    hbcofa=ds_var.variables['hyam'].values, 
                    hbcofb=ds_var.variables['hybm'].values, 
                    plevo=p_levels, 
                    psfc=ds_p.sel(
                        time=slice(begin_time, end_time),
                        lat=slice(min_lat, max_lat),
                        lon=slice(min_lon, max_lon)
                        ).variables[pvar].values, 
                    intyp=2, 
                    p0=ds_var.variables['P0'].values/100., 
                    ii=1, 
                    kxtrp=True)
    
    ds = xr.Dataset({var: (('time', 'lev', 'lat', 'lon'), var_interp)},
                    coords={'time': ds_var.sel(
                        time=slice(begin_time, end_time)).variables['time'].values,
                           'lat': ds_var.lat.values[ds_var.lat.values >= min_lat],
                           'lon': ds_var.lon.values,
                           'lev': p_levels})
    return ds
#    ds.to_netcdf(out_file)
