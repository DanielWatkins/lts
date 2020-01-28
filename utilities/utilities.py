"""Helper functions to access files via the CESM variable names.
get_filenames(variable, source, freq): return a full list of file names 
     that match the specified variable, source, and frequency.



"""
import numpy as np
import os
import pandas as pd

def get_filenames(variable, source, freq='monthly', **kwargs):
    """Returns list of file names matching the ensemble
    number, variable, and frequency.
    
    variable=variable name in the source data set
    freq=daily or monthly or 6-hourly
    domain=atmosphere or ice, probably not used for ERA
    source=ERAI, ERA5, CESMLE, CESM2-CAM6, CESM2-WACCM
    ens=string eg '001'
    """
    
    def get_cesmle_file_names(variable, freq, **kwargs):
        """Reads filenames for cesm large ensemble from glade
        and returns as a master list."""

        domain = 'atm' # Can add a step here when other variables become available
        model = 'cam'
        
        if variable == 'ICEFRAC':
            model = 'cice'
            
        
        
        
        assert freq == 'monthly', 'Frequency must be monthly, for now.'

        file_loc = '/'.join(['/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE',
                             domain, 'proc/tseries', freq, variable])
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
            variable += '_nh'

        ensembles = [str(i).zfill(3) for i in range(1,36)]+[str(i) 
                                           for i in range(101,108)]
        flist = os.listdir(file_loc)
        files = []
        for ens in ensembles:
            prefix = ['.'.join(['b.e11.B20TRC5CNBDRD.f09_g16', 
                                ens, model, fcode, variable]),
                  '.'.join(['b.e11.BRCP85C5CNBDRD.f09_g16', 
                            ens, model, fcode, variable])]
            files += [f for f in flist if f.startswith(prefix[0])] 
            files += [f for f in flist if f.startswith(prefix[1])]
        flist =  [file_loc + '/' + f for f in files]
        flist.sort()
        return flist
    
    def get_cesm2_file_names(variable, freq, **kwargs):
        """Reads filenames for cesm2 from glade
        and returns as a master list."""

        if kwargs['atm_model'] == 'cam6':
            model = 'CESM2'
        elif kwargs['atm_model'] == 'waccm':
            model = 'CESM2-WACCM'
        
        assert freq == 'monthly', 'Frequency must be monthly, for now.'

        if variable in ['siconca']:
            domain = 'ice'
            dfreq = 'SImon'
        else:
            domain = 'atm'
            dfreq = 'Amon'
            
        file_loc = '/'.join(['/glade/collections/cmip/CMIP6/CMIP/NCAR',
                             model, 'historical'])
        ensembles = os.listdir(file_loc) 
        flist = []
        for ens in ensembles:
            floc = '/'.join([file_loc, ens, dfreq,
                                         variable,'gn','latest'])
            
            files = os.listdir(floc)
            
            flist += [floc + '/' + f for f in files]# if f.endswith('.nc')]
            

        flist.sort()
        return flist
    
    
    
    ###################### MAIN PART ##############################
    if source=='cesm-le':
        varname = get_varnames(variable, source, freq)
        return get_cesmle_file_names(varname, freq)
    
    elif source=='cesm2-cam6':
        varname = get_varnames(variable, source, freq)
        return get_cesm2_file_names(varname, freq, atm_model='cam6')
    
    elif source=='cesm2-waccm':
        varname = get_varnames(variable, source, freq)
        return get_cesm2_file_names(varname, freq, atm_model='waccm')
    
    elif source=='era-i':
        varname = get_varnames(variable, source, freq)
        return get_erai_file_names(varname, freq)
    
    elif source=='era5':
        varname = get_varnames(variable, source, freq)
        return get_era5_file_names(varname, freq)



# def get_cesm2_file_names(ens, variable, freq, domain):

# def get_erai_file_names(variable, freq, domain):

# def get_era5_file_names(variable, freq, domain):

def get_varnames(variable, source, freq):
    """Given a variable name and a source, provides the name of 
    the corresponding variable in the source data if it exists.
    Including frequency as an input because I think the ERA data
    are different at different time resolutions sometimes.
    """
    
    cesmle = {'latitude': 'LAT',
              'longitude': 'LON',
              'land_mask': 'LANDFRAC',
              'sea_ice_concentration': 'ICEFRAC',
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
              
    cmip6 =  {'latitude': 'lat',
              'longitude': 'lon',
              'land_mask': np.nan,
              'sea_ice_concentration': 'siconca',
              'total_cloud_cover': 'clt',
              'low_cloud_fraction': np.nan,
              '2m_temperature': 'tas',
              'air_temperature': 'ta',
              'sensible_heat_flux': 'hfls',
              'latent_heat_flux': 'hfss',
              '10m_wind_speed': 'sfcWind',
              'sea_level_pressure': 'psl',
              'surface_downward_longwave': 'rlds',
              'surface_upward_longwave': 'rlus',
              'surface_downward_shortwave': 'rsds',
              'surface_upward_shortwave': 'rsus',
              'net_surface_longwave': np.nan,
              'net_surface_shortwave': np.nan}
              
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
              'net_surface_shortwave': 'SSR_GDS4_SFC_120'}   
              
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
               'era-i': erai,
               'era5': era5}
    
    if source not in sources.keys():
        print('Source must be one of the following: \n----------------------- \n' + '\n'.join(list(sources.keys())))
    elif variable not in cesmle.keys():
        print('Variable must be one of the following: \n----------------------- \n'  + '\n'.join(list(cesmle.keys())))
    else:
        return sources[source][variable]


