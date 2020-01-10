"""Translates variable names between source datasets."""
import numpy as np

def varname(name, source):
    """Given a variable name and a source, provides the name of the corresponding variable in the 
    source data if it exists."""
    
    cesmle = {'latitude': 'LAT',
              'longitude': 'LON',
              'sea_ice_concentration': 'SIC',
              'total_cloud_cover': 'TCC',
              'low_cloud_fraction': 'LCC',
              '2m_temperature': 'TREFHT',
              'sensible_heat_flux': 'SHFLX',
              'latent_heat_flux': 'LHFLX',
              '10m_wind_speed': 'U10',
              'sea_level_pressure': 'PSL',
              'surface_downward_longwave': np.nan,
              'surface_upward_longwave': np.nan,
              'surface_downward_shortwave': np.nan,
              'surface_upward_shortwave': np.nan,
              'net_surface_longwave': 'FLNS',
              'net_surface_shortwave': 'FSNS'}
              
    cmip6 =  {'latitude': 'lat',
              'longitude': 'lon',
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
              
    

