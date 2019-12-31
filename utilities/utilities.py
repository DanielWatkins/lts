"""Utility functions supporting paper on LTS.

Contents

Global variables:
months: list of month names

Miscellaneous functions:
check_time: tests whether the dates listed in the filename overlap the study period

Data access helpers:
get_filenames(variable, freq, source, domain, ens): returns list of filenames. Sources accepted include
    CESMLE, *ERA-I, *ERA5, *CESM2-CAM6, and *CESM2-WACCM.
    
*get_varnames(variable, source): returns source-specific variable names based on CESM variable name.




get_cesm_file_names(ens, variable, freq, domain): returns a list of filenames 
                                                  matching the ensemble number, variable, frequency, and domain.
*get_erai_file_names(var, daterange)
*get_era5_file_names(var, daterange)

*monthly_maps()

area_weights: check difference between my method and the 1/cos phi method.


*planned

(c) Daniel Watkins 2019
"""

months = ['January', 'February', 'March', 'April', 'May', 'June', 'July',
          'August', 'September', 'October', 'November', 'December']

def check_time(file, st, et):
    """Returns true if the dates in file overlap the range [st, et]."""
    import pandas as pd
    s_e = file.split('.')[-2].split('-')
    s = pd.to_datetime(s_e[0], format='%Y%M')
    e = pd.to_datetime(s_e[1], format='%Y%M') # Different for daily output
    st = pd.to_datetime(st)
    et = pd.to_datetime(et)
    t1 = (s <= st) & (st <= e)
    t2 = (s <= et) & (et <= e)
    t3 = (st <= s) & (e <= et)
    return (t1 | t2) | t3

def get_filenames(variable, freq, source, domain=None, ens=None):
    """Returns list of file names matching the ensemble
    number, variable, and frequency.
    
    variable=variable name in the source data set
    freq=daily or monthly or 6-hourly
    domain=atmosphere or ice, probably not used for ERA
    source=ERAI, ERA5, CESMLE, CESM2-CAM6, CESM2-WACCM
    ens=string eg '001'
    """
    
    if source=='CESMLE':
        return get_cesm_file_names(ens, variable, freq, domain)
    elif source=='CESM2-CAM6':
        return get_cesm2_file_names(ens, variable, freq, domain)
    elif source=='CESM2-WACCM':
        return get_cesm2_file_names(ens, variable, freq, domain)
    elif source=='ERAI':
        return get_erai_file_names(variable, freq, domain)
    elif source=='ERA5':
        return get_era5_file_names(variable, freq, domain)

def get_erai_file_names(variable):
    """Gets the ERA-I file name list for variables corresponding to the
    CESM variable name standard."""
    
    prefix_mdfa = 'sfc'.join(['/gpfs/fs1/collections/rda/data/ds627.1/ei.mdfa.fc12hr.',
                          '/ei.mdfa.fc12hr.',
                          '.regn128sc.'])

    prefix_moda_sfc = 'sfc'.join(['/gpfs/fs1/collections/rda/data/ds627.1/ei.moda.an.',
                      '/ei.moda.an.',
                      '.regn128sc.'])

    prefix_moda_pl = 'pl'.join(['/gpfs/fs1/collections/rda/data/ds627.1/ei.moda.an.',
                      '/ei.moda.an.',
                      '.regn128sc.'])

    erai_prefix =   {'ICEFRAC': prefix_moda_sfc,
                    'CLDTOT': prefix_moda_sfc,
                    'CLDLOW': prefix_moda_sfc,
                    'TREFHT': prefix_moda_sfc,
                    'U10': prefix_moda_sfc,
                    'SHFLX': prefix_mdfa,
                    'LHFLX': prefix_mdfa,
                    'FLNS': prefix_mdfa,
                    'FSNS': prefix_mdfa,
                    'PSL': prefix_mdfa,
                    'T850': prefix_moda_pl}

    erai_varnames = {'ICEFRAC': 'CI_GDS4_SFC_S123',
                    'CLDTOT': 'TCC_GDS4_SFC_S123',
                    'CLDLOW': 'LCC_GDS4_SFC_S123',
                    'TREFHT': '2T_GDS4_SFC_S123',
                    'U10': '10SI_GDS4_SFC_S123',
                    'SHFLX': 'SSHF_GDS4_SFC_120',
                    'LHFLX': 'SLHF_GDS4_SFC_120',
                    'FLNS': 'STR_GDS4_SFC_120',
                    'FSNS': 'SSR_GDS4_SFC_120',
                    'PSL': 'MSL_GDS4_SFC_S123',
                    'T850': 'T_GDS4_ISBL_S123'}
    
    dates = pd.date_range('1979-01-01', '2018-12-31', freq='MS')
    fnames = [erai_prefix[variables[0]] + d.strftime('%Y%m%d00') for d in dates]
    
    
    
def gather_erai_data(variables, saveloc):
    """Selects the ERA-I data for 60-90 latitude for the variables in variables and
    writes netCDF files for each variable in saveloc. Assumes that all the variables
    in the list come from the same source data file. Currently only pulls the isobaric
    data for temperature at 850 hPa.
    
    Recommended usage:
    
    varset1 = ['T850']
    varset2 = ['ICEFRAC', 'CLDTOT', 'CLDLOW', 'TREFHT', 'U10', 'PSL']
    varset3 = ['SHFLX', 'LHFLX', 'FLNS', 'FSNS']
    saveloc = <path with / at the end>
    
    gather_erai_data(varset1, saveloc)
    gather_erai_data(varset2, saveloc)
    gather_erai_data(varset3, saveloc)
    
    The grouping of the variables above is based on which GRIB file contains the
    variables, since the entire file must be loaded at once.
    
    # Todo: copy attributes, make it CF compliant, add variable name
    """
    
    import pandas as pd
    import xarray as xr
    import numpy as np
    import Nio
    import os
    
    
    #### Set up the file locations and variable names ####
    prefix_mdfa = 'sfc'.join(['/gpfs/fs1/collections/rda/data/ds627.1/ei.mdfa.fc12hr.',
                          '/ei.mdfa.fc12hr.',
                          '.regn128sc.'])

    prefix_moda_sfc = 'sfc'.join(['/gpfs/fs1/collections/rda/data/ds627.1/ei.moda.an.',
                      '/ei.moda.an.',
                      '.regn128sc.'])

    prefix_moda_pl = 'pl'.join(['/gpfs/fs1/collections/rda/data/ds627.1/ei.moda.an.',
                      '/ei.moda.an.',
                      '.regn128sc.'])

    erai_prefix =   {'ICEFRAC': prefix_moda_sfc,
                    'CLDTOT': prefix_moda_sfc,
                    'CLDLOW': prefix_moda_sfc,
                    'TREFHT': prefix_moda_sfc,
                    'U10': prefix_moda_sfc,
                    'SHFLX': prefix_mdfa,
                    'LHFLX': prefix_mdfa,
                    'FLNS': prefix_mdfa,
                    'FSNS': prefix_mdfa,
                    'PSL': prefix_mdfa,
                    'T850': prefix_moda_pl}

    erai_varnames = {'ICEFRAC': 'CI_GDS4_SFC_S123',
                    'CLDTOT': 'TCC_GDS4_SFC_S123',
                    'CLDLOW': 'LCC_GDS4_SFC_S123',
                    'TREFHT': '2T_GDS4_SFC_S123',
                    'U10': '10SI_GDS4_SFC_S123',
                    'SHFLX': 'SSHF_GDS4_SFC_120',
                    'LHFLX': 'SLHF_GDS4_SFC_120',
                    'FLNS': 'STR_GDS4_SFC_120',
                    'FSNS': 'SSR_GDS4_SFC_120',
                    'PSL': 'MSL_GDS4_SFC_S123',
                    'T850': 'T_GDS4_ISBL_S123'}
    
    dates = pd.date_range('1979-01-01', '2018-12-31', freq='MS')
    fnames = [erai_prefix[variables[0]] + d.strftime('%Y%m%d00') for d in dates]
    
    #### Get latitude, longitude, and pressure level indices ####

    if variables[0] == 'T850':
        latname = 'g4_lat_1'
        lonname = 'g4_lon_2'
    else:
        latname = 'g4_lat_0'
        lonname = 'g4_lon_1'
        
    ds = Nio.open_file(fnames[0], 'r')
    if variables[0] == 'T850':
        idx_level = np.where(ds.variables['lv_ISBL0'][:] == 850)[0][0]

    idx_lat = np.where(ds.variables[latname][:] >= 59)[0]
    min_lat = min(idx_lat).squeeze()
    max_lat = max(idx_lat).squeeze()

    lat = ds.variables[latname][min_lat:max_lat]
    lon = ds.variables[lonname][:]
    ds.close()

    #### Loop through the variables, convert to numpy arrays, then save to netcdf ####
    data_dict = {var:[] for var in variables}

    if variables[0] == 'T850':
        for fname in fnames:
            ds = Nio.open_file(fname, 'r')
            vardata = ds.variables[erai_varnames['T850']][idx_level,min_lat:max_lat,:]
            data_dict['T850'].append(vardata)
            ds.close()
            
        data = np.array(data_dict['T850'])
        da = xr.DataArray(data=data, dims={'time':dates, 'lat':lat, 'lon':lon})
        da.to_netcdf(saveloc + 'erai.monthly.T850.nc')
        return
    
    else:
        for fname in fnames:
            try:
                ds = Nio.open_file(fname, 'r')
                for var in variables:
                    vardata = ds.variables[erai_varnames[var]][min_lat:max_lat,:]
                    data_dict[var].append(vardata)
                ds.close()
            except:
                print(fname)
                
        for var in variables:
            data = np.array(data_dict[var])
            da = xr.DataArray(data={var: data}, dims={'time':dates, 'lat':lat, 'lon':lon})
            da.coords['lat'] = lat
            da.coords['lon'] = lon
            da.to_netcdf(saveloc + 'erai.monthly.' + var + '.nc')
        return        
        
def gather_cesmle_data(variables, save_loc):
    """Selects data north of 59 and for the 1979-2018 time range and 
    saves it in saveloc."""
    
    import numpy as np
    import xarray as xr
    import pandas as pd

    ensembles = [str(i).zfill(3) for i in range(1,36)]+[str(i) for i in range(101,106)] 
    domain = 'atm'
    freq = 'monthly'
    start_time = '1979-01-01'
    end_time = '2018-12-31'
    
    for variable in variables:
        da_all = []
        for ens in ensembles:
            files = get_cesm_file_names(ens, variable, freq, domain=domain)
            ds = []
            for f in files:
                if check_time(f, start_time, end_time):
                    with xr.open_dataset(f) as d:
                        # Important to interpolate first, because slight differences in the 
                        # roundoff error in lat/lon lead to duplication of coordinates
                        # in xarray's concat
                        new_lat = np.round(d.lat.values, 3)
                        new_lon = np.round(d.lon.values, 3)
                        new_lat = new_lat[new_lat > 59]
                        data = d.sel(lat=slice(55,90), time=slice(start_time, end_time)).load()
                        data = data.interp(lat=new_lat, lon=new_lon, kwargs={'fill_value':None})
                        ds.append(data)
                        del data
            ds_full = xr.concat(ds, dim='time')
            da_all.append(ds_full[variable])
            del ds
            del ds_full
        da = xr.concat(da_all, dim=pd.Index(ensembles, name='ens'))
        ds = xr.Dataset(data_vars={variable: (['ens', 'time', 'lat', 'lon'], da.values)},
                   coords={'ens': ensembles, 
                           'time': pd.date_range(start_time, end_time, freq='MS'),
                           'lat': new_lat,
                           'lon': new_lon})
        ds.to_netcdf(save_loc + 'cesmle.monthly.' + variable + '.nc')
        print(variable)
        
        

def get_cesm_file_names(ens, variable, freq, domain):
    """Returns list of file names matching the ensemble
    number, variable, and frequency."""
    import os
    import pandas as pd
    
    file_loc = '/'.join(['/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE',
                         domain, 'proc/tseries', freq, variable])

    if domain=='atm':
        model = 'cam'
    elif domain=='ice':
        model = 'cice'
    else:
        print('Domain must be atm or ice')
    
    if domain=='atm':
        if freq=='daily':
            fcode = 'h1'
        elif freq == 'monthly':
            fcode = 'h0'
    elif domain == 'ice':
        if freq == 'monthly':
            fcode = 'h'
        if freq == 'daily':
            # variable for the folder is var + _d
            # variables for northern hemisphere add _nh
            fcode = 'h1'
            file_loc += '_d'
        variable += '_nh'
        
    prefix = ['.'.join(['b.e11.B20TRC5CNBDRD.f09_g16', ens, model, fcode, variable]),
              '.'.join(['b.e11.BRCP85C5CNBDRD.f09_g16', ens, model, fcode, variable])]
    
    flist = os.listdir(file_loc)
    files = [f for f in flist if f.startswith(prefix[0])] 
    files += [f for f in flist if f.startswith(prefix[1])]
    flist =  [file_loc + '/' + f for f in files]
    flist.sort()
    return flist

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
                        time=slice(begin_time,end_time)).variables['time'].values,
                           'lat': ds_var.lat.values[ds_var.lat.values >= min_lat],
                           'lon': ds_var.lon.values,
                           'lev': p_levels})

    return ds
    #ds.to_netcdf(out_file)
    
    
def weights2(lats, lons, min_lat, fraction=None):
    """Compute a weight for each lat/lon grid box. Assumes that
    the polar grid points are offset (so values are centered halfway
    between 90 and the next highest lat), and the rest of the grid points
    represent grid centers.

    Fraction (optional) is an array specifying the portion of the grid cell that "counts".
    For example, LANDFRAC for land, 1-LANDFRAC-ICEFRAC for ocean, and ICEFRAC. 
    """

    import numpy as np
    
    dlat = np.diff(lats)[0]/2
    nlon = len(lons)
    lats = lats.copy()
    lats[0] = lats[0] + dlat
    lats[-1] = lats[-1] - dlat
    lats = lats[lats >= min_lat]

    sin_lat1 = np.sin((lats + dlat) * np.pi / 180.)
    sin_lat2 = np.sin((lats - dlat) * np.pi / 180.)
    sin_lat = np.abs(sin_lat1 - sin_lat2)
    weights = (np.zeros((len(lats), nlon)) + np.array([sin_lat / (np.sum(sin_lat)*nlon)]).T)

    if fraction is not None:
        weights = fraction * weights
        
    weights = weights / np.nansum(weights)
    return weights


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
    
def average(data, dim=None, weights=None):
    """
    weighted average for xray objects

    Code by Github user mathause, from xarray issue #422
    https://github.com/pydata/xarray/issues/422
    accessed Aug 14, 2019

    Parameters
    ----------
    data : Dataset or DataArray
        the xarray object to average over
    dim : str or sequence of str, optional
        Dimension(s) over which to apply average.
    weights : DataArray
        weights to apply. Shape must be broadcastable to shape of data.

    Returns
    -------
    reduced : Dataset or DataArray
        New xarray object with average applied to its data 
        and the   indicated dimension(s) removed.

    """
    import xarray as xr
    
    def average_da(self, dim=None, weights=None):
        """
        weighted average for DataArrays

        Parameters
        ----------
        dim : str or sequence of str, optional
        Dimension(s) over which to apply average.
        weights : DataArray
            weights to apply. Shape must be broadcastable to shape of self.

        Returns
        -------
        reduced : DataArray
            New DataArray with average applied to its data and the indicated
            dimension(s) removed.

        """

        if weights is None:
            return self.mean(dim)
        else:
            if not isinstance(weights, xr.DataArray):
                raise ValueError("weights must be a DataArray")

            # if NaNs are present, we need individual weights
            if self.notnull().any():
                total_weights = weights.where(self.notnull()).sum(dim=dim)
            else:
                total_weights = weights.sum(dim)

            return (self * weights).sum(dim) / total_weights

    def average_ds(self, dim=None, weights=None):
        """
        weighted average for Datasets

        Parameters
        ----------
        dim : str or sequence of str, optional
            Dimension(s) over which to apply average.
        weights : DataArray
            weights to apply. Shape must be broadcastable to shape of data.

        Returns
        -------
        reduced : Dataset
            New Dataset with average applied to its data and the indicated
            dimension(s) removed.

        """

        if weights is None:
            return self.mean(dim)
        else:
            return self.apply(average_da, dim=dim, weights=weights)



    if isinstance(data, xr.Dataset):
        return average_ds(data, dim, weights)
    elif isinstance(data, xr.DataArray):
        return average_da(data, dim, weights)
    else:
        raise ValueError("date must be an xarray Dataset or DataArray")
