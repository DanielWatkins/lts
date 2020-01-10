"""Helper functions to access files via the CESM variable names.
get_filenames(variable, freq, source, domain=None, ens=None)
get_cesmle_file_names(ens, variable, freq, domain)



"""
import os
import pandas as pd

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
        return get_cesmle_file_names(ens, variable, freq, domain)
    elif source=='CESM2-CAM6':
        return get_cesm2_file_names(ens, variable, freq, domain)
    elif source=='CESM2-WACCM':
        return get_cesm2_file_names(ens, variable, freq, domain)
    elif source=='ERAI':
        return get_erai_file_names(variable, freq, domain)
    elif source=='ERA5':
        return get_era5_file_names(variable, freq, domain)


def get_cesmle_file_names(ens, variable, freq, domain):
    """Returns list of file names matching the ensemble
    number, variable, and frequency."""

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

# def get_cesm2_file_names(ens, variable, freq, domain):

# def get_erai_file_names(variable, freq, domain):

# def get_era5_file_names(variable, freq, domain):
