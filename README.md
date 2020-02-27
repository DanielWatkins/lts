# lts
Code for analyzing lower tropospheric stability on the NCAR supercomputing system. It provides tools to access datasets stored on glade, and requires that the user be running code on a computer with access to glade. This package provides the code required to reproduce the results of my paper. 

Data used:
* CESM1-CAM5
* CESM1-WACCM
* CESM2-CAM6
* CESM2-WACCM

Planned
* ERA-Interim
* ERA5
* CESM Large Ensemble

# Usage
The first step is to create a Parameters object for each of the configurations desired. The Parameters object is passed around to various parts of the code, and keeps track of begin/end times, source dataset, lat/lon subsets, and variables needed.

````
source = 'cesm1-cam5' # or cesm1-waccm, cesm2-cam6, cesm2-waccm

params = lts.Parameters(source)
````

To see which variables are available, you can use the method `params.list_available_variables()`. (not yet implemented)
You can see the currently selected parameters by calling 
`params.display()`.
Choose a save location large enough to store the subsetted data. I made directories
for each of the sources, so I use
`params.save_location = '/glade/scratch/dwatkins/' + params1.source + '/'`

The variables parameter must be a list, and has the near-natural-language names of variables in it. 
`params.variables = ['air_temperature',
                     '2m_temperature',
                     'surface_downward_longwave',
                     'eastward_wind',
                     'northward_wind',
                     'condensed_ice_path',
                     'condensed_water_path',
                     'sensible_heat_flux',
                     'sea_level_pressure']
                     `

Next, by calling get_data(), datasets are read from glade, subsetted, and saved locally. 

## To Do
Within get_data or load_data, I want lts to get computed. 
There could also be an option to warn if overwriting existing data.
I need to test all the variable options.

Next steps: 
- build a function to load locally saved data
- make sure I can load in the cloud data
- make functions for histograms and timeseries analysis
- make sure all the 4 datasets can load