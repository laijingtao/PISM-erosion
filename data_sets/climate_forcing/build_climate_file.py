#!/usr/bin/env python
# Copyright (C) 2017 Jingtao Lai

import numpy as np
from netCDF4 import Dataset

infile = '../bed_dem/pism_Synthetic_v1.nc'
outfile = 'constant_climate.nc'

# dimensions of Olympic
x_dim = 200
y_dim = 180
dim_value = {'x': np.arange(x_dim),
             'y': np.arange(y_dim)}
if not(infile is None):
    indata = Dataset(infile, 'r')
    x_dim, = indata.variables['x'].shape
    y_dim, = indata.variables['y'].shape
    dim_value['x'] = indata.variables['x'][:]
    dim_value['y'] = indata.variables['y'][:]
    indata.close()

air_temp_mean_annual = 9.5
air_temp_mean_july = 15.5
precipitation = 1500.0

outdata = Dataset(outfile, 'w')

dim_len = {'x': x_dim, 'y': y_dim, 'time': None}
for dim_name in ['x', 'y', 'time']:
    try:
        outdata.createDimension(dim_name, dim_len[dim_name])
        outdata.createVariable(dim_name, np.float64, (dim_name,))
        if not(dim_len[dim_name] is None):
            outdata.variables[dim_name][:] = dim_value[dim_name]
    except:
        continue

var_name_list = ['air_temp_mean_annual',
                 'air_temp_mean_july',
                 'precipitation']
var_value = {'air_temp_mean_annual': air_temp_mean_annual,
             'air_temp_mean_july': air_temp_mean_july,
             'precipitation': precipitation}
var_unit = {'air_temp_mean_annual': 'celsius',
             'air_temp_mean_july': 'celsius',
             'precipitation': 'mm yr-1'}
for var_name in var_name_list:
    try:
        var = outdata.variables[var_name]
    except:
        var = outdata.createVariable(var_name, np.float64, ('y', 'x',))
    var[:] = var_value[var_name]
    var.units = var_unit[var_name]

outdata.close()
