#!/usr/bin/env python
# Copyright (C) 2017 Jingtao Lai

import numpy as np
from netCDF4 import Dataset

infile = 'test.nc'
outfile = 'test_climate.nc'

# dimensions of Olympic
x_dim = 200
y_dim = 180

air_temp_mean_annual = 1.0
air_temp_mean_july = 2.0
precipitation = 3.0

outdata = Dataset('test_climate.nc', 'w')

dim_len = {'x': x_dim, 'y': y_dim, 'time': None}
for dim_name in ['x', 'y', 'time']:
    try:
        outdata.createDimension(dim_name, dim_len[dim_name])
        outdata.createVariable(dim_name, np.float64, (dim_name,))
        outdata.variables[dim_name][:] = np.arange(dim_len[dim_name])
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
             'precipitation': 'kg m-2 yr-1'}
for var_name in var_name_list:
    try:
        var = outdata.variables[var_name]
    except:
        var = outdata.createVariable(var_name, np.float64, ('y', 'x',))
    var[:] = var_value[var_name]
    var.units = var_unit[var_name]
