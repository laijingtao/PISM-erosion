#!/usr/bin/env python
# Copyright (C) 2017 Jingtao Lai

import numpy as np
from netCDF4 import Dataset

infile = 'test.nc'
outfile = 'test_climate.nc'

outdata = Dataset('test_climate.nc', 'a')

# dimensions of Olympic
x_dim = 200
y_dim = 180

dim_len = {'x': x_dim, 'y': y_dim, 'time': None}
for dim_name in ['x', 'y', 'time']:
    try:
        outdata.createDimension(dim_name, dim_len[dim_name])
        outdata.createVariable(dim_name, np.float64, (dim_name,))
        outdata.variables[dim_name][:] = np.arange(dim_len[dim_name])
    except:
        continue

var_value = {'air_temp_mean_annual': air_temp_mean_annual,
             'air_temp_mean_july': air_temp_mean_july,
             'precipitation': precipitation}

