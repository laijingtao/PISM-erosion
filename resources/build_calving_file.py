#!/usr/bin/env python
# Copyright (C) 2017 Jingtao Lai

import os
import numpy as np
from netCDF4 import Dataset
import subprocess as sub

def build_ocean_kill_file(infile=None, outfile=None, *args, **kwargs):
    try:
        thk = kwargs['thk']
    except:
        thk = 0
    
    write_data(infile, outfile, thk=thk)

def write_data(infile=None, outfile=None, *args, **kwargs):
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
    outdata = Dataset(outfile, 'w')

    thk = kwargs['thk']

    dim_len = {'x': x_dim, 'y': y_dim, 'time': None}
    for dim_name in ['x', 'y', 'time']:
        try:
            outdata.createDimension(dim_name, dim_len[dim_name])
            outdata.createVariable(dim_name, np.float64, (dim_name,))
            if not(dim_len[dim_name] is None):
                outdata.variables[dim_name][:] = dim_value[dim_name]
        except:
            continue

    var_name_list = ['thk']
    var_value = {'thk': thk}
    var_unit = {'thk': 'm'}
    for var_name in var_name_list:
        try:
            var = outdata.variables[var_name]
        except:
            var = outdata.createVariable(var_name, np.float64, ('y', 'x',))
        var[:] = var_value[var_name]
        var.units = var_unit[var_name]

    outdata.variables['thk'].standard_name = 'land_ice_thickness'

    outdata.close()

if __name__=='__main__':
    infile = '../data_sets/bed_dem/test_dem.nc'
    outfile = 'test_ocean_kill.nc'
    build_ocean_kill_file(infile=infile, outfile=outfile) 
