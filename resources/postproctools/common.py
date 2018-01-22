#!/usr/bin/env python
# Copyright (C) 2017 Jingtao Lai

import os
import sys
import numpy as np
try:
    import subprocess32 as sub
except:
    import subprocess as sub
from netCDF4 import Dataset

FILL_VALUE = -2e9

def dist2lonlat(easting=None, northing=None, *args, **kwargs):
    # x/y --> lon/lat
    # default projection is for olympic mountains

    from pyproj import Proj

    try:
        projection = kwargs['projection']
    except:
        projection = 'epsg:26710' 

    proj = Proj(init=projection)
    ee, nn = np.meshgrid(easting, northing)
    lon, lat = proj(ee, nn, inverse=True)
    return lon, lat

def nc_reorder_xy_nco(infile=None, outfile=None):
    # In pism 0.7, the dimension order is (x, y)
    # This function will change it to (y, x)
    if infile is None:
        sys.exit('Must provide an input!')
    if outfile is None:
        outfile = infile
    cmd = ['ncpdq', '-a', 'y,x', '-O', infile, outfile]
    sub.call(cmd)

def nc_copy_dim(infile=None, outfile=None):
    # copy dimensions
    indata = Dataset(infile, 'r')
    if os.path.isfile(outfile):
        outdata = Dataset(outfile, 'a')
    else:
        outdata = Dataset(outfile, 'w')

    for dim_name, dim in indata.dimensions.iteritems():
        try:
            outdata.createDimension(dim_name, len(dim) if not dim.isunlimited() else None)
        except:
            continue
        try:
            in_dim_var = indata.variables[dim_name]
        except:
            continue
        out_dim_var = outdata.createVariable(dim_name, in_dim_var.datatype,
                                             in_dim_var.dimensions)
        out_dim_var.setncatts(
            {k: in_dim_var.getncattr(k) for k in in_dim_var.ncattrs()})
        out_dim_var[:] = in_dim_var[:]
    indata.close()
    outdata.close()

def nc_copy_var(infile=None, outfile=None, var=None):
    # copy variables
    if var is None:
        sys.exit('Must provide at least one variable name!')
    if type(var) is str:
        var = [var]
    indata = Dataset(infile, 'r')
    outdata = Dataset(outfile, 'a')
  
    for var_name in var:
        if var_name in indata.dimensions.keys():
            continue
        in_var = indata.variables[var_name]
        out_var = outdata.createVariable(var_name, in_var.datatype,
                                         in_var.dimensions)
        out_var.setncatts(
            {k: in_var.getncattr(k) for k in in_var.ncattrs()})
        out_var[:] = in_var[:]
    '''
    # copy all variables
    for var_name, in_var in indata.variables.iteritems():
        if var_name in indata.dimensions.keys():
            continue
        out_var = outdata.createVariable(var_name, in_var.datatype,
                                         in_var.dimensions)
        out_var.setncatts(
            {k: in_var.getncattr(k) for k in in_var.ncattrs()})
        out_var[:] = in_var[:]
    '''
    indata.close()
    outdata.close()

def nc_regrid_dem(infile=None, outfile=None, grid=None, *args, **kwargs):
    # Based on gdalwarp and nco
    try:
        fill_value = kwargs['fill_value']
    except:
        fill_value = -2.e9
    try:
        resample = kwargs['resample']
    except:
        resample = 'average'
    # if infile's dimension has units, gdalwarp won't work
    cmd = ['ncatted', '-a', 'units,x,d,,', '-a', 'units,y,d,,',
           infile, 'tmp.nc']
    sub.call(cmd)
    cmd = ['gdalwarp', '-of', 'netCDF', '-overwrite', 
           '-tr', str(grid), str(grid), '-srcnodata', str(fill_value),
           '-dstnodata', str(fill_value), '-r', resample, 'tmp.nc', outfile]
    print ' '.join(cmd)
    sub.call(cmd)
    cmd = ['ncrename', '-v', 'Band1,topg', '-v', 'lon,x', '-v', 'lat,y', outfile]
    print ' '.join(cmd)
    sub.call(cmd)
    cmd = ['ncrename', '-d', 'lon,x', '-d', 'lat,y', outfile]
    print ' '.join(cmd)
    sub.call(cmd)
    cmd = ['ncatted', '-a', 'units,x,o,c,m', '-a', 'units,y,o,c,m',
           '-a', 'standard_name,x,o,c,projection_x_coordinate', '-a', 'long_name,x,d,,', 
           '-a', 'standard_name,y,o,c,projection_y_coordinate', '-a', 'long_name,y,d,,', 
           '-a', 'long_name,topg,d,,', '-a', '_FillValue,topg,o,f,{}'.format(fill_value),
           '-a', 'standard_name,topg,o,c,bedrock_altitude',
           outfile]
    print ' '.join(cmd)
    sub.call(cmd)
    cmd = ['rm', 'tmp.nc']
    sub.call(cmd)

def get_grid_size(infile):
    indata = Dataset(infile, 'r')
    try:
        indata.variables['x']
    except:
        sys.exit('Couldn\'t find x dimension')
    try:
        indata.variables['y']
    except:
        sys.exit('Couldn\'t find y dimension')

    dx = abs(indata.variables['x'][0]-indata.variables['x'][1])
    dy = abs(indata.variables['y'][0]-indata.variables['y'][1])

    return dx, dy
