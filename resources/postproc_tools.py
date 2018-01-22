#!/usr/bin/env python
# Copyright (C) 2017 Jingtao Lai

import numpy as np
import os
import sys
import itertools
try:
    import subprocess32 as sub
except:
    import subprocess as sub
from collections import OrderedDict
from argparse import ArgumentParser
from netCDF4 import Dataset


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

def calc_erosion_nco(infile=None, outfile=None):
    if infile is None:
        sys.exit('Must provide an input file!')
    if outfile is None:
        outfile = infile
    
    cmd = ['ncap2', '-O', '-s',
           'erosion_1=1e-4*velbase_mag;erosion_2=2.7e-7*velbase_mag^2.02;',
           infile, outfile]
    #print sub.list2cmdline(cmd)
    sub.call(cmd)

def calc_erosion(infile=None, outfile=None):
    fill_value = -2.0e9
    if infile is None:
        sys.exit('Must provide an input file!')
    overwrite = False
    if outfile is None:
        outfile = infile
        overwrite = True
    else:
        nc_copy_dim(infile, outfile)

    if overwrite:
        indata = Dataset(infile, 'a')
        outdata = indata
    else:
        indata = Dataset(infile, 'r')
        outdata = Dataset(outfile, 'a')

    time = indata.variables['time'][:]
    grid = abs(indata.variables['x'][0]-indata.variables['x'][1])
    col = len(indata.variables['x'])
    row = len(indata.variables['y'])

    for erosion_name in ['erosion_1', 'erosion_2']:
        if erosion_name=='erosion_1':
            erosion = 1e-4*indata.variables['velbase_mag'][:]
        if erosion_name=='erosion_2':
            erosion = 2.7e-7*np.power(indata.variables['velbase_mag'][:], 2.02)
        erosion[np.where(erosion<=0)] = fill_value
        erosion = np.ma.masked_values(erosion, fill_value)
        try:
            erosion_var = outdata.createVariable(
                erosion_name, np.float64, ('time', 'y', 'x',),
                fill_value=fill_value)
        except:
            erosion_var = outdata.variables[erosion_name]
            erosion_var._FillValue = fill_value
        erosion_var[:] = erosion
        erosion_var.units = 'm year-1'

    indata.close()
    if not overwrite:
        outdata.close()

def calc_total_erosion_nco(infile=None, outfile=None):
    if infile is None:
        sys.exit('Must provide an input file!')
    if outfile is None:
        outfile = infile
    dx, dy = get_grid_size(infile)
    cmd = ['ncap2', '-s',
           'total_erosion_1=erosion_1.total($x,$y)*'+str(dx)+'*'+str(dx)+\
           ';total_erosion_2=erosion_2.total($x,$y)*'+str(dx)+'*'+str(dx),
           '-A', infile, outfile]
    sub.call(cmd)
    cmd = ['ncatted', '-a', 'units,total_erosion_1,o,c,"m3 year-1"',
           '-a', 'units,total_erosion_2,o,c,"m3 year-1"', outfile]
    sub.call(cmd)

def calc_erosion_space_averaged(infile=None, outfile=None):
    if infile is None:
        sys.exit('Must provide an input file!')
    overwrite = False
    if outfile is None:
        outfile = infile
        overwrite = True
    else:
        nc_copy_dim(infile, outfile)

    if overwrite:
        indata = Dataset(infile, 'a')
        outdata = indata
    else:
        indata = Dataset(infile, 'r')
        outdata = Dataset(outfile, 'a')

    time = indata.variables['time'][:]
    grid = abs(indata.variables['x'][0]-indata.variables['x'][1])
    col = len(indata.variables['x'])
    row = len(indata.variables['y'])

    for erosion_name in ['erosion_1', 'erosion_2']:
        erosion = indata.variables[erosion_name][:]
        total_erosion = np.array(
            [erosion_slice.sum() for erosion_slice in erosion])*grid*grid
        erosion_space_averaged = total_erosion/(grid*grid*row*col)
        try:
            erosion_space_averaged_var = outdata.createVariable(
                erosion_name+'_space_averaged', np.float64, ('time',))
        except:
            erosion_space_averaged_var = outdata.variables[erosion_name+'_space_averaged']
        erosion_space_averaged_var[:] = erosion_space_averaged
        erosion_space_averaged_var.units = 'm year-1'

    indata.close()
    if not overwrite:
        outdata.close()

def calc_percent_of_time_covered_by_ice(infile=None, outfile=None):
    if infile is None:
        sys.exit('Must provide an input file!')
    overwrite = False
    if outfile is None:
        outfile = infile
        overwrite = True
    else:
        nc_copy_dim(infile, outfile)

    if overwrite:
        indata = Dataset(infile, 'a')
        outdata = indata
    else:
        indata = Dataset(infile, 'r')
        outdata = Dataset(outfile, 'a')
    
    time = indata.variables['time'][:]/(365*24*3600.)
   
    var_list = ['velbase_mag', 'erosion_1', 'erosion_2']
    flag = True
    for var in var_list:
        if var in indata.variables.keys():
            flag = False
            break
    if flag:
        sys.exit('No usable variable')

    percent = np.zeros(indata.variables[var][:][0].shape)
    for i in range(len(time)):
        #print i
        var_slice = indata.variables[var][:][i]
        #percent[np.where(np.logical_not(var_slice.mask))] += 1.
        percent[np.where(var_slice.data>0.)] += 1.
    percent = percent/float(len(time))*100.

    try:
        fill_value = indata.variables[var]._FillValue
    except:
        fill_value = -2.e9
    percent[np.where(percent<=0.0)] = fill_value
    try:
        percent_var = outdata.createVariable('percent_of_time_covered_by_ice',
                                             np.float64,
                                             ('y', 'x',),
                                             fill_value=fill_value)
    except:
        percent_var = outdata.variables['percent_of_time_covered_by_ice']
    percent_var[:] = percent

    indata.close()
    if not overwrite:
        outdata.close()
        
def calc_time_averaged(infile=None, outfile=None, var_name=None, *args, **kwargs):
    fill_value = -2e9
    if var_name is None:
        sys.exit('Must provide a var name!')
    if infile is None:
        sys.exit('Must provide an input file!')
    overwrite = False
    if outfile is None:
        outfile = infile
        overwrite = True
    else:
        nc_copy_dim(infile, outfile)

    if overwrite:
        indata = Dataset(infile, 'a')
        outdata = indata
    else:
        indata = Dataset(infile, 'r')
        outdata = Dataset(outfile, 'a')

    time = indata.variables['time'][:]/(365*24*3600.)
    '''
    for start_index in range(len(time)):
        if time[start_index]>=time[0]+unstable_time:
            break
    '''
    start_index = 0
    end_index = len(time)-1

    var = indata.variables[var_name][:]
    if np.ma.is_masked(var):
        var = var.data
        var[np.where(var == fill_value)] = 0.0
    var_time_averaged = np.zeros(var[0].shape)
    for i in range(start_index, end_index+1):
        tmp_var_slice = var[i]
        if i==start_index:
            time_step = (time[i+1]-time[i])/2
        elif i==end_index:
            time_step = (time[i]-time[i-1])/2
        else:
            time_step = (time[i+1]-time[i-1])/2
        var_time_averaged = var_time_averaged+tmp_var_slice*time_step
    var_time_averaged = var_time_averaged/(time[end_index]-time[start_index])
    var_time_averaged[np.where(var_time_averaged<=0.)] = fill_value
    var_time_averaged = np.ma.masked_values(var_time_averaged,
                                            fill_value)
    try:
        var_time_averaged_var = outdata.createVariable(
            '{}_time_averaged'.format(var_name), np.float64, ('y','x',),
            fill_value=fill_value)
    except:
        var_time_averaged_var = outdata.variables['{}_time_averaged'.format(var_name)]
        var_time_averaged_var._FillValue = fill_value
    var_time_averaged_var[:] = var_time_averaged
    try:
        var_time_averaged_var.units = indata.variables[var_name].units
    except:
        pass

    indata.close()
    if not overwrite:
        outdata.close()

def calc_time_percent(infile=None, outfile=None, var_name=None, *args, **kwargs):
    # Percent% of the time the value of var is lower than the return value

    fill_value = -2e9

    if var_name is None:
        sys.exit('Must provide a var name!')
    try:
        percent = kwargs['percent']
    except:
        sys.exit('Must provide a percent!')
    
    if infile is None:
        sys.exit('Must provide an input file!')
    overwrite = False
    if outfile is None or outfile==infile:
        outfile = infile
        overwrite = True
    else:
        nc_copy_dim(infile, outfile)

    if overwrite:
        indata = Dataset(infile, 'a')
        outdata = indata
    else:
        indata = Dataset(infile, 'r')
        outdata = Dataset(outfile, 'a')

    time = indata.variables['time'][:]/(365*24*3600.)

    var_sorted = indata.variables[var_name][:].copy()
    var_sorted.sort(axis=0)

    ncols, nrows = var_sorted[0].shape
    var_percent_value = np.zeros((ncols, nrows))
    for i in range(ncols):
        for j in range(nrows):
            try:
                n_time = var_sorted[:, i, j].count()
            except:
                n_time = len(var_sorted[:, i, j])
            if n_time==0:
                var_percent_value[i, j] = fill_value
            else:
                index = int(percent/100.0*n_time-1)
                var_percent_value[i, j] = var_sorted[:, i, j][index]
    var_percent_value[np.where(var_percent_value<=0)] = fill_value
    var_percent_value = np.ma.masked_values(var_percent_value, fill_value)

    try:
        var_percent = outdata.createVariable(
            '{}_time_{}_percent'.format(var_name, percent), np.float64, ('y','x',),
            fill_value=fill_value)
    except:
        var_percent = outdata.variables['{}_time_{}_percent'.format(var_name, percent)]
    var_percent[:] = var_percent_value
    var_percent.units = indata.variables[var_name].units

    indata.close()
    if not overwrite:
        outdata.close()
    
