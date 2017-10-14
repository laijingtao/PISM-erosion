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
#sys.path.append('../../resources/')
#from resources import *

def generate_file_names_synthetic():
    domain = 'synthetic'
    grid = '1000'
    start = -125000
    end = 0
    stress_balance = 'ssa+sia'
    climate = 'paleo'
    delta_T_values = [1.5, 0.0, -1.5]
    frac_P_values = [0.5, 1.0, 1.5]
    combinations = list(itertools.product(delta_T_values, frac_P_values))

    file_name_list = []
    for n, combination in enumerate(combinations):
        delta_T, frac_P = combination

        name_options = OrderedDict()
        name_options['sb'] = stress_balance
        name_options['delta_T'] = delta_T
        name_options['frac_P'] = frac_P
        experiment =  '_'.join([climate, '_'.join(['_'.join([k, str(v)]) for k, v
            in name_options.items()])])
        file_name =\
        '{domain}_g{grid}m_{experiment}_{start}_{end}a.nc'.format(domain=domain.lower(),
                                                                  grid=grid,
                                                                  experiment=experiment,
                                                                  start=int(start),
                                                                  end=int(end))
        file_name_list.append(file_name)
        
    return file_name_list


def generate_file_names():
    domain = 'olympics'
    grid = '1000'
    start = -125000
    end = 0
    stress_balance = 'ssa+sia'
    climate = 'paleo'

    sia_e_values = [1.0, 3.0]
    q_values = [0.33, 0.5]
    phi_values = [20, 30]
    combinations = list(itertools.product(sia_e_values, q_values,
                                          phi_values))

    file_name_list = []
    for n, combination in enumerate(combinations):
        sia_e, q, phi = combination

        name_options = OrderedDict()
        name_options['sb'] = stress_balance
        name_options['sia_e'] = sia_e
        name_options['q'] = q
        name_options['phi'] = phi
        experiment =  '_'.join([climate, '_'.join(['_'.join([k, str(v)]) for k, v
            in name_options.items()])])
        file_name =\
        '{domain}_g{grid}m_{experiment}_{start}_{end}a.nc'.format(domain=domain.lower(),
                                                                  grid=grid,
                                                                  experiment=experiment,
                                                                  start=int(start),
                                                                  end=int(end))
        file_name_list.append(file_name)
        
    return file_name_list

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

def nc_regrid_dem(infile=None, outfile=None, grid=None):
    # Based on gdalwarp and nco
    # if infile's dimension has units, gdalwarp won't work
    cmd = ['ncatted', '-a', 'units,x,d,,', '-a', 'units,y,d,,',
           infile, 'tmp.nc']
    sub.call(cmd)
    cmd = ['gdalwarp', '-of', 'netCDF', '-overwrite', 
           '-tr', str(grid), str(grid), 'tmp.nc', outfile]
    print ' '.join(cmd)
    sub.call(cmd)
    cmd = ['ncrename', '-v', 'Band1,topg', '-v', 'lon,x', '-v', 'lat,y', outfile]
    print ' '.join(cmd)
    sub.call(cmd)
    cmd = ['ncrename', '-d', 'lon,x', '-d', 'lat,y', outfile]
    print ' '.join(cmd)
    sub.call(cmd)
    cmd = ['ncatted', '-a', 'units,x,o,c,m', '-a', 'units,y,o,c,m',
           '-a', 'standard_name,x,d,,', '-a', 'long_name,x,d,,', 
           '-a', 'standard_name,y,d,,', '-a', 'long_name,y,d,,', 
           '-a', 'long_name,topg,d,,', '-a', '_FillValue,topg,o,f,-2.0e9',
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
    if infile is None:
        sys.exit('Must provide an input file!')
    overwrite = False
    if outfile is None:
        outfile = infile
        overwrite = True
    else:
        nc_copy_dim(infile, outfile)

    indata = Dataset(infile, 'a')
    if overwrite:
        outdata = indata
    else:
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
        try:
            erosion_var = outdata.createVariable(
                erosion_name, np.float64, ('time', 'y', 'x',),
                fill_value=-2000000000.0)
        except:
            erosion_var = outdata.variables[erosion_name]
        erosion_var[:] = erosion
        erosion_var.units = 'm year-1'

    indata.close()
    if not overwrite:
        outdata.close()


def calc_total_erosion(infile=None, outfile=None):
    if infile is None:
        sys.exit('Must provide an input file!')
    overwrite = False
    if outfile is None:
        outfile = infile
        overwrite = True
    else:
        nc_copy_dim(infile, outfile)

    indata = Dataset(infile, 'a')
    if overwrite:
        outdata = indata
    else:
        outdata = Dataset(outfile, 'a')

    time = indata.variables['time'][:]
    grid = abs(indata.variables['x'][0]-indata.variables['x'][1])
    erosion_1 = indata.variables['erosion_1'][:]
    erosion_2 = indata.variables['erosion_2'][:]
    
    total_erosion_1 = np.array(
        [erosion_slice.sum() for erosion_slice in erosion_1])*grid*grid
    total_erosion_2 = np.array(
        [erosion_slice.sum() for erosion_slice in erosion_2])*grid*grid

    total_erosion_1_var = outdata.createVariable('total_erosion_1', 
        np.float64, ('time',))
    total_erosion_1_var[:] = total_erosion_1
    total_erosion_1_var.units = 'm3 year-1'
    total_erosion_2_var = outdata.createVariable('total_erosion_2', 
        np.float64, ('time',))
    total_erosion_2_var[:] = total_erosion_2
    total_erosion_2_var.units = 'm3 year-1'
    
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

    indata = Dataset(infile, 'a')
    if overwrite:
        outdata = indata
    else:
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

def calc_erosion_time_averaged(infile=None, outfile=None):
    if infile is None:
        sys.exit('Must provide an input file!')
    overwrite = False
    if outfile is None:
        outfile = infile
        overwrite = True
    else:
        nc_copy_dim(infile, outfile)

    indata = Dataset(infile, 'a')
    if overwrite:
        outdata = indata
    else:
        outdata = Dataset(outfile, 'a')

    time = indata.variables['time'][:]/(365*24*3600.)

    for erosion_name in ['erosion_1', 'erosion_2']:
        erosion = indata.variables[erosion_name][:]
        erosion_time_averaged = np.zeros(erosion[0].shape)
        for i in range(len(erosion)):
            tmp_erosion_slice = erosion[i].data
            tmp_erosion_slice[np.where(erosion[i].mask)] = 0.
            if i==0:
                time_step = (time[i+1]-time[i])/2
            elif i==len(time)-1:
                time_step = (time[i]-time[i-1])/2
            else:
                time_step = (time[i+1]-time[i-1])/2
            erosion_time_averaged = erosion_time_averaged+tmp_erosion_slice*time_step
        erosion_time_averaged = erosion_time_averaged/(time[-1]-time[0])
        erosion_time_averaged[np.where(erosion_time_averaged<=0.)] =\
            indata.variables[erosion_name]._FillValue
        try:
            erosion_time_averaged_var = outdata.createVariable(
                erosion_name+'_time_averaged', np.float64, ('y','x',),
                fill_value=indata.variables[erosion_name]._FillValue)
        except:
            erosion_time_averaged_var = outdata.variables[erosion_name+'_time_averaged']
        erosion_time_averaged_var[:] = erosion_time_averaged
        erosion_time_averaged_var.units = 'm year-1'

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

    indata = Dataset(infile, 'a')
    if overwrite:
        outdata = indata
    else:
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

