#!/usr/bin/env python
# Copyright (C) 2017 Jingtao Lai

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import itertools
try:
    import subprocess32 as sub
except:
    import subprocess as sub
from collections import OrderedDict
from argparse import ArgumentParser
#sys.path.append('../../resources/')
#from resources import *


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
    combinations = list(itertools.product(sia_e_values,
                                          q_values,
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

def get_grid_size(infile):
    from netCDF4 import Dataset

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

def calc_total_erosion(infile=None, outfile=None):
    from netCDF4 import Dataset

    if infile is None:
        sys.exit('Must provide an input file!')
    overwrite = False
    if outfile is None:
        outfile = infile
        overwrite = True

    indata = Dataset(infile, 'a')

    time = indata.variables['time'][:]
    grid = abs(indata.variables['x'][0]-indata.variables['x'][1])
    erosion_1 = indata.variables['erosion_1'][:]
    erosion_2 = indata.variables['erosion_2'][:]
    
    total_erosion_1 = np.array(
        [erosion_slice.sum() for erosion_slice in erosion_1])*grid*grid
    total_erosion_2 = np.array(
        [erosion_slice.sum() for erosion_slice in erosion_2])*grid*grid
    
    if not overwrite:
        outdata = Dataset(outfile, 'w', format='NETCDF4')
        time_dim = outdata.createDimension('time', None)
        time_var = outdata.createVariable('time', np.float64, ('time',))
        time_var[:] = time
        time_var.units = 'seconds since 1-1-1'
        time_var.calendar = '365_day'
        time_var.long_name = 'time'
        time_var.axis = 'T'
    else:
        outdata = indata

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
