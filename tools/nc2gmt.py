#!/usr/bin/env python

# nc2gmt.py
# Auther: Jingtao Lai, Geology@UIUC
# Description: convert input to output that is ready for GMT


import sys
import numpy as np
import scipy.interpolate as interp
from pyproj import Proj
from argparse import ArgumentParser

# try different netCDF modules
try:
    from netCDF4 import Dataset
except:
    print "netCDF4 is not installed!"
    sys.exit(1)

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    if v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        print('Boolean value expected.')
        exit(0)

# Set up the option parser
parser = ArgumentParser()
parser.description = '''Convert netCDF file to xyz file that is ready for GMT'''
parser.add_argument("-i", "--input", dest="infile",
                    help="Input file.", default=None)
parser.add_argument("-o", "--output", dest="outfile",
                    help="Output file.", default="out.nc")
parser.add_argument("-v", "--var", dest="var",
                    help="Name of variable.", default=None)
parser.add_argument("--srs", dest="srs",
                    help='''a valid proj4 string describing describing the projection''',
                    default=None)
parser.add_argument("--interp", dest="interp_flag",
                    help='''Turn on/off interpolation''',
                    default=False)
parser.add_argument("--interp_factor", dest="interp_factor",
                    help='''Interpolation factor''',
                    default=5)

options = parser.parse_args()
infile = options.infile
outfile = options.outfile
var = options.var
srs = options.srs
interp_flag = str2bool(options.interp_flag)
interp_factor = options.interp_factor

if infile is None:
    print('Input file required.')
    parser.print_help()
    exit(0)

if var is None:
    print('Variable name required.')
    parser.print_help()
    exit(0)


def get_projection_from_file(nc):

    # First, check if we have a global attribute 'proj4'
    # which contains a Proj4 string:
    try:
        p = Proj(str(nc.proj4))
        print(
            'Found projection information in global attribute proj4, using it')
    except:
        try:
            p = Proj(str(nc.projection))
            print(
                'Found projection information in global attribute projection, using it')
        except:
            try:
                # go through variables and look for 'grid_mapping' attribute
                for var in nc.variables.keys():
                    if hasattr(nc.variables[var], 'grid_mapping'):
                        mappingvarname = nc.variables[var].grid_mapping
                        print(
                            'Found projection information in variable "%s", using it' % mappingvarname)
                        break
                var_mapping = nc.variables[mappingvarname]
                p = Proj(proj="stere",
                         ellps=var_mapping.ellipsoid,
                         datum=var_mapping.ellipsoid,
                         units="m",
                         lat_ts=var_mapping.standard_parallel,
                         lat_0=var_mapping.latitude_of_projection_origin,
                         lon_0=var_mapping.straight_vertical_longitude_from_pole,
                         x_0=var_mapping.false_easting,
                         y_0=var_mapping.false_northing)
            except:
                print('No mapping information found, exiting.')
                sys.exit(1)

    return p


indata = Dataset(infile, 'r')

if not(var in indata.variables.keys()):
    print('Didn\'t find variable \''+var+'\'.')
    exit(0)

x_in = indata.variables['x'][:]
y_in = indata.variables['y'][:]
z_in = indata.variables[var][:]

if srs:
    # use projection from command line
    try:
        proj = Proj(init=srs)
    except:
        proj = Proj(srs)
else:
    # Get projection from file
    proj = get_projection_from_file(indata)

indata.close()

z_in = np.squeeze(z_in)

if interp_flag:
    print('nc2gmt.py: Interpolation is on. Interpolation factor is '+str(interp_factor))
    nx = len(x_in)
    ny = len(y_in)

    easting = np.zeros(nx*interp_factor)
    for i in range(nx-1):
        index = np.linspace(i*interp_factor, (i+1)*interp_factor, interp_factor+1, dtype='int')
        tmp = np.linspace(x_in[i], x_in[i+1], interp_factor+1, dtype='int')
        easting[index[0:len(index)-1]] = tmp[0:len(index)-1]
    easting[nx*interp_factor-1] = x_in[nx-1]

    northing = np.zeros(ny*interp_factor)
    for i in range(ny-1):
        index = np.linspace(i*interp_factor, (i+1)*interp_factor, interp_factor+1, dtype='int')
        tmp = np.linspace(y_in[i], y_in[i+1], interp_factor+1, dtype='int')
        northing[index[0:len(index)-1]] = tmp[0:len(index)-1]
    northing[ny*interp_factor-1] = y_in[ny-1]

    nn, ee = np.meshgrid(northing, easting)
    yy, xx = np.meshgrid(y_in, x_in)
    xx = np.reshape(xx, nx*ny)
    yy = np.reshape(yy, nx*ny)
    z_in = np.reshape(z_in, nx*ny)
    z_out = interp.griddata((xx, yy), z_in, (ee, nn), method='linear')
else:
    easting = x_in
    northing = y_in
    z_out = z_in
    nn, ee = np.meshgrid(northing, easting)

lon, lat = proj(ee, nn, inverse=True)

z_out[np.where(z_out==0)] = np.nan

with open(outfile, 'w') as f:
   n, m = z_out.shape
   for i in range(n):
       for j in range(m):
           f.write('{} {} {}\n'.format(lon[i, j], lat[i, j], z_out[i, j]))
