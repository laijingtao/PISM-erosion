#!/usr/bin/env python

# nc2gmt.py
# Auther: Jingtao Lai, Geology@UIUC
# Description: convert input to output that is ready for GMT


import sys
import numpy as np
from argparse import ArgumentParser

# try different netCDF modules
try:
    from netCDF4 import Dataset
except:
    print "netCDF4 is not installed!"
    sys.exit(1)


# Set up the option parser
parser = ArgumentParser()
parser.description = '''Extract the lon, lat and z inforation from the \
                        original file, and write to a xyz file.'''
parser.add_argument("-i", "--input", dest="infile",
                    help="Input file.", default=None)
parser.add_argument("-o", "--output", dest="outfile",
                    help="Output file.", default="out.xyz")
parser.add_argument("-v", "--var", dest="var",
                    help="Name of variable.", default=None)

options = parser.parse_args()
infile = options.infile
outfile = options.outfile
var = options.var

if infile is None:
    print('Input file required.')
    parser.print_help()
    exit(0)

if var is None:
    print('Variable name required.')
    parser.print_help()
    exit(0)

indata = Dataset(infile, 'r')

if not(var in indata.variables.keys()):
    print('Didn\'t find variable \''+var+'\'.')
    exit(0)
x = indata.variables['lon'][:]
y = indata.variables['lat'][:]
z = indata.variables[var][:]
z = np.squeeze(z)
z = z.astype(float)

indata.close()

z[np.where(np.isnan(z))] = 0.0
z[np.where(z<=0.)] = np.nan

with open(outfile, 'w') as f:
   n, m = x.shape
   for i in range(n):
       for j in range(m):
           f.write('{} {} {}\n'.format(x[i, j], y[i, j], z[i, j]))


