#!/usr/bin/env python
# Copyright (C) 2017 Jingtao Lai

import os
import sys
try:
    import subprocess32 as sub
except:
    import subprocess as sub
from argparse import ArgumentParser
import numpy as np
sys.path.append('./tools')
import nc2gmt

inpath = './'
outpath = '/mnt/d/jlai11/pism-olympics/plot'

if not os.path.isdir(inpath):
    sys.exit('Invalid input path!')
if not os.path.isdir(outpath):
    sys.exit('Invalid output path!')

infile = 'test.nc'
outfile = 'test.jpg'
var = 'topg'
#region = '-124.8/-122.6/47/48.4' 
region = '-124.5/-122.75/47.25/48.25'
cbartitle = 'Elevation'
cbarunit = 'm'
cbarrange = 'Auto'
#cbarrange = '0/2200/100'

if cbarrange == 'Auto':
    from netCDF4 import Dataset
    indata = Dataset(infile, 'r')
    z_in = indata.variables[var][:]
    z_in = np.squeeze(z_in)
    z_in = z_in.astype(float)
    cbarrange = str(int(z_in.min()))+'/'+str(int(z_in.max()))+'/'\
        +str((z_in.max()-z_in.min())/10.0)

# convert nc data to xyz data
print 'Converting data ...'
nc2gmt.main(infile=infile, outfile='tmp.xyz', var=var,
            srs='+init=epsg:26710', interp_flat='True')
#cmd = ['python', './tools/nc2gmt.py', '-i', infile, '-o', 'tmp.xyz',
#        '-v', var, '--srs', '+init=epsg:26710', '--interp', 'True']
#sub.call(cmd)

with open('tmp.gmt', 'w') as f:
    cmd = ['gmt', 'xyz2grd', 'tmp.xyz', '-Gtmp.nc', '-R'+region, '-I0.5m']
    f.write(sub.list2cmdline(cmd)+'\n')
    cmd = ['gmt', 'psbasemap', '-R'+region, '-Jm6', '-Ba0.5f0.25', '-V', '-P',
           '-K', '-X1.5', '-Y2', '>>', 'tmp.ps']
    f.write(sub.list2cmdline(cmd)+'\n')
    cmd = ['gmt', 'makecpt', '-T'+cbarrange, '-Cjet', '-Z', '>tmp.cpt']
    f.write(sub.list2cmdline(cmd)+'\n')
    cmd = ['gmt', 'grdimage', 'tmp.nc', '-Ctmp.cpt', '-Jm', '-E300', '-nb',
           '-Q', '-P', '-O', '-K', '>>', 'tmp.ps']
    f.write(sub.list2cmdline(cmd)+'\n')
    cmd = ['gmt', 'psscale', '-Dx12.5c/0.11c+w8c/0.4c', '-Ctmp.cpt', '-Baf',
           '-Bx+l'+cbartitle, '-By+l'+cbarunit, '-O', '>>', 'tmp.ps']
    f.write(sub.list2cmdline(cmd)+'\n')
    cmd = ['gmt', 'psconvert', 'tmp.ps', '-A', '-Tj']
    f.write(sub.list2cmdline(cmd)+'\n')

cmd = ['bash', 'tmp.gmt']
sub.call(cmd)
cmd = ['mv', 'tmp.jpg', outfile]
sub.call(cmd)
cmd = ['rm', 'gmt.history', 'tmp.nc', 'tmp.xyz', 'tmp.ps', 'tmp.cpt', 'tmp.gmt']
sub.call(cmd)
print 'GMT finished.'

print 'Moving data to '+outpath
cmd = ['mv', outfile, outpath]
sub.call(cmd)
