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

def main(*args, **kwargs):
    try:
        infile = kwargs['infile']
    except:
        infile = 'test.nc'
    try:
        outfile = kwargs['outfile']
    except:
        outfile = 'out.jpg'
    try:
        outpath = kwargs['outpath']
    except:
        outpath = '/mnt/c/Users/jtlai/Work/Research/glacier/plot/'
    try:
        var = kwargs['var']
    except:
        var = 'topg'
    try:
        region = kwargs['region']
    except:
        #region = '-124.8/-122.6/47/48.4' 
        region = '-124.5/-122.75/47.25/48.25'
    try:
        cbartitle = kwargs['cbartitle']
        if cbartitle is None:
            cbartitle = var
    except:
        cbartitle = var
    try:
        cbarunit = kwargs['cbarunit']
        if cbarunit is None:
            cbarunit = ' '
    except:
        cbarunit = ' '
    try:
        cbarrange = kwargs['cbarrange']
    except:
        cbarrange = 'auto'
        #cbarrange = '0/2200/100'
   
    if infile is None:
        sys.exit('Input file needed.')
    if var is None:
        sys.exit('Variable name needed.')
    if not os.path.isdir(outpath):
        sys.exit('Invalid output path!')

    if cbarrange == 'auto':
        from netCDF4 import Dataset
        indata = Dataset(infile, 'r')
        z_in = indata.variables[var][:]
        z_in = np.squeeze(z_in)
        z_in = z_in.astype(float)
        cbarrange = str(z_in.min())+'/'+str(z_in.max())+'/'\
            +str((z_in.max()-z_in.min())/10.0)
        indata.close()

    # convert nc data to xyz data
    print 'Converting data ...'
    nc2gmt.main(infile=infile, outfile='tmp.xyz', var=var,
                srs='+init=epsg:26710', interp_flag='True')

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

if __name__=='__main__':
    # Set up the option parser
    parser = ArgumentParser()
    parser.description = '''pythonized GMT plotting'''
    parser.add_argument('-i', '--input', dest='infile',
                        help='Input file.', default=None)
    parser.add_argument('-o', '--output', dest='outfile',
                        help='Output file.', default='out.jpg')
    parser.add_argument('--o_dir', dest='outpath',
                        help='Output directory.',
                        default='/mnt/c/Users/jtlai/Work/Research/glacier/plot/')
    parser.add_argument('-v', '--var', dest='var',
                        default=None)
    parser.add_argument('-r', '--region', dest='region',
                        default='-124.5/-122.75/47.25/48.25')
    parser.add_argument('--cbartitle', dest='cbartitle',
                        default=None)
    parser.add_argument('--cbarunit', dest='cbarunit',
                        default=None)
    parser.add_argument('--cbarrange', dest='cbarrange',
                        default='auto')
    
    options = parser.parse_args()
    infile = options.infile
    outfile = options.outfile
    outpath = options.outpath
    var = options.var
    region = options.region
    cbartitle = options.cbartitle
    cbarunit = options.cbarunit
    cbarrange = options.cbarrange

    main(infile=infile, outfile=outfile, outpath=outpath, var=var,
         region=region, cbartitle=cbartitle, cbarunit=cbarunit, 
         cbarrange=cbarrange)
