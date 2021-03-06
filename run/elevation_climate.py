#!/usr/bin/env python
# Copyright (C) 2016 Andy Aschwanden
# Modified by Jingtao Lai. 2017 Apr.

import numpy as np
import itertools
from collections import OrderedDict
import os
try:
    import subprocess32 as sub
except:
    import subprocess as sub
from argparse import ArgumentParser
import sys
sys.path.append('../resources/')
from resources import *


# set up the option parser
parser = ArgumentParser()
parser.description = "Generating scripts for prognostic simulations."
parser.add_argument("-n", '--n_procs', dest="n", type=int,
                    help='''number of cores/processors. default=2.''', default=2)
parser.add_argument("-w", '--wall_time', dest="walltime",
                    help='''walltime. default: 24:00:00.''', default="24:00:00")
parser.add_argument("-q", '--queue', dest="queue", choices=list_queues(),
                    help='''queue. default=t1standard.''', default='normal')
parser.add_argument("--climate", dest="climate",
                    choices=['elev', 'paleo', 'present'],
                    help="Climate", default='elev')
parser.add_argument("-d", "--domain", dest="domain",
                    choices=['olympics', 'olympics_mtns', 'synthetic'],
                    help="sets the modeling domain", default='synthetic')
parser.add_argument("--start_year", dest="start", type=float,
                    help="Start year", default=0)
parser.add_argument("--duration", dest="duration", type=float,
                    help="Duration", default=10)
parser.add_argument("--exstep", dest="exstep", type=float,
                    help="Spatial time series writing interval", default=100)
parser.add_argument("-f", "--o_format", dest="oformat",
                    choices=['netcdf3', 'netcdf4_parallel', 'pnetcdf'],
                    help="output format", default='netcdf3')
parser.add_argument("-g", "--grid", dest="grid", type=int,
                    choices=accepted_resolutions(),
                    help="horizontal grid resolution", default=1000)
parser.add_argument("-i", "--input_file", dest="input_file",
                    help="Input file to restart from", default=None)
parser.add_argument("--bootstrap", dest="bootstrap",
                    help="Turn on/off bootstrap when -i is provided", default=False)
parser.add_argument("--o_dir", dest="odir",
                    help="output directory. Default: current directory", default='test')
parser.add_argument("--o_size", dest="osize",
                    choices=['small', 'medium', 'big', 'big_2d'],
                    help="output size type", default='medium')
parser.add_argument("-s", "--system", dest="system",
                    choices=list_systems(),
                    help="computer system to use.", default='debug')
parser.add_argument("--bed_version", dest="bed_version",
                    help="Version of bed DEM.", default='1')
parser.add_argument("--stress_balance", dest="stress_balance",
                    choices=['sia', 'ssa+sia', 'ssa'],
                    help="stress balance solver", default='sia')


options = parser.parse_args()


nn = options.n
odir = options.odir
oformat = options.oformat
osize = options.osize
queue = options.queue
walltime = options.walltime
system = options.system

if system in ['keeling']:
    pism_data_dir = os.environ['PISM_DATA_DIR']

bed_version = options.bed_version
climate = options.climate
duration = options.duration
grid = options.grid
stress_balance = options.stress_balance

start = options.start
end  = start + options.duration
exstep = options.exstep
domain = options.domain
pism_exec = generate_domain(domain)
input_file = options.input_file
bootstrap = options.bootstrap
if input_file is None:
    bootstrap = True
pism_dataname = 'pism_{domain}_{grid}m_v{version}.nc'.format(domain=domain.capitalize(),
                                                             grid=grid,
                                                             version=bed_version)
if system in ['keeling']:
    pism_dataname = pism_data_dir+'/bed_dem/'+pism_dataname
state_dir = 'state'
scalar_dir = 'scalar'
spatial_dir = 'spatial'
if system in ['keeling']:
    odir = pism_data_dir+odir
if not os.path.isdir(odir):
    os.mkdir(odir)
for tsdir in (scalar_dir, spatial_dir, state_dir):
    if not os.path.isdir(os.path.join(odir, tsdir)):
        os.mkdir(os.path.join(odir, tsdir))
odir_tmp = '_'.join([odir, 'tmp'])
if not os.path.isdir(odir_tmp):
    os.mkdir(odir_tmp)

# Configuration File Setup
pism_config = 'olympics_config'
pism_config_nc = '.'.join([pism_config, 'nc'])
if system in ['keeling']:
    pism_config_nc = pism_data_dir+pism_config_nc
pism_config_cdl = os.path.join('../config', '.'.join([pism_config, 'cdl']))
# Anaconda libssl problem
if system in ['chinook', 'keeling', 'debug']:
    ncgen = '/usr/bin/ncgen'
else:
    ncgen = 'ncgen'
cmd = [ncgen, '-o',
       pism_config_nc, pism_config_cdl]
sub.call(cmd)


hydrology = 'diffuse'

# ########################################################
# set up model initialization
# ########################################################


# Model Parameters
ssa_n = (3.0)
ssa_e = (1.0)

# Model Parameters for Sensitivity Study
ela_values = [2200, 2000, 1800]
mb_min_values = [-3.]
mb_max_values = [3.]
sia_e_values = [3.0]
ppq_values = [0.50]
tefo_values = [0.020]
phi_min_values = [15]
phi_max_values = [45]
topg_min_values = [0]
topg_max_values = [3000]
temp_lapse_rate_values = [5.0, 6.0]
combinations = list(itertools.product(ela_values,
                                      mb_min_values,
                                      mb_max_values,
                                      sia_e_values,
                                      ppq_values,
                                      tefo_values,
                                      phi_min_values,
                                      phi_max_values,
                                      topg_min_values,
                                      topg_max_values))

tsstep = 'yearly'

scripts = []
scripts_post = []
batch_scripts_dir = './'
if system in ['keeling']:
    batch_scripts_dir = './batch_scripts/'
if not os.path.isdir(batch_scripts_dir):
    os.mkdir(batch_scripts_dir)

for n, combination in enumerate(combinations):

    ela, mb_min, mb_max, sia_e, ppq, tefo, phi_min, phi_max, topg_min, topg_max = combination

    ttphi = '{},{},{},{}'.format(phi_min, phi_max, topg_min, topg_max)

    name_options = OrderedDict()
    name_options['sb'] = stress_balance
    name_options['ela'] = ela
    name_options['mb_min'] = mb_min
    name_options['mb_max'] = mb_max
    experiment =  '_'.join([climate, '_'.join(['_'.join([k, str(v)]) for k, v in name_options.items()])])

    script = 'cc_{}_g{}m_{}.sh'.format(domain.lower(), grid, experiment)
    scripts.append(script)
    script_post = 'cc_{}_g{}m_{}_post.sh'.format(domain.lower(), grid, experiment)
    scripts_post.append(script_post)

    for filename in (script):
        try:
            os.remove(filename)
        except OSError:
            pass

    batch_header, batch_system = make_batch_header(system, nn, walltime, queue)

    with open(batch_scripts_dir+script, 'w') as f:

        f.write(batch_header)
        
        outfile = '{domain}_g{grid}m_{experiment}_{start}_{end}a.nc'.format(domain=domain.lower(),
                                                                            grid=grid,
                                                                            experiment=experiment,
                                                                            start=int(start),
                                                                            end=int(end))

        prefix = generate_prefix_str(pism_exec)

        # Setup General Parameters
        general_params_dict = OrderedDict()
        if input_file is None:
            general_params_dict['i'] = pism_dataname
            general_params_dict['bootstrap'] = ''
        else:
            general_params_dict['i'] = input_file
            if bootstrap:
                general_params_dict['bootstrap'] = ''
        general_params_dict['ys'] = start
        general_params_dict['ye'] = end
        general_params_dict['o'] = os.path.join(odir, state_dir, outfile)
        general_params_dict['o_format'] = oformat
        general_params_dict['o_size'] = osize
        general_params_dict['config_override'] = pism_config_nc

        if input_file is None:
            grid_params_dict = generate_grid_description(grid, accepted_resolutions(), domain)
        else:
            if bootstrap:
                grid_params_dict = generate_grid_description(grid, accepted_resolutions(), 
                                                             domain, restart=False)
            else:
                grid_params_dict = generate_grid_description(grid, accepted_resolutions(), 
                                                             domain, restart=True)

        # Setup Stress Balance Paramters
        sb_params_dict = OrderedDict()
        sb_params_dict['sia_e'] = sia_e
        sb_params_dict['ssa_e'] = ssa_e
        sb_params_dict['ssa_n'] = ssa_n
        sb_params_dict['pseudo_plastic_q'] = ppq
        sb_params_dict['till_effective_fraction_overburden'] = tefo
        sb_params_dict['topg_to_phi'] = ttphi
        sb_params_dict['ssa_method'] = 'fd'

        stress_balance_params_dict = generate_stress_balance(stress_balance, sb_params_dict)

        # Setup Climate Forcing
        climate_params_dict = generate_climate(climate,
                                               ice_surface_temp='0,0,-100,5000',
                                               climatic_mass_balance='{},{},0,{},2500'.format(mb_min, mb_max, ela))

        # Setup Ocean Forcing
        ocean_params_dict = generate_ocean('null')

        # Setup Hydrology Model
        hydro_params_dict = generate_hydrology(hydrology)

        # Setup Carving Model
        calving_params_dict = generate_calving('float_kill')


        # Setup Scalar and Spatial Time Series Reporting
        exvars = default_spatial_ts_vars()
        spatial_ts_dict = generate_spatial_ts(outfile, exvars, exstep, odir=odir_tmp, split=True)
        scalar_ts_dict = generate_scalar_ts(outfile, tsstep, odir=os.path.join(odir, scalar_dir))

        # Merge All Parameter Dictionaries
        all_params_dict = merge_dicts(general_params_dict, 
                                      grid_params_dict, 
                                      stress_balance_params_dict, 
                                      climate_params_dict, 
                                      ocean_params_dict, 
                                      hydro_params_dict, 
                                      calving_params_dict, 
                                      spatial_ts_dict, 
                                      scalar_ts_dict)
        all_params = ' '.join([' '.join(['-' + k, str(v)]) for k, v in all_params_dict.items()])

        if system in ('debug'):
            cmd = ' '.join([batch_system['mpido'], prefix, all_params, 
                           '2>&1 | tee {outdir}/job.${batch}'.format(outdir=odir, batch=batch_system['job_id'])])
        else:
            cmd = ' '.join([batch_system['mpido'], prefix, all_params, 
                           '> {outdir}/job.${batch}  2>&1'.format(outdir=odir, batch=batch_system['job_id'])])

        f.write(cmd)
        f.write('\n')
        f.write('\n')
        f.write('echo \"Done\"\n')
        f.write('\n')
        f.write('bash {}\n'.format(script_post))
        f.write('\n')

    post_header = make_batch_post_header(system)

    with open(batch_scripts_dir+script_post, 'w') as f:

        f.write(post_header)

        extra_file = spatial_ts_dict['extra_file']
        #myfiles = ' '.join(['{}_{:.3f}.nc'.format(extra_file, k) for k in np.arange(start + exstep, end, exstep)])
        myfiles = ' '.join(['{}-{:012.3f}.nc'.format(extra_file, k) \
                            for k in np.arange(start + exstep, end, exstep)])
        myoutfile = extra_file + '.nc'
        myoutfile = os.path.join(odir, spatial_dir, os.path.split(myoutfile)[-1])
        cmd = ' '.join(['ncrcat -O -6 -h', myfiles, myoutfile, '\n'])
        f.write(cmd)
        cmd = ' '.join(['ncks -O -4', os.path.join(odir, state_dir, outfile),
                        os.path.join(odir, state_dir, outfile), '\n'])
        f.write(cmd)
        # pism0.7 has x,y dimensions, we want y,x
        cmd = ' '.join(['ncpdq', '-a', 'y,x', '-O', myoutfile, myoutfile, '\n'])
        f.write(cmd)
        cmd = ' '.join(['ncpdq', '-a', 'y,x', '-O', 
                        os.path.join(odir, state_dir, outfile),
                        os.path.join(odir, state_dir, outfile), '\n'])
        f.write(cmd)


scripts = uniquify_list(scripts)
scripts_post = uniquify_list(scripts_post)
print '\n'.join([script for script in scripts])
print('\nwritten\n')
print '\n'.join([script for script in scripts_post])
print('\nwritten\n')

if system in ['keeling']:
    with open('submit_batch.sh', 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write('cd ./batch_scripts\n')
        for script in scripts:
            f.write('sbatch {}\n'.format(script))
    sub.call(['chmod', 'u+x', 'submit_batch.sh'])

