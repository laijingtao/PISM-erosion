#!/usr/bin/env python
# Copyright (C) 2017 Jingtao Lai

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
from build_climate_file import build_constant_climate, build_paleo_modifier
from build_calving_file import build_ocean_kill_file
from read_params import read_params


# set up the option parser
parser = ArgumentParser()
parser.description = "Generating scripts for simulations."
parser.add_argument("-n", '--n_procs', dest="n", type=int,
                    help='''number of cores/processors. default=2.''', default=2)
parser.add_argument("-w", '--wall_time', dest="walltime",
                    help='''walltime. default: 24:00:00.''', default="24:00:00")
parser.add_argument("-q", '--queue', dest="queue", choices=list_queues(),
                    help='''queue. default=t1standard.''', default='normal')
parser.add_argument("--climate", dest="climate",
                    choices=['elev', 'paleo', 'present'],
                    help="Climate", default='paleo')
parser.add_argument("-d", "--domain", dest="domain",
                    choices=['olympics', 'olympics_mtns', 'synthetic', 'fluvial'],
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
parser.add_argument("-p", "--params", dest="params_list",
                    help="Comma-separated list with params for sensitivity", default=None)
parser.add_argument("--params_file", dest="params_file",
                    help="Parameters file", default=None)
parser.add_argument("--batch_scripts_dir", dest="batch_scripts_dir",
                    help="directory to save batch scripts", default=None)

options = parser.parse_args()

nn = options.n
odir = options.odir
oformat = options.oformat
osize = options.osize
queue = options.queue
walltime = options.walltime
system = options.system
bed_version = options.bed_version
climate = options.climate
duration = options.duration
grid = options.grid
stress_balance = options.stress_balance

if system in ['keeling']:
    pism_data_dir = os.environ['PISM_DATA_DIR']
    pism_work_dir = '/data/keeling/a/jlai11/glacier/pism-erosion/'
else:
    pism_data_dir = './'
    pism_work_dir = '../'

start = options.start
end  = start + options.duration
#start = 0
#end = 125000/10
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
pism_dataname = pism_data_dir+'/bed_dem/'+pism_dataname
perf_dir = 'performance'
state_dir = 'state'
scalar_dir = 'scalar'
spatial_dir = 'spatial'
postproc_dir = 'postproc'
initdata_dir = 'init_data'
odir = os.path.join(pism_data_dir, odir)
if not os.path.isdir(odir):
    os.mkdir(odir)
for tsdir in [perf_dir, scalar_dir, spatial_dir, state_dir, postproc_dir, initdata_dir]:
    if not os.path.isdir(os.path.join(odir, tsdir)):
        os.mkdir(os.path.join(odir, tsdir))
if input_file is None:
    cmd = ['cp', pism_dataname, 
            os.path.join(odir, initdata_dir, 'initial_condition.nc')]
else:
    cmd = ['cp', input_file, 
            os.path.join(odir, initdata_dir, 'initial_condition.nc')]
sub.call(cmd)
odir_tmp = '_'.join([odir, 'tmp'])
if not os.path.isdir(odir_tmp):
    os.mkdir(odir_tmp)
batch_scripts_dir = options.batch_scripts_dir
if batch_scripts_dir is None:
    if system in ['keeling']:
        batch_scripts_dir = options.odir
    else:
        batch_scripts_dir = './'
if not os.path.isdir(batch_scripts_dir):
    os.mkdir(batch_scripts_dir)

# Configuration File Setup
pism_config = 'olympics_config'
pism_config_nc = '.'.join([pism_config, 'nc'])
pism_config_nc = pism_data_dir+pism_config_nc
pism_config_cdl = os.path.join('../config', '.'.join([pism_config, 'cdl']))
# Anaconda libssl problem
if system in ['keeling', 'debug']:
    ncgen = '/usr/bin/ncgen'
else:
    ncgen = 'ncgen'
cmd = [ncgen, '-o', pism_config_nc, pism_config_cdl]
sub.call(cmd)

# Check if perform sensitivity study
params_list = options.params_list
do_delta_T = False
do_frac_P = False
if params_list is not None:
    params = params_list.split(',')
    if 'delta_T' in params:
        do_delta_T = True
    if 'frac_P' in params:
        do_frac_P = True

params_file = options.params_file

# ########################################################
# set up model initialization
# ########################################################

# Model Parameters
ssa_n = (3.0)
ssa_e = (1.0)

# Model Parameters for Sensitivity Study
wind_direction = 220
precip_scale_factor_values = [0.07]
ela_values = [1500]
mb_min_values = [-5.]
mb_max_values = [5.]
sia_e_values = [3.0]
ppq_values = [0.50]
tefo_values = [0.020]
phi_min_values = [20]
phi_max_values = [30]
topg_min_values = [-500]
topg_max_values = [4000]
temp_lapse_rate_values = [6.0]
if do_delta_T:
    delta_T_values = read_params(params_file, param='delta_T')
else:
    delta_T_values = [0.0]
if do_frac_P:
    frac_P_values = read_params(params_file, param='frac_P')
else:
    frac_P_values = [1.0]
combinations = list(itertools.product(precip_scale_factor_values,
                                      sia_e_values,
                                      ppq_values,
                                      tefo_values,
                                      phi_min_values,
                                      phi_max_values,
                                      topg_min_values,
                                      topg_max_values,
                                      temp_lapse_rate_values,
                                      delta_T_values,
                                      frac_P_values))

tsstep = 'yearly'

# Setup climate file
print pism_dataname
build_constant_climate(
    infile=pism_dataname,
    outfile=os.path.join(odir, initdata_dir, 'constant_climate.nc'),
    air_temp_mean_annual=read_params(params_file, param='air_temp_mean_annual'),
    air_temp_mean_july=read_params(params_file, param='air_temp_mean_july'),
    precipitation=read_params(params_file, param='precipitation'))
build_paleo_modifier(
    delta_T=delta_T_values,
    frac_P=frac_P_values,
    climate_forcing_dir=os.path.join(pism_work_dir, 'data_sets/climate_forcing'),
    out_dir=os.path.join(odir, initdata_dir))

# Setup calving file
ocean_kill_file = os.path.join(odir, initdata_dir, 'ocean_kill_file.nc')
build_ocean_kill_file(infile=pism_dataname,
                      outfile=ocean_kill_file,
                      thk=0)

scripts = []
scripts_post = []
outfile_names = []

# Generate pism script
for n, combination in enumerate(combinations):

    precip_scale_factor, sia_e, ppq, tefo, phi_min, phi_max, topg_min, topg_max, temp_lapse_rate, delta_T, frac_P = combination

    phi_max = phi_min
    ttphi = '{},{},{},{}'.format(phi_min, phi_max, topg_min, topg_max)

    name_options = OrderedDict()
    name_options['sb'] = stress_balance
    name_options['delta_T'] = delta_T
    name_options['frac_P'] = frac_P
    experiment =  '_'.join(
            [climate, '_'.join(['_'.join([k, str(v)]) for k, v in name_options.items()])])

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

    with open(os.path.join(batch_scripts_dir, script), 'w') as f:

        f.write(batch_header)
        
        outfile = '{domain}_g{grid}m_{experiment}_{start}_{end}a.nc'.format(domain=domain.lower(),
                                                                            grid=grid,
                                                                            experiment=experiment,
                                                                            start=int(start),
                                                                            end=int(end))
        outfile_names.append(outfile)

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
            grid_params_dict = generate_grid_description(grid, accepted_resolutions(), domain,
                                                         dem_file=pism_dataname)
        else:
            if bootstrap:
                grid_params_dict = generate_grid_description(grid, accepted_resolutions(), 
                                                             domain, dem_file=pism_dataname,
                                                             restart=False)
            else:
                grid_params_dict = generate_grid_description(grid, accepted_resolutions(), 
                                                             domain, dem_file=pism_dataname,
                                                             restart=True)

        # Setup Stress Balance Paramters
        sb_params_dict = OrderedDict()
        sb_params_dict['sia_e'] = sia_e
        sb_params_dict['ssa_e'] = ssa_e
        sb_params_dict['ssa_n'] = ssa_n
        sb_params_dict['pseudo_plastic_q'] = ppq
        sb_params_dict['till_effective_fraction_overburden'] = tefo
        #sb_params_dict['topg_to_phi'] = ttphi
        sb_params_dict['plastic_phi'] = phi_min
        sb_params_dict['ssafd_ksp_divtol'] = 1e300
        sb_params_dict['cfbc'] = ''
        sb_params_dict['ssa_method'] = 'fd'

        stress_balance_params_dict = generate_stress_balance(stress_balance, sb_params_dict)

        # Setup Climate Forcing
        climate_file = os.path.join(odir, initdata_dir, 'constant_climate.nc')
        atmosphere_paleo_file = os.path.join(
                odir, initdata_dir, 'paleo_modifier_T_{}_P_{}.nc'.format(delta_T, frac_P))
        climate_params_dict = generate_climate(
            climate,
            **{'atmosphere_yearly_cycle_file': climate_file,
               'atmosphere_lapse_rate_file': climate_file,
               'temp_lapse_rate': temp_lapse_rate,
               'atmosphere_delta_T_file': atmosphere_paleo_file,
               #'atmosphere_paleo_precip_file': atmosphere_paleo_file,
               #'atmosphere.precip_exponential_factor_for_temperature': precip_scale_factor,
               'atmosphere_frac_P_file': atmosphere_paleo_file})

        # Setup Ocean Forcing
        ocean_params_dict = generate_ocean('null')

        # Setup Hydrology Model
        hydro_params_dict = generate_hydrology('diffuse')

        # Setup Carving Model
        calving_params_dict = generate_calving('float_kill')
        '''
        calving_params_dict = generate_calving('ocean_kill',
                                               ocean_kill_file=ocean_kill_file)
        '''

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
                           '2>&1 | tee {outdir}/job.${batch}'.format(outdir=odir, 
                                                                     batch=batch_system['job_id'])])
        else:
            cmd = ' '.join([batch_system['mpido'], prefix, all_params, 
                           '> {outdir}/job.${batch}  2>&1'.format(outdir=odir, 
                                                                  batch=batch_system['job_id'])])

        f.write(cmd)
        f.write('\n')
        f.write('\n')
        f.write('echo \"Done\"\n')
        f.write('\n')
        f.write('bash {}\n'.format(script_post))
        f.write('\n')

    post_header = make_batch_post_header(system)

    with open(os.path.join(batch_scripts_dir, script_post), 'w') as f:

        f.write(post_header)

        extra_file = spatial_ts_dict['extra_file']
        myfiles = ' '.join(['{}_{:.3f}.nc'.format(extra_file, k) \
                            for k in np.arange(start + exstep, end, exstep)])
        '''
        myfiles = ' '.join(['{}-{:012.3f}.nc'.format(extra_file, k) \
                            for k in np.arange(start + exstep, end, exstep)])
        '''
        myoutfile = extra_file + '.nc'
        myoutfile = os.path.join(odir, spatial_dir, os.path.split(myoutfile)[-1])
        cmd = ' '.join(['ncrcat -O -6 -h', myfiles, myoutfile, '\n'])
        f.write(cmd)
        cmd = ' '.join(['ncks -O -4', os.path.join(odir, state_dir, outfile),
                        os.path.join(odir, state_dir, outfile), '\n'])
        f.write(cmd)
        # pism0.7 has x,y dimensions, we want y,x
        '''
        cmd = ' '.join(['ncpdq', '-a', 'y,x', '-O', myoutfile, myoutfile, '\n'])
        f.write(cmd)
        cmd = ' '.join(['ncpdq', '-a', 'y,x', '-O', 
                        os.path.join(odir, state_dir, outfile),
                        os.path.join(odir, state_dir, outfile), '\n'])
        f.write(cmd)
        '''


scripts = uniquify_list(scripts)
scripts_post = uniquify_list(scripts_post)
print '\n'.join([script for script in scripts])
print('\nwritten\n')
print '\n'.join([script for script in scripts_post])
print('\nwritten\n')

with open(os.path.join(odir, 'file_name_list.txt'), 'w') as f:
    for outfile in outfile_names:
        f.write(outfile+'\n')

if system in ['keeling']:
    with open('submit_batch.sh', 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write('cd '+batch_scripts_dir+'\n')
        for script in scripts:
            f.write('sbatch {}\n'.format(script))
    sub.call(['chmod', 'u+x', 'submit_batch.sh'])

