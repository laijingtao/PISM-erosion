"""
resources
=========

Provides:
  - general resources such as grid constructors, calving, hydrology, etc.
    for the Olympic Peninsula

"""

from collections import OrderedDict
import os


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    if v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        print('Boolean value expected.')

def accepted_resolutions():

    return (50, 100, 200, 250, 500, 1000, 2000, 5000)


def generate_prefix_str(pism_exec):
    '''
    Generate prefix string.

    Returns: string
    '''

    try:
        p = os.environ['PISM_PREFIX']  + pism_exec
    except:
        p  = pism_exec
    
    return p


def generate_domain(domain):
    '''
    Generate domain specific options

    Returns: string
    '''
    
    if domain.lower() in ('olympics', 'olympics_mtns', 'synthetic', 'fluvial'):
        pism_exec = 'pismr'
    else:
        print('Domain {} not recognized, exiting'.format(domain))
        import sys
        sys.exit(0)

    return pism_exec


def default_spatial_ts_vars():
    '''
    Returns a list of commonly-used extra vars
    '''
    
    exvars = ['air_temp_snapshot',
              'beta',
              'bmelt',
              'cell_area',
              'climatic_mass_balance',
              'dHdt',
              'diffusivity',
              'effective_air_temp',
              'effective_precipitation',
              'ice_mass',
              'ice_surface_temp',
              'mask',
              'lat',
              'lon',
              'saccum',
              'smelt',
              'srunoff',
              'taub_mag',
              'tauc',
              'taud_mag',
              'tempicethk_basal',
              'temppabase',
              'tempsurf',
              'thk',
              'topg',
              'usurf',
              'velbase_mag',
              'velsurf_mag']
    
    return exvars



def generate_spatial_ts(outfile, exvars, step, start=None, end=None, split=None, odir=None):
    '''
    Return dict to generate spatial time series

    Returns: OrderedDict
    '''

    # check if list or comma-separated string is given.
    try:
        exvars = ','.join(exvars)
    except:
        pass

    params_dict = OrderedDict()
    if split is True:
        outfile, ext = os.path.splitext(outfile)
        params_dict['extra_split'] = ''
    if odir is None:
        params_dict['extra_file'] = 'ex_' + outfile
    else:
        params_dict['extra_file'] = os.path.join(odir, 'ex_' + outfile)
    params_dict['extra_vars'] = exvars
        
    if step is None:
        step = 'yearly'

    if (start is not None and end is not None):
        times = '{start}:{step}:{end}'.format(start=start, step=step, end=end)
    else:
        times = step
        
    params_dict['extra_times'] = times
        
  
    return params_dict


def generate_scalar_ts(outfile, step, start=None, end=None, odir=None):
    '''
    Return dict to create scalar time series

    Returns: OrderedDict
    '''

    params_dict = OrderedDict()
    if odir is None:
        params_dict['ts_file'] = 'ts_' + outfile
    else:
        params_dict['ts_file'] = os.path.join(odir, 'ts_' + outfile)
    
    if step is None:
        step = 'yearly'

    if (start is not None and end is not None):
        times = '{start}:{step}:{end}'.format(start=start, step=step, end=end)
    else:
        times = step
    params_dict['ts_times'] = times

    return params_dict


def generate_snap_shots(outfile, times, odir=None):
    '''
    Return dict to generate snap shots

    Returns: OrderedDict
    '''

    params_dict = OrderedDict()
    if odir is None:
        params_dict['save_file'] = 'save_' + outfile.split('.nc')[0]
    else:
        params_dict['save_file'] = os.path.join(odir, 'save_' + outfile.split('.nc')[0])

    params_dict['save_times'] = ','.join(str(e) for e in times)
    params_dict['save_split'] = ''
    params_dict['save_force_output_times'] = ''

    return params_dict


def generate_grid_description(grid_resolution, accepted_resolutions, domain,
                              dem_file=None, restart=False):
    '''
    Generate grid description dict

    Returns: OrderedDict
    '''

    if domain.lower() in ('olympics'):
        mx_max = 4000
        my_max = 3600
    elif domain.lower() in ('olympics_mtns'):
        mx_max = 2400
        my_max = 2000
    elif domain.lower() in ['synthetic', 'fluvial']:
        if dem_file is None:
            import sys
            sys.exit('Must provide a dem file')
        from netCDF4 import Dataset
        dem_data = Dataset(dem_file, 'r')
        x_len = abs(dem_data.variables['x'][0]-dem_data.variables['x'][-1])\
                + abs(dem_data.variables['x'][0]-dem_data.variables['x'][1])
        y_len = abs(dem_data.variables['y'][0]-dem_data.variables['y'][-1])\
                + abs(dem_data.variables['y'][0]-dem_data.variables['y'][1])
        dem_data.close()
        # mx_max = 40km/50m
        # my_max = 20km/50m
        mx_max = int(x_len/50.)
        my_max = int(y_len/50.)
    else:
        print('domain {} not recongnized'.format(domain))
    resolution_max = 50
    

    try:
        grid_resolution in accepted_resolutions
        pass
    except:
        print('grid resolution {}m not recognized'.format(grid_resolution))

    grid_div = (grid_resolution / resolution_max)
              
    mx = mx_max / grid_div
    my = my_max / grid_div

    horizontal_grid = OrderedDict()
    horizontal_grid['Mx'] = mx
    horizontal_grid['My'] = my

    if grid_resolution < 200:
        skip_max = 500
        mz = 101
        mzb = 21
    elif (grid_resolution >= 200) and (grid_resolution <= 500):
        skip_max = 250
        mz = 51
        mzb = 11
    else:
        skip_max = 100
        mz = 26
        mzb = 6

    vertical_grid = OrderedDict()
    vertical_grid['Lz'] = 2000
    vertical_grid['Lbz'] = 2000
    vertical_grid['z_spacing'] = 'equal'
    vertical_grid['Mz'] = mz
    vertical_grid['Mbz'] = mzb

    grid_options = {}
    grid_options['skip'] = ''
    grid_options['skip_max'] = skip_max

    grid_dict = merge_dicts(horizontal_grid, vertical_grid, grid_options)

    if restart is True:
        return grid_options
    else:
        return grid_dict


def merge_dicts(*dict_args):
    '''
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.

    Returns: OrderedDict
    '''
    result = OrderedDict()
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def uniquify_list(seq, idfun=None):
    '''
    Remove duplicates from a list, order preserving.
    From http://www.peterbe.com/plog/uniqifiers-benchmark
    '''

    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen:
            continue
        seen[marker] = 1
        result.append(item)
    return result


def generate_stress_balance(stress_balance, additional_params_dict):
    '''
    Generate stress balance params

    Returns: OrderedDict
    '''

    accepted_stress_balances = ('sia', 'ssa+sia')

    if stress_balance not in accepted_stress_balances:
        print('{} not in {}'.format(stress_balance, accepted_stress_balances))
        print('available stress balance solvers are {}'.format(accepted_stress_balances))
        import sys
        sys.exit(0)

    params_dict = OrderedDict()
    params_dict['stress_balance'] = stress_balance
    #if stress_balance in ('ssa+sia'):
    if stress_balance in ['ssa+sia']:
        params_dict['options_left'] = ''
        # params_dict['cfbc'] = ''
        #params_dict['sia_flow_law'] = 'gpbld3'
        params_dict['pseudo_plastic'] = ''
        params_dict['bed_smoother_range'] = 50

    return merge_dicts(additional_params_dict, params_dict)


def generate_hydrology(hydro, **kwargs):
    '''
    Generate hydrology params

    Returns: OrderedDict
    '''
    
    params_dict = OrderedDict()
    if hydro in ('null'):
        params_dict['hydrology'] = 'null'
    elif hydro in ('diffuse'):
        params_dict['hydrology'] = 'null'
        params_dict['hydrology_null_diffuse_till_water'] = ''
    elif hydro in ('routing'):
        params_dict['hydrology'] = 'routing'
    elif hydro in ('distributed'):
        params_dict['hydrology'] = 'distributed'
    else:
        print('hydrology {} not recognized, exiting'.format(hydro))
        import sys
        sys.exit(0)

    return merge_dicts(params_dict, kwargs)


def generate_calving(calving, **kwargs):
    '''
    Generate calving params

    Returns: OrderedDict
    '''

    params_dict = OrderedDict()
    if calving in ('ocean_kill'):
        params_dict['calving'] = calving
    elif calving in ('eigen_calving', 'vonmises_calving'):
        params_dict['calving'] = '{},thickness_calving'.format(calving)
    elif calving in ('hybrid_calving'):
        params_dict['calving'] = 'eigen_calving,vonmises_calving,thickness_calving'
    elif calving in ('float_kill'):
        params_dict['calving'] = calving
        params_dict['float_kill_margin_only'] = ''
    else:
        print('calving {} not recognized, exiting'.format(calving))
        import sys
        sys.exit(0)

    return merge_dicts(params_dict, kwargs)


def generate_climate(climate, **kwargs):
    '''
    Generate climate params

    Returns: OrderedDict
    '''
    
    params_dict = OrderedDict()
    ice_density = 910.
    if climate in ('elev'):
        params_dict['surface'] = 'elevation'
        params_dict['ice_surface_temp'] = '0,0,-100,5000'
        params_dict['climatic_mass_balance'] = '-3.,3,0,800,2500'
    elif climate in ('present'):
        params_dict['atmosphere'] = 'yearly_cycle,lapse_rate'
        params_dict['surface.pdd.factor_ice'] = 4.59 / ice_density  # Shea et al (2009)
        params_dict['surface.pdd.factor_snow'] = 3.04 / ice_density  # Shea et al (2009)
        params_dict['surface.pdd.refreeze'] = 0
        if 'atmosphere_yearly_cycle_file' not in kwargs:
            params_dict['atmosphere_yearly_cycle_file'] = 'olympics_climate_1000m.nc'
        else:
            params_dict['atmosphere_yearly_cycle_file'] = kwargs['atmosphere_yearly_cycle_file']
        if 'temp_lapse_rate' not in kwargs:
            params_dict['temp_lapse_rate'] = 4.5
        else:
            params_dict['temp_lapse_rate'] = kwargs['temp_lapse_rate']
        if 'atmosphere_lapse_rate_file' not in kwargs:
            params_dict['atmosphere_lapse_rate_file'] = 'olympics_climate_1000m.nc'
        else:
            params_dict['atmosphere_lapse_rate_file'] = kwargs['atmosphere_lapse_rate_file']
        params_dict['surface'] = 'pdd'
    elif climate in ('paleo'):
        params_dict['atmosphere'] = 'yearly_cycle,lapse_rate,delta_T,frac_P'
        params_dict['surface.pdd.factor_ice'] = 4.59 / ice_density  # Shea et al (2009)
        params_dict['surface.pdd.factor_snow'] = 3.04 / ice_density  # Shea et al (2009)
        params_dict['surface.pdd.refreeze'] = 0
        if 'atmosphere_yearly_cycle_file' not in kwargs:
            params_dict['atmosphere_yearly_cycle_file'] = 'olympics_climate_1000m.nc'
        else:
            params_dict['atmosphere_yearly_cycle_file'] = kwargs['atmosphere_yearly_cycle_file']
        if 'temp_lapse_rate' not in kwargs:
            params_dict['temp_lapse_rate'] = 4.5
        else:
            params_dict['temp_lapse_rate'] = kwargs['temp_lapse_rate']
        if 'atmosphere_lapse_rate_file' not in kwargs:
            params_dict['atmosphere_lapse_rate_file'] = 'olympics_climate_1000m.nc'
        else:
            params_dict['atmosphere_lapse_rate_file'] = kwargs['atmosphere_lapse_rate_file']
        if 'atmosphere_delta_T_file' not in kwargs:
            params_dict['atmosphere_delta_T_file'] = 'paleo_modifier.nc'
        else:
            params_dict['atmosphere_delta_T_file'] = kwargs['atmosphere_delta_T_file']
        if 'atmosphere_delta_T_file' not in kwargs:
            params_dict['atmosphere_delta_T_file'] = 'paleo_modifier.nc'
        else:
            params_dict['atmosphere_delta_T_file'] = kwargs['atmosphere_delta_T_file']
        '''
        if 'atmosphere_paleo_precip_file' not in kwargs:
            params_dict['atmosphere_paleo_precip_file'] = 'paleo_modifier.nc'
        else:
            params_dict['atmosphere_paleo_precip_file'] = kwargs['atmosphere_paleo_precip_file']
        '''
        params_dict['surface'] = 'pdd'
    elif climate in ('constant'):
        params_dict['atmosphere'] = 'yearly_cycle,lapse_rate,delta_T,frac_P'
        params_dict['surface.pdd.factor_ice'] = 4.59 / ice_density  # Shea et al (2009)
        params_dict['surface.pdd.factor_snow'] = 3.04 / ice_density  # Shea et al (2009)
        params_dict['surface.pdd.refreeze'] = 0
        if 'atmosphere_yearly_cycle_file' not in kwargs:
            params_dict['atmosphere_yearly_cycle_file'] = 'constant_climate.nc'
        else:
            params_dict['atmosphere_yearly_cycle_file'] = kwargs['atmosphere_yearly_cycle_file']
        if 'temp_lapse_rate' not in kwargs:
            params_dict['temp_lapse_rate'] = 4.5
        else:
            params_dict['temp_lapse_rate'] = kwargs['temp_lapse_rate']
        if 'atmosphere_lapse_rate_file' not in kwargs:
            params_dict['atmosphere_lapse_rate_file'] = 'constant_climate.nc'
        else:
            params_dict['atmosphere_lapse_rate_file'] = kwargs['atmosphere_lapse_rate_file']
        params_dict['surface'] = 'pdd'
    else:
        print('climate {} not recognized, exiting'.format(climate))
        import sys
        sys.exit(0)
        
    return merge_dicts(params_dict, kwargs)

        
def generate_ocean(ocean, **kwargs):
    '''
    Generate ocean params

    Returns: OrderedDict
    '''

    params_dict = OrderedDict()
    if ocean in ('null'):
        pass
    elif ocean in ('const'):
        params_dict['ocean'] = 'constant'
    else:
        print('ocean {} not recognized, exiting'.format(ocean))
        import sys
        sys.exit(0)

    return merge_dicts(params_dict, kwargs)


def list_systems():

    '''
    Return a list of supported systems.
    '''
    
    list = ['debug',
            'chinook',
            'keeling']
    
    return list


def list_queues():

    '''
    Return a list of supported queues.
    '''
    
    list = ['debug',
            'devel',
            'gpu',
            'gpu_long',
            'normal',
            'long',
            'standard',
            'standard_16',
            't2standard',
            't2small']
    
    return list


def make_batch_header(system, cores, walltime, queue):
    '''
    Generate header file for different HPC system.

    Returns: String
    '''
    
    systems = {}
    mpido = 'mpiexec -n {cores}'.format(cores=cores)
    systems['debug'] = {'mpido' : mpido,
                        'submit': 'echo',
                        'job_id' : 'PBS_JOBID'}
    mpido = 'mpirun -np {cores} -machinefile ./nodes_$SLURM_JOBID'.format(cores=cores)                         
    systems['chinook'] = {'mpido' : mpido,
                          'submit' : 'sbatch',
                          'work_dir' : 'SLURM_SUBMIT_DIR',
                          'job_id' : 'SLURM_JOBID',
                          'queue' : {
                              't2standard' : 24,
                              't2small' : 24,
                              'debug' : 24}}
    mpido = 'mpirun -np {cores}'.format(cores=cores)
    systems['keeling'] = {'mpido' : mpido,
                          'submit' : 'sbatch',
                          'work_dir' : 'SLURM_SUBMIT_DIR',
                          'job_id' : 'SLURM_JOBID',
                          'queue' : {
                              'normal' : 12 }}

    assert system in systems.keys()
    if system not in 'debug':
        assert queue in systems[system]['queue'].keys()
        assert cores > 0

        ppn = systems[system]['queue'][queue]
        nodes = cores / ppn

    if system in ('debug'):

        header = ''
        
    elif system in ('chinook'):
        
        header = """#!/bin/sh
#SBATCH --partition={queue}
#SBATCH --ntasks={cores}
#SBATCH --tasks-per-node={ppn}
#SBATCH --time={walltime}
#SBATCH --mail-user=aaschwanden@alaska.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=pism.%j

module list

cd $SLURM_SUBMIT_DIR

# Generate a list of compute node hostnames reserved for this job,
# this ./nodes file is necessary for slurm to spawn mpi processes
# across multiple compute nodes
srun -l /bin/hostname | sort -n | awk \'{{print $2}}\' > ./nodes_$SLURM_JOBID

""".format(queue=queue, walltime=walltime, nodes=nodes, ppn=ppn, cores=cores)
    elif system in ('keeling'):
        
        header = """#!/bin/bash
#SBATCH -n {cores}
#SBATCH --mem-per-cpu=4096
#SBATCH --time={walltime}
#SBATCH --mail-user=jlai11@illinois.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --export=PATH,LD_LIBRARY_PATH

module list

cd $SLURM_SUBMIT_DIR

""".format(walltime=walltime, cores=cores)
    else:
        header = """#!/bin/bash
#PBS -q {queue}
#PBS -l walltime={walltime}
#PBS -l nodes={nodes}:ppn={ppn}
#PBS -j oe

module list

cd $PBS_O_WORKDIR

""".format(queue=queue, walltime=walltime, nodes=nodes, ppn=ppn, cores=cores)

    return header, systems[system]


def make_batch_post_header(system):

    if system in ('chinook'):
        header = """#!/bin/bash
#PBS -q t2small
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe

module list

cd $PBS_O_WORKDIR

"""
    else:
        header = """#!/bin/bash
set -x -e

"""
    return header
