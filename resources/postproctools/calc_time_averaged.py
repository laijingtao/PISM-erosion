from .common import *

def calc_time_averaged(infile=None, outfile=None, var_name=None, *args, **kwargs):
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
        var[np.where(var == FILL_VALUE)] = 0.0
    var_time_averaged = np.zeros(var[0].shape)
    for i in range(start_index, end_index+1):
        tmp_var_slice = var[i]
        if i == start_index:
            time_step = (time[i+1]-time[i])/2
        elif i == end_index:
            time_step = (time[i]-time[i-1])/2
        else:
            time_step = (time[i+1]-time[i-1])/2
        var_time_averaged = var_time_averaged+tmp_var_slice*time_step
    var_time_averaged = var_time_averaged/(time[end_index]-time[start_index])
    var_time_averaged[np.where(var_time_averaged<=0.)] = FILL_VALUE
    var_time_averaged = np.ma.masked_values(var_time_averaged,
                                            FILL_VALUE)
    try:
        var_time_averaged_var = outdata.createVariable(
            '{}_time_averaged'.format(var_name), np.float64, ('y','x',),
            fill_value=FILL_VALUE)
    except:
        var_time_averaged_var = outdata.variables['{}_time_averaged'.format(var_name)]
        var_time_averaged_var._FillValue = FILL_VALUE
    var_time_averaged_var[:] = var_time_averaged
    try:
        var_time_averaged_var.units = indata.variables[var_name].units
    except:
        pass

    indata.close()
    if not overwrite:
        outdata.close()
