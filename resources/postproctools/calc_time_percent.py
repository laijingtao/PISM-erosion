from .common import *

def calc_time_percent(infile=None, outfile=None, var_name=None, *args, **kwargs):
    """Percent% of the time the value of var is lower than the return value
    """
    if var_name is None:
        sys.exit('Must provide a var name!')
    try:
        percent = kwargs['percent']
    except:
        sys.exit('Must provide a percent!')
    
    if infile is None:
        sys.exit('Must provide an input file!')
    overwrite = False
    if outfile is None or outfile==infile:
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

    var_sorted = indata.variables[var_name][:].copy()
    var_sorted.sort(axis=0)

    ncols, nrows = var_sorted[0].shape
    var_percent_value = np.zeros((ncols, nrows))
    for i in range(ncols):
        for j in range(nrows):
            try:
                n_time = var_sorted[:, i, j].count()
            except:
                n_time = len(var_sorted[:, i, j])
            if n_time==0:
                var_percent_value[i, j] = FILL_VALUE
            else:
                index = int(percent/100.0*n_time-1)
                var_percent_value[i, j] = var_sorted[:, i, j][index]
    var_percent_value[np.where(var_percent_value<=0)] = FILL_VALUE
    var_percent_value = np.ma.masked_values(var_percent_value, FILL_VALUE)

    try:
        var_percent = outdata.createVariable(
            '{}_time_{}_percent'.format(var_name, percent), np.float64, ('y','x',),
            fill_value=FILL_VALUE)
    except:
        var_percent = outdata.variables['{}_time_{}_percent'.format(var_name, percent)]
    var_percent[:] = var_percent_value
    var_percent.units = indata.variables[var_name].units

    indata.close()
    if not overwrite:
        outdata.close()