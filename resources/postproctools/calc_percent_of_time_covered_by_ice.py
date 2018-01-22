from .common import *

def calc_percent_of_time_covered_by_ice(infile=None, outfile=None):
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
   
    var_list = ['velbase_mag', 'erosion_1', 'erosion_2']
    flag = True
    for var in var_list:
        if var in indata.variables.keys():
            flag = False
            break
    if flag:
        sys.exit('No usable variable')

    percent = np.zeros(indata.variables[var][:][0].shape)
    for i in range(len(time)):
        #print i
        var_slice = indata.variables[var][:][i]
        #percent[np.where(np.logical_not(var_slice.mask))] += 1.
        percent[np.where(var_slice.data>0.)] += 1.
    percent = percent/float(len(time))*100.

    try:
        fill_value = indata.variables[var]._FillValue
    except:
        fill_value = FILL_VALUE
    percent[np.where(percent<=0.0)] = fill_value
    try:
        percent_var = outdata.createVariable('percent_of_time_covered_by_ice',
                                             np.float64,
                                             ('y', 'x',),
                                             fill_value=fill_value)
    except:
        percent_var = outdata.variables['percent_of_time_covered_by_ice']
    percent_var[:] = percent

    indata.close()
    if not overwrite:
        outdata.close()