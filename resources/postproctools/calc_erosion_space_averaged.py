from .common import *

def calc_erosion_space_averaged(infile=None, outfile=None):
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

    time = indata.variables['time'][:]
    grid = abs(indata.variables['x'][0]-indata.variables['x'][1])
    col = len(indata.variables['x'])
    row = len(indata.variables['y'])

    for erosion_name in ['erosion_1', 'erosion_2']:
        erosion = indata.variables[erosion_name][:]
        total_erosion = np.array(
            [erosion_slice.sum() for erosion_slice in erosion])*grid*grid
        erosion_space_averaged = total_erosion/(grid*grid*row*col)
        try:
            erosion_space_averaged_var = outdata.createVariable(
                erosion_name+'_space_averaged', np.float64, ('time',))
        except:
            erosion_space_averaged_var = outdata.variables[erosion_name+'_space_averaged']
        erosion_space_averaged_var[:] = erosion_space_averaged
        erosion_space_averaged_var.units = 'm year-1'

    indata.close()
    if not overwrite:
        outdata.close()
        