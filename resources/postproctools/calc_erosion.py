from .common import *

def calc_erosion(infile=None, outfile=None):
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

    for erosion_name in ['erosion_1', 'erosion_2']:
        if erosion_name == 'erosion_1':
            erosion = 1e-4*indata.variables['velbase_mag'][:]
        if erosion_name == 'erosion_2':
            erosion = 2.7e-7*np.power(indata.variables['velbase_mag'][:], 2.02)
        erosion[np.where(erosion <= 0)] = FILL_VALUE
        erosion = np.ma.masked_values(erosion, FILL_VALUE)
        try:
            erosion_var = outdata.createVariable(
                erosion_name, np.float64, ('time', 'y', 'x',),
                fill_value=FILL_VALUE)
        except:
            erosion_var = outdata.variables[erosion_name]
            erosion_var._FillValue = FILL_VALUE
        erosion_var[:] = erosion
        erosion_var.units = 'm year-1'

    indata.close()
    if not overwrite:
        outdata.close()
        