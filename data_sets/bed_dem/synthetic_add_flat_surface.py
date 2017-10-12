import numpy as np
from netCDF4 import Dataset

def add_flat_surface(infile=None, outfile=None, width=20000.0):
    indata = Dataset(infile, 'r')
    outdata = Dataset(outfile, 'w')
    
    grid = abs(indata.variables['x'][1]-indata.variables['x'][0])
    ncols = len(indata.variables['x'][:])
    nrows = len(indata.variables['y'][:])

    z_in_var = indata.variables['topg']
    z_in = z_in_var[:]
    z_flat = np.zeros((z_in.shape[0], int(width/grid)))
    z_flat[:] = z_in.min()
    z_out = np.concatenate((z_in, z_flat), axis=1)
    z_out[np.where(z_out==z_in_var._FillValue)] = np.nan
    z_out = np.ma.array(z_out, mask=np.isnan(z_out))
    z_out -= z_out.min()

    outdata.createDimension('x', int(ncols+width/grid))
    outdata.createDimension('y', int(nrows))
    x_var = outdata.createVariable('x', np.float64, ('x',))
    y_var = outdata.createVariable('y', np.float64, ('y',))
    z_var = outdata.createVariable('topg', np.float64, ('y', 'x',),
                                   fill_value=z_in_var._FillValue)

    x_var[:] = np.arange(indata.variables['x'][:].min(),
                         indata.variables['x'][:].max()+width+grid/2.0,
                         grid)
    x_var.units = 'm'
    y_var[:] = indata.variables['y'][:]
    y_var.units = 'm'
    z_var[:] = z_out
    z_var.units = 'm'

    indata.close()
    outdata.close()


if __name__ == '__main__':
    add_flat_surface('test_dem.nc', 'test_dem_flat.nc')

