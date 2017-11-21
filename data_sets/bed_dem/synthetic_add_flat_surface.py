import numpy as np
from netCDF4 import Dataset
import sys

def add_flat_surface(infile=None, outfile=None, width=20000.0, loc='both', flat_elev=None):
    if not(loc in ['right', 'left', 'both']):
        sys.exit('loc is invalid')
    indata = Dataset(infile, 'r')
    outdata = Dataset(outfile, 'w')
    
    grid = abs(indata.variables['x'][1]-indata.variables['x'][0])
    ncols = len(indata.variables['x'][:])
    nrows = len(indata.variables['y'][:])

    z_in_var = indata.variables['topg']
    z_in = z_in_var[:]
    z_flat = np.zeros((z_in.shape[0], int(width/grid)))
    if flat_elev is None:
        flat_elev = z_in.min()
    z_flat[:] = flat_elev
    if loc == 'left':
        z_out = np.concatenate((z_flat, z_in), axis=1)
    elif loc == 'right':
        z_out = np.concatenate((z_in, z_flat), axis=1)
    elif loc == 'both':
        z_tmp = np.concatenate((z_in, z_flat), axis=1)
        z_out = np.concatenate((z_flat, z_tmp), axis=1)
    
    fill_value = -2.0e9
    #fill_value = z_in_var._FillValue
    z_out[np.where(z_out==fill_value)] = np.nan
    z_out = np.ma.array(z_out, mask=np.isnan(z_out))
    #z_out -= z_out.min()

    if loc =='left':
        outdata.createDimension('x', int(ncols+width/grid))
        x_var = outdata.createVariable('x', np.float64, ('x',))
        x_var[:] = np.arange(-width+grid/2.0,
                             indata.variables['x'][:].max()+grid/2.0,
                             grid)
    elif loc == 'right':
        outdata.createDimension('x', int(ncols+width/grid))
        x_var = outdata.createVariable('x', np.float64, ('x',))
        x_var[:] = np.arange(indata.variables['x'][:].min(),
                             indata.variables['x'][:].max()+width+grid/2.0,
                             grid)
    elif loc == 'both':
        outdata.createDimension('x', int(ncols+width*2/grid))
        x_var = outdata.createVariable('x', np.float64, ('x',))
        x_var[:] = np.arange(-width+grid/2.0,
                             indata.variables['x'][:].max()+width+grid/2.0,
                             grid)
    x_var.units = 'm'
    x_var.standard_name = 'projection_x_coordinate'

    outdata.createDimension('y', int(nrows))
    y_var = outdata.createVariable('y', np.float64, ('y',))
    y_var[:] = indata.variables['y'][:]
    y_var.units = 'm'
    y_var.standard_name = 'projection_y_coordinate'

    z_var = outdata.createVariable('topg', np.float64, ('y', 'x',),
                                   fill_value=fill_value)
    z_var[:] = z_out
    z_var.units = 'm'
    z_var.standard_name = 'bedrock_altitude'

    indata.close()
    outdata.close()
    
def add_ramp(infile=None, outfile=None, width=20000.0, loc='both', ramp_depth=200.0):
    if not(loc in ['right', 'left', 'both']):
        sys.exit('loc is invalid')
    indata = Dataset(infile, 'r')
    outdata = Dataset(outfile, 'w')
    
    grid = abs(indata.variables['x'][1]-indata.variables['x'][0])
    ncols = len(indata.variables['x'][:])
    nrows = len(indata.variables['y'][:])

    z_in_var = indata.variables['topg']
    z_in = z_in_var[:]

    z_ramp = np.zeros((z_in.shape[0], int(width/grid)))
    ni, nj = z_ramp.shape
    for j in range(nj):
        z_ramp[:, j] = -ramp_depth*(1.0-float(j)/float(nj))

    np.random.seed()
    z_ramp += np.random.rand(z_ramp.shape[0], z_ramp.shape[1])*5

    if loc == 'left':
        z_out = np.concatenate((z_ramp, z_in), axis=1)
    elif loc == 'right':
        z_out = np.concatenate((z_in, np.flip(z_ramp, 1)), axis=1)
    elif loc == 'both':
        z_tmp = np.concatenate((z_ramp, z_in), axis=1)
        z_out = np.concatenate((z_tmp, np.flip(z_ramp, 1)), axis=1)
    
    fill_value = -2.0e9
    #fill_value = z_in_var._FillValue
    z_out[np.where(z_out==fill_value)] = np.nan
    z_out = np.ma.array(z_out, mask=np.isnan(z_out))
    #z_out -= z_out.min()

    if loc =='left':
        outdata.createDimension('x', int(ncols+width/grid))
        x_var = outdata.createVariable('x', np.float64, ('x',))
        x_var[:] = np.arange(-width+grid/2.0,
                             indata.variables['x'][:].max()+grid/2.0,
                             grid)
    elif loc == 'right':
        outdata.createDimension('x', int(ncols+width/grid))
        x_var = outdata.createVariable('x', np.float64, ('x',))
        x_var[:] = np.arange(indata.variables['x'][:].min(),
                             indata.variables['x'][:].max()+width+grid/2.0,
                             grid)
    elif loc == 'both':
        outdata.createDimension('x', int(ncols+width*2/grid))
        x_var = outdata.createVariable('x', np.float64, ('x',))
        x_var[:] = np.arange(-width+grid/2.0,
                             indata.variables['x'][:].max()+width+grid/2.0,
                             grid)
    x_var.units = 'm'
    x_var.standard_name = 'projection_x_coordinate'

    outdata.createDimension('y', int(nrows))
    y_var = outdata.createVariable('y', np.float64, ('y',))
    y_var[:] = indata.variables['y'][:]
    y_var.units = 'm'
    y_var.standard_name = 'projection_y_coordinate'

    z_var = outdata.createVariable('topg', np.float64, ('y', 'x',),
                                   fill_value=fill_value)
    z_var[:] = z_out
    z_var.units = 'm'
    z_var.standard_name = 'bedrock_altitude'

    indata.close()
    outdata.close()

if __name__ == '__main__':
    add_flat_surface_left('test_dem.nc', 'test_dem_tmp.nc')
    add_flat_surface_right('test_dem_tmp.nc', 'test_dem_flat.nc')

