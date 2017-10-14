import numpy as np
from netCDF4 import Dataset
import subprocess as sub

def cleanup(infile=None, outfile=None, grid=100, *args, **kwargs):
    cmd = ['/usr/bin/gdalwarp', '-of', 'netCDF', '-overwrite',
           '-tr', str(grid), str(grid), infile, 'tmp_fluvial.nc']
    print ' '.join(cmd)
    sub.call(cmd)

    indata = Dataset('tmp_fluvial.nc', 'r')
    outdata = Dataset(outfile, 'w')

    topg_in = indata.variables['Band1'][:].copy()
    topg_in = topg_in.astype(float)
    topg_in[np.where(topg_in<=0)] = np.nan
    
    print 'Cleaning...'
    x_min = 1e9
    x_max = -1
    y_min = 1e9
    y_max = -1

    ny, nx = topg_in.shape
    for i in range(ny):
        for j in range(nx):
            if np.isnan(topg_in[i, j]):
                continue
            else:
                if i<y_min:
                    y_min = i
                if i>y_max:
                    y_max = i
                if j<x_min:
                    x_min = j
                if j>x_max:
                    x_max = j

    try:
        cut_width = kwargs['cut_width']
    except:
        cut_width = 1500.
    x_min += int(cut_width/grid)

    nrows = y_max-y_min+1
    ncols = x_max-x_min+1

    outdata.createDimension('x', ncols)
    outdata.createDimension('y', nrows)
    x_var = outdata.createVariable('x', np.float64, ('x',))
    y_var = outdata.createVariable('y', np.float64, ('y',))
    x_var[:] = np.arange(ncols)*grid+grid/2.0
    y_var[:] = np.arange(nrows)*grid+grid/2.0
    x_var.units = 'm'
    y_var.units = 'm'

    topg_var = outdata.createVariable('topg', np.float64, ('y', 'x',), fill_value=-100)
    topg_out = topg_in[y_min:y_max+1, x_min:x_max+1]
    topg_out = np.ma.array(topg_out, mask=np.isnan(topg_out))
    try:
        zmin, zmax = kwargs['zlim']
    except:
        zmin = None
        zmax = None
    if zmin is None:
        zmin = topg_out.min()
    if zmax is None:
        zmax = topg_out.max()
    topg_out = (topg_out-topg_out.min())/(topg_out.max()-topg_out.min())*(zmax-zmin)+zmin
    topg_var[:] = topg_out
    topg_var.units = 'm'

    indata.close()
    outdata.close()

    cmd = ['rm', 'tmp_fluvial.nc']
    print ' '.join(cmd)
    sub.call(cmd)
    
if __name__=='__main__':
    cleanup(infile='fluvial.tif', outfile='pism_Fluvial_100m_v1.nc', zlim=[0., 3000.])
