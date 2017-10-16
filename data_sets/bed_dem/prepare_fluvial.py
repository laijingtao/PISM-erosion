import subprocess as sub
import sys
sys.path.append('../../resources')
from analysis_tools import nc_regrid_dem
from cleanup_fluvial import cleanup
from synthetic_add_flat_surface import add_flat_surface

infile = 'fluvial.tif'
version = 1

fill_value = -100

cleanup(infile=infile, outfile='pism_Fluvial_v{}.nc'.format(version), 
        grid=50, zlim=[0., None], fill_value=fill_value)

for grid in [100, 250, 500, 1000]:
    print 'Grid resolution =', grid
    print 'Regriding...'
    nc_regrid_dem(infile='pism_Fluvial_v{}.nc'.format(version),
                  outfile='pism_Fluvial_{}m_v{}_original.nc'.format(grid, version),
                  grid=grid, fill_value=fill_value)
    print 'Adding flat surface...'
    add_flat_surface(infile='pism_Fluvial_{}m_v{}_original.nc'.format(grid, version),
                     outfile='pism_Fluvial_{}m_v{}.nc'.format(grid, version))
    cmd = ['rm', 'pism_Fluvial_{}m_v{}_original.nc'.format(grid, version)]
    print ' '.join(cmd)
    sub.call(cmd)
    

