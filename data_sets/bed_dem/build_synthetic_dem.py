#!/usr/bin/env python
# Copyright (C) 2017 Jingtao Lai
# Generate synthetic fluvial bed dem

import numpy as np
from netCDF4 import Dataset

from landlab.components.flow_routing.route_flow_dn_JL import FlowRouter
from landlab.components.stream_power.fastscape_stream_power_JL import FastscapeEroder
from landlab import RasterModelGrid

def build_dem(grid=1000.):
    x_max = 40000.
    y_max = 20000.
    z_max = 3000.

    nrows = int(y_max/grid)
    ncols = int(x_max/grid)
    mg = RasterModelGrid(nrows, ncols, grid)
    mg.add_zeros('node', 'topographic__elevation', units='m')
    z = mg.at_node['topographic__elevation']
    z += (np.max(mg.node_x)-mg.node_x)/x_max*z_max

    return mg

def write_dem(mg, outfile):
    outdata = Dataset(outfile, 'w')

    nrows, ncols = mg.shape
    outdata.createDimension('x', ncols)
    outdata.createDimension('y', nrows)
    x_var = outdata.createVariable('x', np.float64, ('x',))
    y_var = outdata.createVariable('y', np.float64, ('y',))
    x_var[:] = mg.node_x[0:ncols]+mg.dx/2.
    y_var[:] = mg.node_y[0:nrows]+mg.dx/2.

    topg_var = outdata.createVariable('topg', np.float64, ('y', 'x',))
    topg_var[:] = mg.at_node['topographic__elevation']

    outdata.close()

if __name__ == '__main__':
    mg = build_dem()
    write_dem(mg, 'test_dem.nc')
