#!/usr/bin/env python
# Copyright (C) 2017 Jingtao Lai
# Generate synthetic fluvial bed dem

import numpy as np
from netCDF4 import Dataset

from landlab.components.flow_routing.route_flow_dn_JL import FlowRouter
from landlab.components.stream_power.fastscape_stream_power_JL import FastscapeEroder
from landlab.components.diffusion.diffusion import LinearDiffuser
from landlab import RasterModelGrid

def build_dem(grid=1000.):
    x_max = 40000.
    y_max = 20000.
    z_max = 0.
    dx = grid
    
    uplift_rate = 0.005
    runtime = 3000000.
    dt = 500.
    nt = int(runtime/dt)
    K = 1e-5
    D = 0.1

    nrows = int(y_max/grid)
    ncols = int(x_max/grid)
    mg = RasterModelGrid(nrows+2, ncols+2, dx)
    mg.add_zeros('node', 'topographic__elevation', units='m')
    z = mg.at_node['topographic__elevation']
    z += np.random.rand(len(z))/0.005
    z += (np.max(mg.node_x)-mg.node_x)/x_max*z_max
   
    mg.set_closed_boundaries_at_grid_edges(False, True, True, True)

    fr = FlowRouter(mg)
    sp = FastscapeEroder(mg, K_sp=K)
    lin_diffuse = LinearDiffuser(mg, linear_diffusivity=D)

    for i in range(nt):
        mg.at_node['topographic__elevation'][mg.core_nodes] += uplift_rate*dt
        # ramp
        mg.at_node['topographic__elevation'] += \
                0.0*np.sqrt((float(mg.node_x.max())-mg.node_x)/x_max)*uplift_rate*dt
        # valley
        mg.at_node['topographic__elevation'] += \
                0.0*np.absolute(mg.node_y-mg.node_y.mean())/(y_max/2)*uplift_rate*dt
        mg = fr.route_flow(routing_flat=False)
        mg = sp.erode(mg, dt)
        #mg = lin_diffuse.diffuse(dt)
        print 'Building... [{}%]\r'.format(int((i+1)*100.0/nt)),

    return mg

def write_dem(mg, outfile):
    print '\nWriting...'
    outdata = Dataset(outfile, 'w')

    nrows, ncols = mg.shape
    nrows -= 2
    ncols -= 2
    outdata.createDimension('x', ncols)
    outdata.createDimension('y', nrows)
    x_var = outdata.createVariable('x', np.float64, ('x',))
    y_var = outdata.createVariable('y', np.float64, ('y',))
    x_var[:] = mg.node_x[1:ncols+1]-mg.dx/2.
    y_var[:] = mg.node_y[1:nrows+1]-mg.dx/2.

    topg_var = outdata.createVariable('topg', np.float64, ('y', 'x',))
    topg_var[:] = mg.at_node['topographic__elevation'][mg.core_nodes]

    outdata.close()

if __name__ == '__main__':
    mg = build_dem(grid=500)
    write_dem(mg, 'test_dem.nc')
