#!/usr/bin/env python
# Copyright (C) 2017 Jingtao Lai
# Generate synthetic fluvial bed dem

import numpy as np
from netCDF4 import Dataset

from landlab.components.flow_routing.route_flow_dn_JL import FlowRouter
from landlab.components.stream_power.fastscape_stream_power_JL import FastscapeEroder
from landlab.components.diffusion.diffusion import LinearDiffuser
from landlab import RasterModelGrid

def build_dem(grid=1000., z_limit=None, *args, **kwargs):
    x_max = 100000.
    y_max = 100000.
    z_max = 0.
    dx = grid
    
    try:
        uplift_rate = kwargs['uplift_rate']
    except:
        uplift_rate = 0.005
    try:
        runtime = kwargs['runtime']
    except:
        runtime = 5000000.
    try:
        dt = kwargs['dt']
    except:
        dt = 10000.
    nt = int(runtime/dt)
    try:
        K = kwargs['K']
    except:
        K = 1e-5
    try:
        D = kwargs['D']
    except:
        D = 0.01

    nrows = int(y_max/grid)
    ncols = int(x_max/grid)
    mg = RasterModelGrid(nrows+2, ncols+2, dx)
    mg.add_zeros('node', 'topographic__elevation', units='m')
    z = mg.at_node['topographic__elevation']
    z += np.random.rand(len(z))/0.1
    z += (np.max(mg.node_x)-mg.node_x)/x_max*z_max
   
    mg.set_closed_boundaries_at_grid_edges(False, False, False, False)

    fr = FlowRouter(mg)
    sp = FastscapeEroder(mg, K_sp=K, threshold_sp=50.0*K)
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
        mg = lin_diffuse.diffuse(dt)
        print 'Building... [{}%]\r'.format(int((i+1)*100.0/nt)),

    if z_limit is not None:
        z = mg.at_node['topographic__elevation']
        z[np.where(z>z_limit)] = z_limit

    return mg

def extract_basin(mg):
    z = mg.at_node['topographic__elevation']
    node_stack = mg.at_node['flow__upstream_node_order']
    receiver = mg.at_node['flow__receiver_node']

    '''
    if mg.node_x[np.argmax(mg.at_node['drainage_area'])]>mg.node_x.max()/2:
        right_outlet = True
    else:
        right_outlet = False
    '''
    right_outlet = False
    left_edge, = np.where(mg.node_x==0)
    outlet = left_edge[np.argmax(mg.at_node['drainage_area'][left_edge])]
    isbasin = np.zeros(len(z), dtype=bool)
    isbasin[outlet] = True
    for node in node_stack:
        if isbasin[receiver[node]]:
            isbasin[node] = True
    z[np.where(np.logical_not(isbasin))] = np.nan

    x_min = mg.node_x[np.where(isbasin)].min()
    x_max = mg.node_x[np.where(isbasin)].max()
    if right_outlet:
        x_max -= 2000
    else:
        x_min += 2000
    y_min = mg.node_y[np.where(isbasin)].min()
    y_max = mg.node_y[np.where(isbasin)].max()
    new_domain, = np.where(np.logical_and(
        np.logical_and(mg.node_x>=x_min, mg.node_x<=x_max),
        np.logical_and(mg.node_y>=y_min, mg.node_y<=y_max)))

    ncols = int((x_max-x_min)/mg.dx+1)
    nrows = int((y_max-y_min)/mg.dx+1)
    dx = mg.dx
    basin = RasterModelGrid(nrows+2, ncols+2, dx)
    basin.add_zeros('node', 'topographic__elevation', units='m')
    basin_z = basin.at_node['topographic__elevation']
    basin_z[basin.core_nodes] = z[new_domain]
    if right_outlet:
        tmp = basin_z.reshape(nrows+2, ncols+2)
        tmp = np.fliplr(tmp)
        basin_z = tmp.reshape((nrows+2)*(ncols+2))
        basin.at_node['topographic__elevation'] = basin_z
    
    return basin

def write_dem(mg, outfile, zmin=None, zmax=None):
    print '\nWriting...'
    outdata = Dataset(outfile, 'w')

    nrows, ncols = mg.shape
    nrows -= 2
    ncols -= 2
    outdata.createDimension('x', ncols)
    outdata.createDimension('y', nrows)
    x_var = outdata.createVariable('x', np.float64, ('x',))
    y_var = outdata.createVariable('y', np.float64, ('y',))
    x_var[:] = mg.node_x[np.where(mg.node_y==0)][1:ncols+1]-mg.dx/2.
    x_var.units = 'm'
    y_var[:] = mg.node_y[np.where(mg.node_x==0)][1:nrows+1]-mg.dx/2.
    y_var.units = 'm'

    topg_var = outdata.createVariable('topg', np.float64, ('y', 'x',), fill_value=-100.0)
    z = mg.at_node['topographic__elevation'][mg.core_nodes]
    #z[np.where(z>zmax)] = np.nan
    z = np.ma.array(z, mask=np.isnan(z))
    z = z.reshape(nrows, ncols)
    if zmin is None:
        zmin = z.min()
    if zmax is None:
        zmax = z.max()
    z = (z-z.min())/(z.max()-z.min())*(zmax-zmin)+zmin
    topg_var[:] = z
    topg_var.units = 'm'

    outdata.close()

if __name__ == '__main__':
    mg = build_dem(grid=1000, zmax=None)
    mg = extract_basin(mg)
    write_dem(mg, 'test_dem.nc', zmin=0, zmax=3000)
