#!/usr/bin/env python
# Copyright (C) 2017 Jingtao Lai

import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
import itertools
import os
import sys
sys.path.append('../resources/')
from resources import *

pism_data_dir = '/mnt/d/jlai11/pism-olympics/run/2017_07_calib/'

domain = 'olympics'
grid = '1000'
start = -125000
end = 0
stress_balance = 'ssa+sia'
climate = 'paleo'

sia_e_values = [1.0, 3.0]
q_values = [0.33, 0.5]
phi_values = [20, 30]
combinations = list(itertools.product(sia_e_values,
                                      q_values,
                                      phi_values))

for n, combination in enumerate(combinations):
    sia_e, q, phi = combination

    name_options = OrderedDict()
    name_options['sb'] = stress_balance
    name_options['sia_e'] = sia_e
    name_options['q'] = q
    name_options['phi'] = phi
    experiment =  '_'.join([climate, '_'.join(['_'.join([k, str(v)]) for k, v
        in name_options.items()])])
    file_name =\
    '{domain}_g{grid}m_{experiment}_{start}_{end}a.nc'.format(domain=domain.lower(),
                                                              grid = grid,
                                                              experiment = experiment,
                                                              start = int(start),
                                                              end = int(end))
    print(file_name)
