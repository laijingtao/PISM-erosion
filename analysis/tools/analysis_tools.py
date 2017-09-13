#!/usr/bin/env python
# Copyright (C) 2017 Jingtao Lai

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import itertools
try:
    import subprocess32 as sub
except:
    import subprocess as sub
from collections import OrderedDict
from argparse import ArgumentParser
#sys.path.append('../../resources/')
#from resources import *


def generate_file_names():
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

    file_name_list = []
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
                                                                  grid=grid,
                                                                  experiment=experiment,
                                                                  start=int(start),
                                                                  end=int(end))
        file_name_list.append(file_name)
        
    return file_name_list

def calc_erosion(infile=None, outfile=None):
    if infile is None:
        sys.exit('Must provide an input file!')
    if outfile is None:
        outfile = infile
    
    cmd = ['ncap2', '-O', '-s',
           'erosion_1=1e-4*velbase_mag;erosion_2=2.7e-7*velbase_mag^2.02;',
           infile, outfile]
    #print sub.list2cmdline(cmd)
    sub.call(cmd)

