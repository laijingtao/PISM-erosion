#!/usr/bin/env python
# Copyright (C) 2017 Jingtao Lai

import numpy as np
from scipy import signal

def xcf(a=None, b=None, dt=None):
    if a is None or b is None:
        sys.exit('Missing array')
    if len(a)!=len(b):
        sys.exit('Two arrays must have the same dimension')

    n = len(a)
    corr = signal.correlate(a-a.mean(), b-b.mean(), mode='full')/(n*np.std(a)*np.std(b))
    lag = np.arange(len(corr))
    lag = lag-n+1
    
    if dt is not None:
        lag = lag*dt

    return corr, lag
