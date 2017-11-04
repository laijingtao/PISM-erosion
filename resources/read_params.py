#!/usr/bin/env python
# Copyright (C) 2017 J. Lai

import numpy as np
import sys

def read_params(params_file=None, *args, **kwargs):
    if params_file is None:
        sys.exit('Must provide a params_file!')
    
    with open(params_file, 'r') as f:
        file_line_list = f.readlines()
    file_line_list = [file_line.strip() for file_line in file_line_list]

    i = 0
    j = 1
    params = {}
    while i<len(file_line_list):
        key = file_line_list[i].replace(':', '')
        params[key] = file_line_list[j]
        i += 2
        j += 2
    
    for key in params:
        if is_num(params[key]):
            params[key] = str2num(params[key])
        elif is_list(params[key]):
            params[key] = str2list(params[key])

    try:
        param_key=kwargs['param']
    except:
        param_key=None

    if param_key is None:
        return params
    else:
        try:
            value = params[param_key]
            return value
        except:
            sys.exit('Cannot find {} in {}'.format(param_key, params_file))

def is_num(s):
    try:
        float(s)
        return True
    except:
        return False

def str2num(s):
    return float(s)

def is_list(s):
    if s[0] in ['[']:
        return True
    else:
        return False

def str2list(s):
    s = s.replace('[', '')
    s = s.replace(']', '')
    num_list = [float(num_s) for num_s in s.split(',')]
    return num_list
