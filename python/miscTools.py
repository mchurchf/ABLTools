#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 15:39:03 2017

@author: pdoubraw
"""
import numpy as np
#==============================================================================
# 
#==============================================================================
def parse(string):
    """
    Convert strings to float.
    """
    try:
        return float(string)
    except:
        return np.nan